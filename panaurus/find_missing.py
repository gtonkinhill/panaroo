import networkx as nx
from skbio import Sequence, DNA
import io, sys
from skbio.io import read
from skbio.metadata import Interval, IntervalMetadata
from skbio.alignment import local_pairwise_align_ssw
from collections import defaultdict

def find_missing(G, gff_file_handles):

    # find all the cycles shorter than cycle_threshold
    complete_basis = set()
    all_pairs = dict(nx.all_pairs_shortest_path(G, cutoff=2))
    for source in all_pairs:
        for sink in all_pairs[source]:
            if len(all_pairs[source][sink])==3:
                mid_size = G.node[all_pairs[source][sink][1]]['size']
                if (G.node[source]['size'] > mid_size) or (
                    G.node[sink]['size'] > mid_size):
                    path = all_pairs[source][sink]
                    complete_basis.add((min(path[0], path[2]), path[1], max(path[0], path[2])))

    # For each cycle check if it looks like somethings missing
    search_lists = defaultdict(list)
    for b in complete_basis:
        # identify which genomes are missing the smallest supported node
        surrounding_members = set(G.node[b[0]]['members']) & set(G.node[b[2]]['members'])
        for member in surrounding_members:
            if member not in G.node[b[1]]['members']:
                # we have a possible search target. First check if theres another
                # path between the nodes
                paths = list(nx.all_shortest_paths(G, b[0], b[2]))
                good_target=True
                for path in paths:
                    if len(path)<3: continue
                    if member in (set(G.node[path[0]]['members']) & set(G.node[path[1]]['members']) & set(G.node[path[2]]['members'])):
                        good_target=False
                if good_target:
                    search_lists[member].append(b)

        print(G.node[b[0]]['size'])
        print(G.node[b[1]]['size'])
        print(G.node[b[2]]['size'])

    # find position of each search by genome to save reading in the same gff3
    # more than once
    for member in search_lists:
        neighbour_id_list = []
        search_sequence_list = []
        missing = []
        for b in search_lists[member]:
            neighbour_ids = []
            for n in [0,2]:
                for id in G.node[b[n]]['seqIDs']:
                    if id.split("_")[0]==member:
                        neighbour_ids.append(id)
            neighbour_id_list.append(neighbour_ids)
            search_sequence_list.append(max(G.node[b[1]]["dna"].split(";"),
                key=len))
            missing.append(b)

        hits = search_seq_gff(gff_handle=gff_file_handles[int(member)],
            neighbour_id_list=neighbour_id_list,
            search_sequence_list=search_sequence_list,
            missing=missing)

        for b, hit in zip(search_lists[member], hits):
            if hit=="": continue
            G.node[b[1]]['members'] += [member]
            G.node[b[1]]['size'] += 1
            G.node[b[1]]['dna']=";".join(
                set(G.node[b[1]]['dna'].split(";") +
                    [hit]))

    return G


def search_seq_gff(gff_handle, neighbour_id_list, search_sequence_list,
    missing,
    prop_match=0.3):

    # reset file handle to the beginning
    gff_handle.seek(0)
    gff3_str = gff_handle.read().split("##FASTA\n")

    if len(gff3_str)!=2:
        raise NameError("File does not appear to be in GFF3 format!")

    contig_records = defaultdict(dict)
    contig_names = []
    for seq in read(io.StringIO(gff3_str[1]), format='fasta'):
        if getattr(seq, 'metadata')['id'] in contig_records:
            raise NameError("Duplicate contig names!")
        contig_records[getattr(seq, 'metadata')['id']]['seq'] = str(seq)
        contig_names.append(getattr(seq, 'metadata')['id'])

    gff = read(io.StringIO(gff3_str[0]), format='gff3')
    for ann in  gff:
        if ann[0] not in contig_records:
            raise NameError("Mismatch in GFF file!")
        contig_records[ann[0]]['annotations'] = list(ann[1].query([(0,1000000000)]))

    # TODO: for now skip entries with no annotationas we skip them reading in the GFF3. May want to adjust this in the future
    for record in contig_records:
        if 'annotations' not in contig_records[record]:
            contig_names.remove(record)

    hits = []
    for neighbour_ids, seq, mis in zip(neighbour_id_list, search_sequence_list, missing):
        found_dna =  ""
        gene_locations = [(int(tid.split("_")[1]),
                            int(tid.split("_")[2])) for tid in neighbour_ids]
        contigA = contig_names[gene_locations[0][0]]
        contigB = contig_names[gene_locations[1][0]]
        gene_num_A = gene_locations[0][1]
        gene_num_B = gene_locations[1][1]
        metaA = contig_records[contigA]['annotations'][gene_num_A]
        metaB = contig_records[contigB]['annotations'][gene_num_B]

        # determine search area in contigs
        found_dna = ""
        if contigA==contigB:
            # the flanking genes are on the same contig, search inbetween
            l_bound = min(getattr(metaA, 'bounds')[0][1], getattr(metaB, 'bounds')[0][1])
            r_bound = max(getattr(metaA, 'bounds')[0][0], getattr(metaB, 'bounds')[0][0])
            search_sequence = contig_records[contigA]['seq'][l_bound:r_bound]
            found_dna = search_dna(seq, search_sequence, prop_match)
        else:
            # if it looks like the genes are near the terminal ends search here
            if gene_num_A<10:
                # we're at the start of contigA
                search_sequence = contig_records[contigA]['seq'][:min(getattr(metaA, 'bounds')[0])]
                found_dna = search_dna(seq, search_sequence, prop_match)
            if (found_dna=="") and (len(contig_records[contigA]['annotations'])-gene_num_A < 10):
                # we're at the  end of contigA
                search_sequence = contig_records[contigA]['seq'][max(getattr(metaA, 'bounds')[0]):]
                found_dna = search_dna(seq, search_sequence, prop_match)
            if (found_dna=="") and (gene_num_B<10):
                # we're at the start of contigB
                search_sequence = contig_records[contigB]['seq'][:min(getattr(metaB, 'bounds')[0])]
                found_dna = search_dna(seq, search_sequence, prop_match)
            if (found_dna=="") and (len(contig_records[contigB]['annotations'])-gene_num_B<10):
                # we're at the  end of contigB
                search_sequence = contig_records[contigB]['seq'][max(getattr(metaB, 'bounds')[0]):]
                found_dna = search_dna(seq, search_sequence, prop_match)
            if found_dna!="":
                print("SUCCESSS!!")

        # if found_dna=="":
        #     print(neighbour_ids)
        #     print(len(contig_records[contigA]['annotations']), len(contig_records[contigB]['annotations']))
        #     print(len(contig_records[contigA]['seq']), len(contig_records[contigB]['seq']))
        #     print(getattr(metaA, 'bounds')[0], getattr(metaB, 'bounds')[0])
        #     # print(mis)
        #     print(len(search_sequence), len(seq))
        #     # # print(bounds)
        #     # # print(min(bounds), max(bounds))
        #     # print(metaA)
        #     # print(metaB)
        #     # print(search_sequence)
        #     print(seq)

        # add results
        hits.append(found_dna)

    # print(hits)

    return hits

def search_dna(seq, search_sequence, prop_match):
    found_dna = ""
    msa = max([local_pairwise_align_ssw(DNA(seq), DNA(search_sequence)),
            local_pairwise_align_ssw(DNA(seq).reverse_complement(), DNA(search_sequence))],
        key=lambda x: x[1])
    best_hit_bounds = msa[2][1]
    if (abs(best_hit_bounds[1]-best_hit_bounds[0]) > (prop_match*len(seq))):
        # we've found a decent match, return it
        found_dna = search_sequence[best_hit_bounds[0]:best_hit_bounds[1]]

    return found_dna
