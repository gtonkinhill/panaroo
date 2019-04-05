import networkx as nx
from skbio import Sequence, DNA
import io, sys
from skbio.io import read
from skbio.metadata import Interval, IntervalMetadata
from skbio.alignment import local_pairwise_align_ssw
from collections import defaultdict

def find_missing(G, gff_file_handles):

    # atm only check small cycles to greatly simplify things
    cycle_threshold=3

    # find all the cycles shorter than cycle_threshold
    complete_basis = []
    for c in nx.connected_components(G):
        sub_G = G.subgraph(c)
        basis = nx.cycle_basis(sub_G, list(sub_G.nodes())[0])
        complete_basis += [b for b in basis if len(b) == cycle_threshold]

    # For each cycle check if it looks like somethings missing
    search_lists = defaultdict(list)
    for b in complete_basis:
        # find node with smallest support
        smallest_support = 99999999
        for n in b:
            if G.node[n]['size'] < smallest_support:
                smallest = n
                smallest_support = G.node[n]['size']
        b.remove(smallest)

        # identify which genomes are missing the smallest supported node
        surrounding_members = set(G.node[b[0]]['members']) & set(G.node[b[1]]['members'])
        for member in surrounding_members:
            if member not in G.node[smallest]['members']:
                search_lists[member].append((smallest, b))

    # find position of each search by genome to save reading in the same gff3
    # more than once
    for member in search_lists:
        neighbour_id_list = []
        search_sequence_list = []
        for mis in search_lists[member]:
            neighbour_ids = []
            for n in mis[1]:
                for id in G.node[n]['seqIDs']:
                    if id.split("_")[0]==member:
                        neighbour_ids.append(id)
            neighbour_id_list.append(neighbour_ids)
            search_sequence_list.append(max(G.node[mis[0]]["dna"].split(";"),
                key=len))

        hit = search_seq_gff(gff_handle=gff_file_handles[int(member)],
            neighbour_id_list=neighbour_id_list,
            search_sequence_list=search_sequence_list)

    return G


def search_seq_gff(gff_handle, neighbour_id_list, search_sequence_list):

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
    for neighbour_ids, seq in zip(neighbour_id_list, search_sequence_list):
        gene_locations = [(int(tid.split("_")[1]),
                            int(tid.split("_")[2])) for tid in neighbour_ids]
        contigA = contig_names[gene_locations[0][0]]
        contigB = contig_names[gene_locations[1][0]]
        gene_num_A = gene_locations[0][1]
        gene_num_B = gene_locations[1][1]

        # determine search area in contigs
        if contigA==contigB:
            # the flanking genes are on the same contig

            metaA = contig_records[contigA]['annotations'][gene_num_A]
            metaB = contig_records[contigB]['annotations'][gene_num_B]

            bounds = list(getattr(metaA, 'bounds')[0]) + list(getattr(metaB, 'bounds')[0])
            search_sequence = contig_records[contigA]['seq'][min(bounds):max(bounds)]

            # print(bounds)
            # print(min(bounds),max(bounds))
            # print(metaA)
            # print(metaB)
            # print(search_sequence)
            # print(seq)

            msa = local_pairwise_align_ssw(DNA(seq), DNA(search_sequence))
            best_hit_bounds = msa[2][1]
            if abs(best_hit_bounds[1]-best_hit_bounds[0]) > (0.5*len(seq)):
                # we've found a decent match, return it
                found_dna = search_sequence[best_hit_bounds[0]:best_hit_bounds[1]]
            hits.append(found_dna)


        # else if True:
        #     # it looks like the genes are near the terminal ends search here
        #     sg
        # else:
        #     # we play it safe and don't search for the gene as we can't be sure
        #     # of direction to search



    # print(gff_handle.name)
    # print(neighbour_ids)
    # # print(search_sequence)
    #
    # flanking_gene_num = [int(id.split("_")[-1]) for id in neighbour_ids]

    return