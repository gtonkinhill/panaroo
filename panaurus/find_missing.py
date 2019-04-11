import networkx as nx
from skbio import Sequence, DNA, Protein
import io, sys
from skbio.io import read
from skbio.metadata import Interval, IntervalMetadata
from skbio.alignment import local_pairwise_align_ssw, StripedSmithWaterman
from collections import defaultdict
import numpy as np
from Bio.Seq import translate, reverse_complement
from cdhit import align_dna_cdhit


def find_missing(G, gff_file_handles, dna_seq_file, prot_seq_file, temp_dir, n_cpu):

    # find all the cycles shorter than cycle_threshold
    print("defining basis...")
    complete_basis = set()
    all_pairs = dict(nx.all_pairs_shortest_path(G, cutoff=2))
    for source in all_pairs:
        for sink in all_pairs[source]:
            if len(all_pairs[source][sink]) == 3:
                mid_size = G.node[all_pairs[source][sink][1]]['size']
                if (G.node[source]['size'] > mid_size) or (G.node[sink]['size']
                                                           > mid_size):
                    path = all_pairs[source][sink]
                    complete_basis.add((min(path[0], path[2]), path[1],
                                        max(path[0], path[2])))

    # For each cycle check if it looks like somethings missing
    print("identify missing nodes...")
    search_count=0
    search_lists = defaultdict(list)
    for b in complete_basis:
        # identify which genomes are missing the smallest supported node
        surrounding_members = set(G.node[b[0]]['members']) & set(
            G.node[b[2]]['members'])
        for member in surrounding_members:
            if member not in G.node[b[1]]['members']:
                # we have a possible search target. First check if theres
                # another path between the nodes
                paths = list(nx.all_shortest_paths(G, b[0], b[2]))
                good_target = True
                for path in paths:
                    if len(path) < 3: continue
                    if member in (set(G.node[path[0]]['members']) & set(
                            G.node[path[1]]['members']) & set(
                                G.node[path[2]]['members'])):
                        good_target = False
                if good_target:
                    search_lists[member].append(b)
                    search_count+=1

    print("num searches", search_count)
    print("search for missing nodes...")
    # find position of each search by genome to save reading in the same gff3
    # more than once
    n_found = 0
    hit_dict = {}
    for member in search_lists:
        neighbour_id_list = []
        search_sequence_list = []
        missing = []
        for b in search_lists[member]:
            neighbour_ids = []
            for n in [0, 2]:
                for id in G.node[b[n]]['seqIDs']:
                    if id.split("_")[0] == member:
                        neighbour_ids.append(id)
            neighbour_id_list.append(neighbour_ids)
            search_sequence_list.append(
                max(G.node[b[1]]["dna"].split(";"), key=len))
            missing.append(b)

        hit_dict[member] = search_seq_gff(
            gff_handle=gff_file_handles[int(member)],
            neighbour_id_list=neighbour_id_list,
            search_sequence_list=search_sequence_list,
            missing=missing,
            temp_dir=temp_dir,
            n_cpu=n_cpu)

    print("update output...")
    with open(dna_seq_file, 'a') as dna_out:
        with open(prot_seq_file, 'a') as prot_out:
            for member in search_lists:
                for b, hit in zip(search_lists[member], hit_dict[member]):
                    if hit == "": continue
                    G.node[b[1]]['members'] += [member]
                    G.node[b[1]]['size'] += 1
                    G.node[b[1]]['dna'] = ";".join(
                        set(G.node[b[1]]['dna'].split(";") + [hit]))
                    dna_out.write(">" + str(member) + "_refound_" +
                                  str(n_found) + "\n" + hit + "\n")
                    # hit_protein = translate_to_match(
                    #     hit, max(G.node[b[1]]["protein"].split(";"), key=len))
                    # G.node[b[1]]['protein'] = ";".join(
                    #     set(G.node[b[1]]['protein'].split(";") +
                    #         [hit_protein]))
                    # prot_out.write(">" + str(member) + "_refound_" +
                    #                str(n_found) + "\n" + hit_protein + "\n")
                    G.node[b[1]]['seqIDs'] += [
                        str(member) + "_refound_" + str(n_found)
                    ]
                    n_found += 1

    return G


def search_seq_gff(gff_handle,
                   neighbour_id_list,
                   search_sequence_list,
                   missing,
                   temp_dir,
                   prop_match=0.2,
                   pairwise_id_thresh=0.95,
                   n_cpu=1):

    # reset file handle to the beginning
    gff_handle.seek(0)
    gff3_str = gff_handle.read().split("##FASTA\n")

    if len(gff3_str) != 2:
        raise NameError("File does not appear to be in GFF3 format!")

    contig_records = defaultdict(dict)
    contig_names = []
    for seq in read(io.StringIO(gff3_str[1]), format='fasta'):
        if getattr(seq, 'metadata')['id'] in contig_records:
            raise NameError("Duplicate contig names!")
        contig_records[getattr(seq, 'metadata')['id']]['seq'] = str(seq)
        contig_names.append(getattr(seq, 'metadata')['id'])

    gff = read(io.StringIO(gff3_str[0]), format='gff3')
    for ann in gff:
        if ann[0] not in contig_records:
            raise NameError("Mismatch in GFF file!")
        annotations = list(ann[1].query([(0, 1000000000)]))
        annotations = [
            a for a in annotations if getattr(a, 'metadata')['type'] == 'CDS'
        ]
        if len(annotations) > 0:
            contig_records[ann[0]]['annotations'] = annotations

    # TODO: for now skip entries with no annotationas we skip them reading in
    # the GFF3. May want to adjust this in the future
    for record in contig_records:
        if 'annotations' not in contig_records[record]:
            contig_names.remove(record)

    hits = []
    for neighbour_ids, seq, mis in zip(neighbour_id_list, search_sequence_list,
                                       missing):
        found_dna = ""
        gene_locations = [(int(tid.split("_")[1]), int(tid.split("_")[2]))
                          for tid in neighbour_ids]
        # print(neighbour_ids)
        contigA = contig_names[gene_locations[0][0]]
        contigB = contig_names[gene_locations[1][0]]
        gene_num_A = gene_locations[0][1]
        gene_num_B = gene_locations[1][1]
        metaA = contig_records[contigA]['annotations'][gene_num_A]
        metaB = contig_records[contigB]['annotations'][gene_num_B]

        # determine search area in contigs
        found_dna = ""
        if contigA == contigB:
            # the flanking genes are on the same contig, search inbetween
            l_bound = min(
                getattr(metaA, 'bounds')[0][1],
                getattr(metaB, 'bounds')[0][1])
            r_bound = max(
                getattr(metaA, 'bounds')[0][0],
                getattr(metaB, 'bounds')[0][0])
            search_sequence = contig_records[contigA]['seq'][l_bound:r_bound]
            found_dna = search_dna(seq, search_sequence, prop_match,
                                   pairwise_id_thresh, temp_dir, n_cpu)
        else:
            # if it looks like the genes are near the terminal ends search here
            max_search_length = 2*len(seq)
            if gene_num_A < 20:
                # we're at the start of contigA
                search_sequence = contig_records[contigA]['seq'][:min(
                    getattr(metaA, 'bounds')[0])]
                if len(search_sequence)< max_search_length:
                    found_dna = search_dna(seq, search_sequence, prop_match,
                                           pairwise_id_thresh, temp_dir, n_cpu)
            if (found_dna == "") and (len(
                    contig_records[contigA]['annotations']) - gene_num_A < 20):
                # we're at the  end of contigA
                search_sequence = contig_records[contigA]['seq'][max(
                    getattr(metaA, 'bounds')[0]):]
                if len(search_sequence)< max_search_length:
                    found_dna = search_dna(seq, search_sequence, prop_match,
                                           pairwise_id_thresh, temp_dir, n_cpu)
            if (found_dna == "") and (gene_num_B < 20):
                # we're at the start of contigB
                search_sequence = contig_records[contigB]['seq'][:min(
                    getattr(metaB, 'bounds')[0])]
                if len(search_sequence)< max_search_length:
                    found_dna = search_dna(seq, search_sequence, prop_match,
                                           pairwise_id_thresh, temp_dir, n_cpu)
            if (found_dna == "") and (len(
                    contig_records[contigB]['annotations']) - gene_num_B < 20):
                # we're at the  end of contigB
                search_sequence = contig_records[contigB]['seq'][max(
                    getattr(metaB, 'bounds')[0]):]
                if len(search_sequence)< max_search_length:
                    found_dna = search_dna(seq, search_sequence, prop_match,
                                           pairwise_id_thresh, temp_dir, n_cpu)

        # if found_dna=="":
        #     print(neighbour_ids)
        #     print(len(contig_records[contigA]['annotations']),
        #         len(contig_records[contigB]['annotations']))
        #     print(len(contig_records[contigA]['seq']),
        #         len(contig_records[contigB]['seq']))
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


def search_dna(seq, search_sequence, prop_match, pairwise_id_thresh, temp_dir, n_cpu):

    found_dna= align_dna_cdhit(query=search_sequence,
        target=seq,
        id=0.95,
        temp_dir=temp_dir,
        n_cpu=n_cpu,
        aS=0.3)

    return found_dna


def translate_to_match(hit, target_prot):

    # translate in all 6 frames splitting on unknown
    dna_seqs = [hit, reverse_complement(hit)]

    proteins = [
        translate(s[i:].ljust(len(s[i:]) + (3 - len(s[i:]) % 3), 'N'))
        for i in range(3) for s in dna_seqs
    ]

    query = StripedSmithWaterman(target_prot)

    alignments = []
    for target_sequence in proteins:
        alignment = query(target_sequence)
        alignments.append((target_sequence, alignment.optimal_alignment_score))

    prot = max(alignments, key=lambda x: x[1])

    return prot[0]


blosum50 = \
    {
        '*': {'*': 1, 'A': -5, 'C': -5, 'B': -5, 'E': -5, 'D': -5, 'G': -5,
              'F': -5, 'I': -5, 'H': -5, 'K': -5, 'M': -5, 'L': -5,
              'N': -5, 'Q': -5, 'P': -5, 'S': -5, 'R': -5, 'T': -5,
              'W': -5, 'V': -5, 'Y': -5, 'X': -5, 'Z': -5},
        'A': {'*': -5, 'A': 5, 'C': -1, 'B': -2, 'E': -1, 'D': -2, 'G': 0,
              'F': -3, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -2,
              'N': -1, 'Q': -1, 'P': -1, 'S': 1, 'R': -2, 'T': 0, 'W': -3,
              'V': 0, 'Y': -2, 'X': -1, 'Z': -1},
        'C': {'*': -5, 'A': -1, 'C': 13, 'B': -3, 'E': -3, 'D': -4,
              'G': -3, 'F': -2, 'I': -2, 'H': -3, 'K': -3, 'M': -2,
              'L': -2, 'N': -2, 'Q': -3, 'P': -4, 'S': -1, 'R': -4,
              'T': -1, 'W': -5, 'V': -1, 'Y': -3, 'X': -1, 'Z': -3},
        'B': {'*': -5, 'A': -2, 'C': -3, 'B': 6, 'E': 1, 'D': 6, 'G': -1,
              'F': -4, 'I': -4, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 5,
              'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -5, 'V': -3,
              'Y': -3, 'X': -1, 'Z': 1},
        'E': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 6, 'D': 2, 'G': -3,
              'F': -3, 'I': -4, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0,
              'Q': 2, 'P': -1, 'S': -1, 'R': 0, 'T': -1, 'W': -3, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 5},
        'D': {'*': -5, 'A': -2, 'C': -4, 'B': 6, 'E': 2, 'D': 8, 'G': -1,
              'F': -5, 'I': -4, 'H': -1, 'K': -1, 'M': -4, 'L': -4, 'N': 2,
              'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -5, 'V': -4,
              'Y': -3, 'X': -1, 'Z': 1},
        'G': {'*': -5, 'A': 0, 'C': -3, 'B': -1, 'E': -3, 'D': -1, 'G': 8,
              'F': -4, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0,
              'Q': -2, 'P': -2, 'S': 0, 'R': -3, 'T': -2, 'W': -3, 'V': -4,
              'Y': -3, 'X': -1, 'Z': -2},
        'F': {'*': -5, 'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -5,
              'G': -4, 'F': 8, 'I': 0, 'H': -1, 'K': -4, 'M': 0, 'L': 1,
              'N': -4, 'Q': -4, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 1,
              'V': -1, 'Y': 4, 'X': -1, 'Z': -4},
        'I': {'*': -5, 'A': -1, 'C': -2, 'B': -4, 'E': -4, 'D': -4,
              'G': -4, 'F': 0, 'I': 5, 'H': -4, 'K': -3, 'M': 2, 'L': 2,
              'N': -3, 'Q': -3, 'P': -3, 'S': -3, 'R': -4, 'T': -1,
              'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -3},
        'H': {'*': -5, 'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2,
              'F': -1, 'I': -4, 'H': 10, 'K': 0, 'M': -1, 'L': -3, 'N': 1,
              'Q': 1, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -3, 'V': -4,
              'Y': 2, 'X': -1, 'Z': 0},
        'K': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2,
              'F': -4, 'I': -3, 'H': 0, 'K': 6, 'M': -2, 'L': -3, 'N': 0,
              'Q': 2, 'P': -1, 'S': 0, 'R': 3, 'T': -1, 'W': -3, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 1},
        'M': {'*': -5, 'A': -1, 'C': -2, 'B': -3, 'E': -2, 'D': -4,
              'G': -3, 'F': 0, 'I': 2, 'H': -1, 'K': -2, 'M': 7, 'L': 3,
              'N': -2, 'Q': 0, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -1,
              'V': 1, 'Y': 0, 'X': -1, 'Z': -1},
        'L': {'*': -5, 'A': -2, 'C': -2, 'B': -4, 'E': -3, 'D': -4,
              'G': -4, 'F': 1, 'I': 2, 'H': -3, 'K': -3, 'M': 3, 'L': 5,
              'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -1,
              'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3},
        'N': {'*': -5, 'A': -1, 'C': -2, 'B': 5, 'E': 0, 'D': 2, 'G': 0,
              'F': -4, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -4, 'N': 7,
              'Q': 0, 'P': -2, 'S': 1, 'R': -1, 'T': 0, 'W': -4, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 0},
        'Q': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2,
              'F': -4, 'I': -3, 'H': 1, 'K': 2, 'M': 0, 'L': -2, 'N': 0,
              'Q': 7, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -1, 'V': -3,
              'Y': -1, 'X': -1, 'Z': 4},
        'P': {'*': -5, 'A': -1, 'C': -4, 'B': -2, 'E': -1, 'D': -1,
              'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -3,
              'L': -4, 'N': -2, 'Q': -1, 'P': 10, 'S': -1, 'R': -3,
              'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': -1},
        'S': {'*': -5, 'A': 1, 'C': -1, 'B': 0, 'E': -1, 'D': 0, 'G': 0,
              'F': -3, 'I': -3, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1,
              'Q': 0, 'P': -1, 'S': 5, 'R': -1, 'T': 2, 'W': -4, 'V': -2,
              'Y': -2, 'X': -1, 'Z': 0},
        'R': {'*': -5, 'A': -2, 'C': -4, 'B': -1, 'E': 0, 'D': -2, 'G': -3,
              'F': -3, 'I': -4, 'H': 0, 'K': 3, 'M': -2, 'L': -3, 'N': -1,
              'Q': 1, 'P': -3, 'S': -1, 'R': 7, 'T': -1, 'W': -3, 'V': -3,
              'Y': -1, 'X': -1, 'Z': 0},
        'T': {'*': -5, 'A': 0, 'C': -1, 'B': 0, 'E': -1, 'D': -1, 'G': -2,
              'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0,
              'Q': -1, 'P': -1, 'S': 2, 'R': -1, 'T': 5, 'W': -3, 'V': 0,
              'Y': -2, 'X': -1, 'Z': -1},
        'W': {'*': -5, 'A': -3, 'C': -5, 'B': -5, 'E': -3, 'D': -5,
              'G': -3, 'F': 1, 'I': -3, 'H': -3, 'K': -3, 'M': -1, 'L': -2,
              'N': -4, 'Q': -1, 'P': -4, 'S': -4, 'R': -3, 'T': -3,
              'W': 15, 'V': -3, 'Y': 2, 'X': -1, 'Z': -2},
        'V': {'*': -5, 'A': 0, 'C': -1, 'B': -3, 'E': -3, 'D': -4, 'G': -4,
              'F': -1, 'I': 4, 'H': -4, 'K': -3, 'M': 1, 'L': 1, 'N': -3,
              'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 5,
              'Y': -1, 'X': -1, 'Z': -3},
        'Y': {'*': -5, 'A': -2, 'C': -3, 'B': -3, 'E': -2, 'D': -3,
              'G': -3, 'F': 4, 'I': -1, 'H': 2, 'K': -2, 'M': 0, 'L': -1,
              'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -1, 'T': -2, 'W': 2,
              'V': -1, 'Y': 8, 'X': -1, 'Z': -2},
        'X': {'*': -5, 'A': -1, 'C': -1, 'B': -1, 'E': -1, 'D': -1,
              'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1,
              'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1,
              'T': -1, 'W': -1, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1},
        'Z': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 5, 'D': 1, 'G': -2,
              'F': -4, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0,
              'Q': 4, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -2, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 5}}
