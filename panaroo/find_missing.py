import networkx as nx
import io, sys
from collections import defaultdict
import numpy as np
from Bio.Seq import translate, reverse_complement, Seq
from Bio import SeqIO
from panaroo.cdhit import align_dna_cdhit
from joblib import Parallel, delayed
import os
import gffutils as gff
from io import StringIO
from .merge_nodes import delete_node, remove_member_from_node
import pyopa 


def get_all_paths(G, length=3):
    all_paths = []
    for node in G.nodes():
        neighs = [v for u, v in nx.bfs_edges(G, source=node, depth_limit=length)]
        for neigh in neighs:
            for p in nx.all_simple_paths(G, source=node, target=neigh, cutoff=2):
                if len(p)==length:
                    all_paths.append(p)
    return(all_paths)


def find_missing(G, gff_file_handles, dna_seq_file, prot_seq_file, temp_dir,
                 n_cpu, remove_by_consensus=False):

    # find all the cycles shorter than cycle_threshold
    print("defining basis...")
    complete_basis = set()
    all_paths = get_all_paths(G, 3)

    for path in all_paths:
        mid_size = G.node[path[1]]['size']
        if (G.node[path[0]]['size'] > mid_size) or (
            G.node[path[2]]['size'] > mid_size):
            complete_basis.add(
                        (min(path[0], path[2]), path[1], max(path[0],
                                                             path[2])))
        
    # For each cycle check if it looks like somethings missing
    print("identify missing nodes...")
    print("len(complete_basis):", len(complete_basis))
    search_count = 0
    search_lists = defaultdict(list)
    # seen_pairs = set()
    for b in complete_basis:
        # identify which genomes are missing the smallest supported node
        surrounding_members = set(G.node[b[0]]['members']) & set(
            G.node[b[2]]['members'])
        for member in surrounding_members:
            if member not in G.node[b[1]]['members']:
                # if (member, b[1]) in seen_pairs: continue
                search_lists[member].append(b)
                # seen_pairs.add((member, b[1]))
                search_count += 1

    print("num searches", search_count)
    print("search for missing nodes...")
    # find position of each search by genome to save reading in the same gff3
    # more than once
    n_found = 0

    print("setting up sample searches")
    neighbour_dict = {}
    search_seq_dict = {}
    missing_dict = {}
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
                G.node[b[1]]["dna"].split(";")[0])
            missing.append(b)

        neighbour_dict[member] = neighbour_id_list
        search_seq_dict[member] = search_sequence_list
        missing_dict[member] = missing

    hit_list_init = Parallel(n_jobs=n_cpu)(delayed(search_seq_gff)(
        member, gff_file_handles[int(member)], neighbour_dict[member],
        search_seq_dict[member], temp_dir, 0.2)
                                      for member in search_lists)

    # only take the best hit for each node/member pair
    hit_list = []
    for member, init_hits in hit_list_init:
        node_mems = defaultdict(list)
        hits = []
        searches = []
        for b, hit in zip(search_lists[member], init_hits):
            node_mems[b[1]].append([hit, b])
        for b1 in node_mems:
            best_hit = max(node_mems[b1], key = lambda x: len(x[0]))
            hits.append(best_hit[0])
            searches.append(best_hit[1])
        search_lists[member] = searches
        hit_list.append([member, hits])


    print("translating found hits...")
    trans_list = []
    for member, hits in hit_list:
        trans_list.append(
            Parallel(n_jobs=n_cpu)(
                delayed(translate_to_match)(
                    hit, G.node[b[1]]["protein"].split(";")[0])
                for b, hit in zip(search_lists[member], hits)))

    additions_by_node = defaultdict(list)
    for mem_hits, trans in zip(hit_list, trans_list):
        member = mem_hits[0]
        hits = mem_hits[1]
        for b, hit, hit_protein in zip(search_lists[member], hits,
                                        trans):
            if hit == "": continue
            additions_by_node[b[1]].append((member, hit, hit_protein))

    bad_nodes = []

    # Check if there's more than one path between nodes i.e we found the same bit of DNA twice
    path_mem_pairs = defaultdict(set)
    for path in complete_basis:
        for mem in G.node[path[1]]['members']:
            path_mem_pairs[(path[0],path[2],mem)].add(path[1])
        for mem, hit, hit_protein in additions_by_node[path[1]]:
            path_mem_pairs[(path[0],path[2],mem)].add(path[1])
    
    for pmp in path_mem_pairs:
        if len(path_mem_pairs[pmp]) > 1:
            
            neigh_mem_count = 0
            for nA in path_mem_pairs[pmp]:
                memfound=0
                for nB in G.neighbors(nA):
                    if pmp[2] in G[nA][nB]['members']:
                        memfound=1
                    break
                neigh_mem_count += memfound
            if neigh_mem_count>1: continue

            # we have two options for a pair for the same member -> delete one
            node_max = -1
            best_node = -1
            for node in path_mem_pairs[pmp]:
                if node in G.nodes():
                    if G.node[node]['size']>node_max:
                        best_node = node
                        node_max = G.node[node]['size']
            # remove member from nodes that arent the best
            for node in path_mem_pairs[pmp]:
                if node==best_node: continue
                G = remove_member_from_node(G, node, pmp[2])
                additions_by_node[node] = [hit for hit in additions_by_node[node] if hit[0]!=pmp[2]]
    # clean up graph by removing empty nodes
    for node in G.nodes():
        if G.node[node]['size']<1:
            bad_nodes.append(node)
    for node in bad_nodes:
        delete_node(G, node)
        del additions_by_node[node]

    # if requested remove nodes that have more refound than in original
    # that is the consensus appears to be it wasn't a good gene most of the
    # time
    if remove_by_consensus:
        for node in additions_by_node:
            if len(additions_by_node[node])>G.node[node]['size']:
                bad_nodes.append(node)
    for node in bad_nodes:
        if node in G.nodes():
            delete_node(G, node)

    # search nodes at the ends with degree 1
    neighbour_dict = defaultdict(list)
    search_seq_dict = defaultdict(list)
    end_search_lists = defaultdict(list)
    for node in G.nodes():
        if G.degree[node]==1:
            neigh = list(G.neighbors(node))[0]
            for member in G.node[neigh]['members']:
                if member not in G.node[node]['members']:
                    for id in G.node[neigh]['seqIDs']:
                        if id.split("_")[0] == member:
                            hitid = id
                    neighbour_dict[member].append([hitid])
                    search_seq_dict[member].append(G.node[node]["dna"].split(";")[0])
                    end_search_lists[member].append(node)
    
    end_mems = neighbour_dict.keys()
    end_hit_list = Parallel(n_jobs=n_cpu)(delayed(search_seq_gff)(
        member, gff_file_handles[int(member)], neighbour_dict[member],
        search_seq_dict[member], temp_dir, 0.2)
                                        for member in end_mems)
        
    for member, hits in end_hit_list:
        trans = Parallel(n_jobs=n_cpu)(
                delayed(translate_to_match)(
                    hit, G.node[node]["protein"].split(";")[0])
                    for node, hit in zip(end_search_lists[member], hits))
        for node, hit, trans_hit in zip(end_search_lists[member], hits, trans):
            if hit=="": continue
            additions_by_node[node].append((member, hit, trans_hit))
    
    print("update output...")
    with open(dna_seq_file, 'a') as dna_out:
        with open(prot_seq_file, 'a') as prot_out:
            for node in additions_by_node:
                if node in bad_nodes: continue
                for member, hit, hit_protein in additions_by_node[node]:
                    G.node[node]['members'] += [member]
                    G.node[node]['size'] += 1
                    G.node[node]['dna'] = ";".join(
                        set(G.node[node]['dna'].split(";") + [hit]))
                    dna_out.write(">" + str(member) + "_refound_" +
                                    str(n_found) + "\n" + hit + "\n")
                    G.node[node]['protein'] = ";".join(
                        set(G.node[node]['protein'].split(";") +
                            [hit_protein]))
                    prot_out.write(">" + str(member) + "_refound_" +
                                    str(n_found) + "\n" + hit_protein + "\n")
                    G.node[node]['seqIDs'] += [
                        str(member) + "_refound_" + str(n_found)
                    ]
                    n_found += 1

    # add edges
    for mem_hits in hit_list:
        member = mem_hits[0]
        hits = mem_hits[1]
        for b, hit in zip(search_lists[member], hits):
            if hit == "": continue
            # remove member from old edge
            is_triplet = True
            for n in b:
                if not G.has_node(n):
                    is_triplet=False
                elif member not in G.node[n]['members']:
                    is_triplet=False
            if not is_triplet: continue
            for edge in [(b[0],b[1]), (b[1], b[2]), (b[0],b[2])]:
                if not G.has_edge(edge[0], edge[1]):
                    is_triplet=False
            if not is_triplet: continue
            if member not in G[b[0]][b[1]]['members']:
                if member not in G[b[1]][b[2]]['members']:
                    if member in G[b[0]][b[2]]['members']:
                        G[b[0]][b[2]]['members'].remove(member)
                        G[b[0]][b[2]]['weight'] -= 1
                        if G[b[0]][b[2]]['weight']==0:
                            G.remove_edge(b[0],b[2])
                        for edge in [(b[0],b[1]), (b[1], b[2])]:
                            G[edge[0]][edge[1]]['members'].append(member)
                            G[edge[0]][edge[1]]['weight'] += 1

    return G


def search_seq_gff(member,
                   gff_handle,
                   neighbour_id_list,
                   search_sequence_list,
                #    missing,
                   temp_dir,
                   prop_match=0.2,
                   pairwise_id_thresh=0.95,
                   n_cpu=1):

    # get matrix needed for alingment.We use a protein modified to be suitable for DNA 
    # are dealing with high identities it does okay on DNA
    env = pyopa.AlignmentEnvironment()
    env.gap_ext = np.array(-4)
    env.gap_open = np.array(-12)
    m = np.full((26,26), -3.0)
    np.fill_diagonal(m , 5.0)
    env.float64_matrix = m
    env.float64_matrix = env.float64_matrix.astype("float64")
    env.create_scaled_matrices()

    # reset file handle to the beginning
    gff_handle.seek(0)
    split = gff_handle.read().split("##FASTA\n")

    if len(split) != 2:
        raise NameError("File does not appear to be in GFF3 format!")

    contig_records = defaultdict(dict)
    contig_names = []
    with StringIO(split[1]) as temp_fasta:
        for record in SeqIO.parse(temp_fasta, 'fasta'):
            if record.id in contig_records:
                raise NameError("Duplicate contig names!")
            contig_records[record.id]['seq'] = str(record.seq)
            contig_names.append(record.id)

    parsed_gff = gff.create_db("\n".join(
        [l for l in split[0].splitlines() if '##sequence-region' not in l]),
                               dbfn=":memory:",
                               force=True,
                               keep_order=True,
                               from_string=True)

    for entry in parsed_gff.all_features(featuretype=()):
        if "CDS" not in entry.featuretype: continue
        if entry.seqid not in contig_records:
            raise NameError("Mismatch in GFF file!")
        if 'annotations' not in contig_records[entry.seqid]:
            contig_records[entry.seqid]['annotations'] = [entry]
        else:
            contig_records[entry.seqid]['annotations'].append(entry)

    # TODO: for now skip entries with no annotationas we skip them reading in
    # the GFF3. May want to adjust this in the future
    for record in contig_records:
        if 'annotations' not in contig_records[record]:
            contig_names.remove(record)

    hits = []
    for neighbour_ids, seq in zip(neighbour_id_list, search_sequence_list,
                                       ):
        found_dna = ""
        
        # printmore=False
        # for tid in neighbour_ids:
        #     if tid in ["0_5_127","1_1_639","2_1_562","3_0_182","2_1_561","1_1_640","1_1_641","2_1_559","2_1_560","0_5_128","3_0_181"]:
        #         printmore=True

        gene_locations = [(int(tid.split("_")[1]), int(tid.split("_")[2]))
                          for tid in neighbour_ids]
        # print(neighbour_ids)
        contigA = contig_names[gene_locations[0][0]]
        if len(gene_locations)>1:
            contigB = contig_names[gene_locations[1][0]]
        gene_num_A = gene_locations[0][1]
        if len(gene_locations)>1:
            gene_num_B = gene_locations[1][1]

        metaA = contig_records[contigA]['annotations'][gene_num_A]
        if len(gene_locations)>1:
            metaB = contig_records[contigB]['annotations'][gene_num_B]

        # if printmore:
        #     print(neighbour_ids)
        #     print(contigA, contigB)
        #     print(gene_num_A, gene_num_B)
        #     print(prop_match)
        #     print(pairwise_id_thresh)

        # determine search area in contigs
        found_dna = ""
        if (len(gene_locations)>1) and ((contigA == contigB) and (abs(gene_num_A - gene_num_B) <= 2)):
            # the flanking genes are on the same contig, search inbetween
            l_bound = min(metaA.end, metaB.end, metaA.start, metaB.start)
            r_bound = max(metaA.start, metaB.start, metaA.end, metaB.end)
            search_sequence = contig_records[contigA]['seq'][l_bound:r_bound]
            found_dna = search_dna(seq, search_sequence, prop_match,
                                   pairwise_id_thresh, env)#temp_dir, n_cpu)#
        else:
            # if it looks like the genes are near the terminal ends search here
            max_search_length = 10 * len(seq)
            if gene_num_A < 20:
                # we're at the start of contigA
                search_sequence = contig_records[contigA]['seq'][:min(
                    metaA.start, metaA.end)]
                if len(search_sequence) < max_search_length:
                    found_dna = search_dna(seq, search_sequence, prop_match,
                                           pairwise_id_thresh, env)#temp_dir, n_cpu)#
            if (found_dna == "") and (len(
                    contig_records[contigA]['annotations']) - gene_num_A < 20):
                # we're at the  end of contigA
                search_sequence = contig_records[contigA]['seq'][max(
                    metaA.start, metaA.end):]
                if len(search_sequence) < max_search_length:
                    found_dna = search_dna(seq, search_sequence, prop_match,
                                           pairwise_id_thresh, env)#temp_dir, n_cpu)#
            if len(gene_locations)>1:
                if (found_dna == "") and (gene_num_B < 20):
                    # we're at the start of contigB
                    search_sequence = contig_records[contigB]['seq'][:min(
                        metaB.start, metaB.end)]
                    if len(search_sequence) < max_search_length:
                        found_dna = search_dna(seq, search_sequence, prop_match,
                                            pairwise_id_thresh, env)#temp_dir, n_cpu)#
                if (found_dna == "") and (len(
                        contig_records[contigB]['annotations']) - gene_num_B < 20):
                    # we're at the  end of contigB
                    search_sequence = contig_records[contigB]['seq'][max(
                        metaB.start, metaB.end):]
                    if len(search_sequence) < max_search_length:
                        found_dna = search_dna(seq, search_sequence, prop_match,
                                            pairwise_id_thresh, env)#temp_dir, n_cpu)#

        # add results
        hits.append(found_dna)

        # if printmore:
        #     print(search_sequence)
        #     print(seq)
        #     print(found_dna)

    # print(hits)

    return (member, hits)


def search_dna(seq, search_sequence, prop_match, pairwise_id_thresh, env):#temp_dir, n_cpu):

    # found_dna = align_dna_cdhit(query=search_sequence,
    #                             target=seq,
    #                             id=pairwise_id_thresh,
    #                             temp_dir=temp_dir,
    #                             n_cpu=n_cpu,
    #                             use_local=True,
    #                             mask=True,
    #                             aS=prop_match)

    # # if nothing found and lots of Ns try searching without masking Ns
    # if found_dna=="":
    #     found_dna = align_dna_cdhit(query=search_sequence,
    #                             target=seq,
    #                             id=pairwise_id_thresh,
    #                             temp_dir=temp_dir,
    #                             n_cpu=n_cpu,
    #                             use_local=True,
    #                             mask=False,
    #                             aS=prop_match)
    gap = np.array("-")
    found_dna = ""

    for sq in [search_sequence, str(Seq(search_sequence).reverse_complement())]:
        query = pyopa.Sequence(sq)
        target = pyopa.Sequence(seq)

        aln = pyopa.align_strings(query, 
            target, env, is_global=False)
        
        a1 = np.array(list(str(aln[0])))
        a2 = np.array(list(str(aln[1])))

        if len(a1)<1: continue

        aln_cols = ((np.cumsum(a1!=gap) * np.cumsum(a2!=gap)) * 
            np.flip(np.cumsum(np.flip(a1)!=gap) * np.cumsum(np.flip(a2)!=gap)))!=0
        Ns = (a1=="N") | (a2=="N")
        pwid = np.sum(((a1==a2) | Ns) & aln_cols)/(1.0*np.sum(aln_cols))

        if (pwid>=pairwise_id_thresh) and (len(a1)/(1.0*len(seq)) >= prop_match):
            keep = (a1!=gap) & (a1!="_")
            found_dna="".join(list(a1[keep]))

    return found_dna


def translate_to_match(hit, target_prot):

    if hit == "": return ""

    # translate in all 6 frames splitting on unknown
    dna_seqs = [hit, reverse_complement(hit)]

    proteins = [
        translate(s[i:].ljust(len(s[i:]) + (3 - len(s[i:]) % 3), 'N'))
        for i in range(3) for s in dna_seqs
    ]

    search_set = set(
        [target_prot[i:i + 3] for i in range(len(target_prot) - 2)])

    alignments = []
    for target_sequence in proteins:
        query_set = set([
            target_sequence[i:i + 3] for i in range(len(target_sequence) - 2)
        ])
        alignments.append(
            (target_sequence, len(search_set.intersection(query_set))))

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
