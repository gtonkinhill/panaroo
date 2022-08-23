import networkx as nx
import io, sys
from collections import defaultdict, Counter
import numpy as np
from Bio.Seq import translate, reverse_complement, Seq
from Bio import SeqIO
from .cdhit import align_dna_cdhit
from .isvalid import del_dups
from joblib import Parallel, delayed
import os
import gffutils as gff
from io import StringIO
import edlib
from .merge_nodes import delete_node, remove_member_from_node
from tqdm import tqdm
import re


# @profile
def find_missing(G,
                 gff_file_handles,
                 dna_seq_file,
                 prot_seq_file,
                 gene_data_file,
                 merge_id_thresh,
                 search_radius,
                 prop_match,
                 pairwise_id_thresh,
                 n_cpu,
                 remove_by_consensus=False,
                 verbose=True):

    # Iterate over each genome file checking to see if any missing accessory genes
    #  can be found.

    # generate mapping between internal nodes and gff ids
    id_to_gff = {}
    with open(gene_data_file, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.split(",")
            if line[2] in id_to_gff:
                raise NameError("Duplicate internal ids!")
            id_to_gff[line[2]] = line[3]

    # identify nodes that have been merged at the protein level
    merged_ids = {}
    for node in G.nodes():
        if (len(G.nodes[node]['centroid']) >
                1) or (G.nodes[node]['mergedDNA']):
            for sid in sorted(G.nodes[node]['seqIDs']):
                merged_ids[sid] = node

    merged_nodes = defaultdict(dict)
    with open(gene_data_file, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.split(",")
            if line[2] in merged_ids:
                mem = int(sid.split("_")[0])
                if merged_ids[line[2]] in merged_nodes[mem]:
                    merged_nodes[mem][merged_ids[line[2]]] = G.nodes[
                        merged_ids[line[2]]]["dna"][G.nodes[merged_ids[
                            line[2]]]['maxLenId']]
                else:
                    merged_nodes[mem][merged_ids[line[2]]] = line[5]

    # iterate through nodes to identify accessory genes for searching
    # these are nodes missing a member with at least one neighbour that has that member
    n_searches = 0
    search_list = defaultdict(lambda: defaultdict(set))
    conflicts = defaultdict(set)
    for node in G.nodes():
        for neigh in G.neighbors(node):
            # seen_mems = set()
            for sid in sorted(G.nodes[neigh]['seqIDs']):
                member = int(sid.split("_")[0])

                conflicts[member].add((neigh, id_to_gff[sid]))
                if member not in G.nodes[node]['members']:
                    if len(G.nodes[node]["dna"][G.nodes[node]
                                                ['maxLenId']]) <= 0:
                        print(G.nodes[node]["dna"])
                        raise NameError("Problem!")
                    search_list[member][node].add(
                        (G.nodes[node]["dna"][G.nodes[node]['maxLenId']],
                         id_to_gff[sid]))

                    n_searches += 1

    if verbose:
        print("Number of searches to perform: ", n_searches)
        print("Searching...")

    all_hits, all_node_locs, max_seq_lengths = zip(*Parallel(n_jobs=n_cpu)(
        delayed(search_gff)(search_list[member],
                            conflicts[member],
                            gff_handle,
                            merged_nodes=merged_nodes[member],
                            search_radius=search_radius,
                            prop_match=prop_match,
                            pairwise_id_thresh=pairwise_id_thresh,
                            merge_id_thresh=merge_id_thresh)
        for member, gff_handle in tqdm(enumerate(gff_file_handles),
                                       disable=(not verbose))))

    if verbose:
        print("translating hits...")

    hits_trans_dict = {}
    for member, hits in enumerate(all_hits):
        hits_trans_dict[member] = Parallel(n_jobs=n_cpu)(
            delayed(translate_to_match)(hit[1], G.nodes[hit[0]]["protein"][0])
            for hit in hits)

    # remove nodes that conflict (overlap)
    nodes_by_size = sorted([(G.nodes[node]['size'], node)
                            for node in G.nodes()],
                           reverse=True)
    nodes_by_size = [n[1] for n in nodes_by_size]
    member = 0
    bad_node_mem_pairs = set()
    bad_nodes = set()
    for node_locs, max_seq_length in zip(all_node_locs, max_seq_lengths):
        seq_coverage = defaultdict(
            lambda: np.zeros(max_seq_length + 2, dtype=bool))

        for node in nodes_by_size:
            if node in bad_nodes: continue
            if node not in node_locs: continue
            contig_id = node_locs[node][0]
            loc = node_locs[node][1]

            if np.sum(seq_coverage[contig_id][loc[0]:loc[1]]) >= (
                    0.5 * (max(G.nodes[node]['lengths']))):
                if member in G.nodes[node]['members']:
                    remove_member_from_node(G, node, member)
                # G.nodes[node]['members'].remove(str(member))
                # G.nodes[node]['size'] -= 1
                bad_node_mem_pairs.add((node, member))
            else:
                seq_coverage[contig_id][loc[0]:loc[1]] = True
        member += 1

    for node in G.nodes():
        if len(G.nodes[node]['members']) <= 0:
            bad_nodes.add(node)
    for node in bad_nodes:
        if node in G.nodes():
            delete_node(G, node)

    # remove by consensus
    if remove_by_consensus:
        if verbose:
            print("removing by consensus...")
        node_hit_counter = Counter()
        for member, hits in enumerate(all_hits):
            for node, dna_hit in hits:
                if dna_hit == "": continue
                if node in bad_nodes: continue
                if (node, member) in bad_node_mem_pairs: continue
                node_hit_counter[node] += 1
        for node in G:
            if node_hit_counter[node] > G.nodes[node]['size']:
                bad_nodes.add(node)
        for node in bad_nodes:
            if node in G.nodes():
                delete_node(G, node)

    if verbose:
        print("Updating output...")

    n_found = 0
    with open(dna_seq_file, 'a') as dna_out:
        with open(prot_seq_file, 'a') as prot_out:
            with open(gene_data_file, 'a') as data_out:
                for member, (hits, node_locs) in enumerate(zip(all_hits, all_node_locs)):
                    i = -1
                    for node, dna_hit in hits:
                        i += 1
                        if dna_hit == "": continue
                        if node in bad_nodes: continue
                        if (node, member) in bad_node_mem_pairs: continue
                        
                        hit_protein = hits_trans_dict[member][i]
                        hit_strand = '+' if node_locs[node][1][2]==0 else '-'
                        G.nodes[node]['members'].add(member)
                        G.nodes[node]['size'] += 1
                        G.nodes[node]['dna'] = del_dups(G.nodes[node]['dna'] +
                                                        [dna_hit])
                        dna_out.write(">" + str(member) + "_refound_" +
                                      str(n_found) + "\n" + dna_hit + "\n")
                        G.nodes[node]['protein'] = del_dups(
                            G.nodes[node]['protein'] + [hit_protein])
                        prot_out.write(">" + str(member) + "_refound_" +
                                       str(n_found) + "\n" + hit_protein +
                                       "\n")
                        data_out.write(",".join([
                            os.path.splitext(
                                os.path.basename(
                                    gff_file_handles[member]))[0], node_locs[node][0],
                            str(member) + "_refound_" + str(n_found),
                            str(member) + "_refound_" +
                            str(n_found), hit_protein, dna_hit, "", 
                            "location:" + str(node_locs[node][1][0]) + '-' +
                            str(node_locs[node][1][1]) + ';strand:' + hit_strand
                        ]) + "\n")
                        G.nodes[node]['seqIDs'] |= set(
                            [str(member) + "_refound_" + str(n_found)])
                        n_found += 1

    if verbose:
        print("Number of refound genes: ", n_found)

    return (G)


def search_gff(node_search_dict,
               conflicts,
               gff_handle_name,
               merged_nodes,
               search_radius=10000,
               prop_match=0.2,
               pairwise_id_thresh=0.95,
               merge_id_thresh=0.7,
               n_cpu=1):

    gff_handle = open(gff_handle_name, 'r')

    # sort sets to fix order
    conflicts = sorted(conflicts)
    for node in node_search_dict:
        node_search_dict[node] = sorted(node_search_dict[node])

    split = gff_handle.read().replace(',', '').split("##FASTA\n")
    node_locs = {}

    if len(split) != 2:
        raise NameError("File does not appear to be in GFF3 format!")

    # load fasta
    contigs = {}
    max_seq_len = 0
    with StringIO(split[1]) as temp_fasta:
        for record in SeqIO.parse(temp_fasta, 'fasta'):
            contigs[record.id] = np.array(list(str(record.seq)))
            max_seq_len = max(max_seq_len, len(contigs[record.id]))

    # load gff annotation
    parsed_gff = gff.create_db("\n".join(
        [l for l in split[0].splitlines() if '##sequence-region' not in l]),
                               dbfn=":memory:",
                               force=True,
                               keep_order=True,
                               from_string=True)

    # mask regions that already have genes and convert back to string
    seen = set()
    for node, geneid in conflicts:
        gene = parsed_gff[geneid]
        start = min(gene.start, gene.end)
        end = max(gene.start, gene.end)

        if node in merged_nodes:
            db_seq = contigs[gene[0]][max(0, (start -
                                              search_radius)):(end +
                                                               search_radius)]
            db_seq = "".join(list(db_seq))

            hit, loc = search_dna(db_seq,
                                  merged_nodes[node],
                                  prop_match=(end - start) /
                                  float(len(merged_nodes[node])),
                                  pairwise_id_thresh=merge_id_thresh,
                                  refind=False)

            # update location
            loc[0] = loc[0] + max(0, (start - search_radius))
            loc[1] = loc[1] + max(0, (start - search_radius))
            node_locs[node] = [gene[0], loc]
        else:
            node_locs[node] = [gene[0], [start - 1, end]]

    for node, geneid in conflicts:
        gene = parsed_gff[geneid]
        start = min(gene.start, gene.end)
        end = max(gene.start, gene.end)
        # contigs[gene[0]][(start - 1):end] = "X"

        if (gene[0], start - 1, end) in seen:
            raise NameError("Duplicate entry!!!")
        seen.add((gene[0], start - 1, end))

    for sid in contigs:
        contigs[sid] = "".join(list(contigs[sid]))

    # search for matches
    hits = []
    for node in node_search_dict:
        best_hit = ""
        best_loc = None
        for search in node_search_dict[node]:
            gene = parsed_gff[search[1]]
            start = min(gene.start, gene.end)
            end = max(gene.start, gene.end)
            db_seq = contigs[gene[0]][max(0, (start -
                                              search_radius)):(end +
                                                               search_radius)]

            hit, loc = search_dna(db_seq,
                                  search[0],
                                  prop_match,
                                  pairwise_id_thresh,
                                  refind=True)
            # update location
            loc[0] = loc[0] + max(0, (start - search_radius))
            loc[1] = loc[1] + max(0, (start - search_radius))

            if len(hit) > len(best_hit):
                best_hit = hit
                best_loc = [gene[0], loc]

        hits.append((node, best_hit))
        if (best_loc is not None) and (best_hit != ""):
            node_locs[node] = best_loc

    gff_handle.close()

    return [hits, node_locs, max_seq_len]


def repl(m):
    return ('X' * len(m.group()))


def search_dna(db_seq, search_sequence, prop_match, pairwise_id_thresh,
               refind):
    found_dna = ""
    start = None
    end = None
    max_hit = 0
    loc = [0, 0]

    # found=False
    # if search_sequence=="":
    #     if refind:
    #         print(">>>>>>>>>>>>>")
    #         print(db_seq)
    #         found=True

    added_E_len = int(len(search_sequence) / 2)

    for i, db in enumerate([db_seq, str(Seq(db_seq).reverse_complement())]):

        # add some Ns at the start and end to deal with fragments at the end of contigs
        db = "E" * added_E_len + db + "E" * added_E_len

        aln = edlib.align(search_sequence,
                          db,
                          mode="HW",
                          task='path',
                          k=10 * len(search_sequence),
                          additionalEqualities=[
                              ('A', 'N'),
                              ('C', 'N'),
                              ('G', 'N'),
                              ('T', 'N'),
                              ('A', 'E'),
                              ('C', 'E'),
                              ('G', 'E'),
                              ('T', 'E'),
                          ])

        # remove trailing inserts
        cig = re.split(r'(\d+)', aln['cigar'])[1:]
        if cig[-1] == "I":
            aln['editDistance'] -= int(cig[-2])
        if cig[1] == "I":
            aln['editDistance'] -= int(cig[0])

        if aln['editDistance'] == -1:
            start = -1
        else:
            # take hit that is closest to the centre of the neighbouring gene
            centre = len(db) / 2.0
            tloc = min(aln['locations'],
                       key=lambda x: min(centre - x[0], centre - x[1]))
            start = tloc[0]
            end = tloc[1] + 1

        # if found:
        #     print(aln)
        #     print(start, end)

        # skip if nothing was found
        if start == -1: continue

        possible_dbs = [db]
        if db.find("NNNNNNNNNNNNNNNNNNNN") != -1:
            possible_dbs += [
                re.sub("^[ACGTEX]{0,}NNNNNNNNNNNNNNNNNNNN", repl, db, 1),
                re.sub("NNNNNNNNNNNNNNNNNNNN[ACGTEX]{0,}$", repl, db, 1)
            ]

        for posdb in possible_dbs:
            # skip if alignment is too short
            n_X = posdb[start:end].count("X")
            n_E = posdb[start:end].count("E")

            aln_length = float(end - start - n_X - n_E)
            if (aln_length / len(search_sequence)) <= prop_match: continue
            if (posdb[start:end].count("A") + posdb[start:end].count("C") +
                    posdb[start:end].count("G") + posdb[start:end].count("T")
                ) / len(search_sequence) <= prop_match:
                continue

            # determine an approximate percentage identity
            pid = 1.0 - (aln['editDistance'] - n_X) / (1.0 * aln_length)

            # skip if identity below threshold
            if pid <= pairwise_id_thresh: continue

            # if found:
            #     print("aln_length:", aln_length)
            #     print("pid:", pid)

            if max_hit < (pid * aln_length):
                found_dna = posdb[start:end]
                max_hit = (pid * aln_length)
                if i == 0:
                    loc = [start, end]
                else:
                    loc = [len(posdb) - tloc[1] - 1, len(posdb) - tloc[0]]
                loc = [
                    max(0,
                        min(loc) - added_E_len),
                    min(max(loc) - added_E_len, len(db_seq)),
                    i
                ]

    # if found:
    #     print(found_dna)
    #     print(loc)
    #     print("<<<<<<<<<<<<<<<<<<")
    seq = found_dna.replace('X', 'N').replace('E', 'N')
    seq = seq.strip('N')

    # if i==1:
    #     seq = str(Seq(seq).reverse_complement())

    return seq, loc


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

    return (prot[0])


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
