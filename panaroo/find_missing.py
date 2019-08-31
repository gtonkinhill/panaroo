import networkx as nx
import io, sys
from collections import defaultdict, Counter
import numpy as np
from Bio.Seq import translate, reverse_complement, Seq
from Bio import SeqIO
from panaroo.cdhit import align_dna_cdhit
from joblib import Parallel, delayed
import os
import gffutils as gff
from io import StringIO
import edlib
from .merge_nodes import delete_node, remove_member_from_node
from tqdm import tqdm


def find_missing(G, gff_file_handles, dna_seq_file, prot_seq_file, 
                gene_data_file,
                n_cpu, remove_by_consensus=False):

    # Iterate over each genome file checking to see if any missing accessory genes
    #  can be found.

    # generate mapping between internal nodes and gff ids
    id_to_gff = {}
    with open(gene_data_file, 'r') as infile:
        next(infile)
        for line in infile:
            line=line.split(",")
            if line[2] in id_to_gff:
                raise NameError("Duplicate internal ids!")
            id_to_gff[line[2]] = line[3]

    # iterate through nodes to identify accessory genes for searching
    # these are nodes missing a member with at least one neighbour that has that member
    n_searches = 0
    search_list = defaultdict(lambda: defaultdict(list))
    conflicts = defaultdict(set)
    for node in G.nodes():
        for neigh in G.neighbors(node):
            for sid in G.node[neigh]['seqIDs']:
                member = sid.split("_")[0]
                conflicts[int(member)].add(id_to_gff[sid])
                if member not in G.node[node]['members']:
                    if len(max(G.node[node]["dna"].split(";"), key=len))<=0:
                        print(G.node[node]["dna"])
                        raise NameError("Problem!")
                    search_list[int(member)][node].append((max(G.node[node]["dna"].split(";"), key=len), 
                                                     id_to_gff[sid]))
                    n_searches += 1

    print("Number of searches to perform: ", n_searches)
    print("Searching...")

    all_hits = Parallel(n_jobs=n_cpu)(
                delayed(search_gff)(
                    search_list[member], conflicts[member], gff_handle)
                    for member, gff_handle in tqdm(enumerate(gff_file_handles)))
    # ,search_radius, prop_match, pairwise_id_thresh, n_cpu)    

    print("translating hits...")

    hits_trans_dict = {}
    for member, hits in enumerate(all_hits):
        hits_trans_dict[member] = Parallel(n_jobs=n_cpu)(
                delayed(translate_to_match)(
                    hit[1], G.node[hit[0]]["protein"].split(";")[0])
                for hit in hits)

    bad_nodes = set()
    if remove_by_consensus:
        print("removing by consensus...")
        node_hit_counter = Counter()
        for hits in all_hits:
            for node, dna_hit in hits:
                if dna_hit=="": continue
                node_hit_counter[node] += 1        
        for node in G:
            if node_hit_counter[node]>G.node[node]['size']:
                bad_nodes.add(node)
        for node in bad_nodes:
            delete_node(G, node)


    print("Updating output...")

    n_found = 0
    with open(dna_seq_file, 'a') as dna_out:
        with open(prot_seq_file, 'a') as prot_out:
            for member, hits in enumerate(all_hits):
                i = -1
                for node, dna_hit in hits:
                    i+=1
                    if dna_hit=="": continue
                    if node in bad_nodes: continue
                    hit_protein = hits_trans_dict[member][i]
                    G.node[node]['members'] += [str(member)]
                    G.node[node]['size'] += 1
                    G.node[node]['dna'] = ";".join(
                        set(G.node[node]['dna'].split(";") + [dna_hit]))
                    dna_out.write(">" + str(member) + "_refound_" +
                                    str(n_found) + "\n" + dna_hit + "\n")
                    G.node[node]['protein'] = ";".join(
                        set(G.node[node]['protein'].split(";") +
                            [hit_protein]))
                    prot_out.write(">" + str(member) + "_refound_" +
                                    str(n_found) + "\n" + hit_protein + "\n")
                    G.node[node]['seqIDs'] += [
                        str(member) + "_refound_" + str(n_found)
                    ]
                    n_found += 1

    print("Number of refound genes: ", n_found)


    return(G)



def search_gff(node_search_dict,
                   conflicts,
                   gff_handle,
                   search_radius=10000,
                   prop_match=0.2,
                   pairwise_id_thresh=0.95,
                   n_cpu=1):

    # reset file handle to the beginning
    gff_handle.seek(0)
    split = gff_handle.read().split("##FASTA\n")

    if len(split) != 2:
        raise NameError("File does not appear to be in GFF3 format!")

    # load fasta
    contigs = {}
    with StringIO(split[1]) as temp_fasta:
        for record in SeqIO.parse(temp_fasta, 'fasta'):
            contigs[record.id] = np.array(list(str(record.seq)))

    # load gff annotation
    parsed_gff = gff.create_db("\n".join(
        [l for l in split[0].splitlines() if '##sequence-region' not in l]),
                               dbfn=":memory:",
                               force=True,
                               keep_order=True,
                               from_string=True)

    # mask regions that already have genes and convert back to string
    for geneid in conflicts:
        gene = parsed_gff[geneid]
        start = min(gene.start, gene.end)
        end = max(gene.start, gene.end)
        contigs[gene[0]][(start-1):end] = "X"
    for sid in contigs:
        contigs[sid] = "".join(list(contigs[sid]))

    # search for matches
    hits = []
    for node in node_search_dict:
        best_hit = ""
        for search in node_search_dict[node]:
            gene = parsed_gff[search[1]]
            start = min(gene.start, gene.end)
            end = max(gene.start, gene.end)
            db_seq = contigs[gene[0]][max(0,(start-search_radius)):(end+search_radius)]

            hit = search_dna(db_seq, search[0], prop_match, pairwise_id_thresh)
            if len(hit)>len(best_hit):
                best_hit = hit
        hits.append((node, best_hit))

    return hits


def search_dna(db_seq, search_sequence, prop_match, pairwise_id_thresh):
    found_dna = ""
    start = None
    end = None

    # found=False
    # if search_sequence=="TTGCGTGGCCTCGTAATCGCTATATCTACTATTATGTCGCCTGAAACCCACTTCGCGGTGGGTTTTTTGTTGTCAGGAGTTTTAATAATGGTGGGGCCGTCTGTGCGCGTGAATGAATGGTTCAGCGCGTATGCGATGGCGGGTATGGCTTACAGCCGTGTGTCGACTTTCTCCGGGGATTATCTCCGCGTAACTGACAACAAGGGAAGGTGCGAATAA":
    #     print(">>>>>>>>>>>>>")
    #     print(db_seq)
    #     found=True

    for db in [db_seq, str(Seq(db_seq).reverse_complement())]:

        aln = edlib.align(search_sequence, 
            db, mode="HW", task='locations', k=10*len(search_sequence),
            additionalEqualities=[('A','N'), ('C','N'), ('G','N'), ('T','N')])
        
        if aln['editDistance']==-1:
            start = -1
        else:
            start = aln['locations'][0][0]
            end = aln['locations'][0][1] + 1

        # if found:
        #     print(aln)
        #     print(start, end)

        # skip if nothing was found
        if start==-1: continue

        # skip if alignment is too short
        n_X = db[start:end].count("X")
        aln_length = float(end-start-n_X)
        if (aln_length/len(search_sequence)) <= prop_match: continue

        # determine an approximate percentage identity
        if (start==0) or  (end==len(db_seq)):
            pid = 1.0-max(0.0, (aln['editDistance'] - n_X - (len(search_sequence)-aln_length)))/aln_length
        else:
            pid = 1.0-(aln['editDistance']-n_X)/(1.0*len(search_sequence))

        # skip if identity below threshold
        if pid <= pairwise_id_thresh: continue

        # if found:
        #     print("aln_length:", aln_length)
        #     print("pid:", pid)
        
        found_dna=db[start:end]

    # if found:
    #     print(found_dna)
    #     print("<<<<<<<<<<<<<<<<<<")

    return found_dna.replace('X','N')

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

