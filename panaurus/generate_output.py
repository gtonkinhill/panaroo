import networkx as nx
from collections import defaultdict
import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import itertools as iter


def generate_roary_gene_presence_absence(G, file_names, dna_file, output_dir):

    # clean file names
    file_names = [
        os.path.splitext(os.path.basename(f.name))[0] for f in file_names
    ]
    file_names = [f.replace(",", "") for f in file_names]

    # calculate node average, min and max lengths
    seq_lengths = {}
    for rec in SeqIO.parse(dna_file, "fasta"):
        seq_lengths[str(rec.id)] = len(rec.seq)
    node_lengths = defaultdict(list)
    for node in G.nodes():
        for seq in G.node[node]['seqIDs']:
            node_lengths[node].append(seq_lengths[seq])

    # generate file
    with open(output_dir + "gene_presence_absence.csv", 'w') as outfile:
        header = [
            "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
            "No. sequences", "Avg sequences per isolate", "Genome Fragment",
            "Order within Fragment", "Accessory Fragment",
            "Accessory Order with Fragment", "QC", "Min group size nuc",
            "Max group size nuc", "Avg group size nuc"
        ] + file_names
        outfile.write(",".join(header) + "\n")

        # Iterate through coponents writing out to file
        used_gene_names = set([""])
        unique_id_count = 0
        frag = 0
        for component in nx.connected_components(G):
            frag += 1
            count = 0
            for node in component:
                count += 1
                name = G.node[node]['annotation'].strip()
                if name not in used_gene_names:
                    entry = [name]
                    used_gene_names.add(name)
                    G.node[node]['name'] = name
                else:
                    G.node[node]['name'] = "group_" + str(unique_id_count)
                    entry = [G.node[node]['name']]
                    unique_id_count += 1
                entry.append(G.node[node]['annotation'])
                entry.append(G.node[node]['description'])
                entry.append(len(G.node[node]['seqIDs']))
                entry.append(len(set(G.node[node]['members'])))
                entry.append((1.0 * len(G.node[node]['seqIDs'])) / len(
                    set(G.node[node]['members'])))
                entry.append(frag)
                entry.append(count)
                entry += ["", "", ""]
                entry.append(np.min(node_lengths[node]))
                entry.append(np.max(node_lengths[node]))
                entry.append(np.mean(node_lengths[node]))
                pres_abs = [""] * len(file_names)
                for seq in G.node[node]['seqIDs']:
                    sample_id = int(seq.split("_")[0])
                    pres_abs[sample_id] = seq
                entry += pres_abs
                outfile.write(",".join([str(e) for e in entry]) + "\n")

    return G


def generate_pan_genome_reference(G, output_dir, split_paralogs=False):

    # need to treat paralogs differently?
    centroids = set()
    records = []

    for node in G.nodes():
        if not split_paralogs and G.node[node]['centroid'] in centroids:
            continue
        records.append(
            SeqRecord(
                Seq(G.node[node]['dna'], generic_dna),
                id=G.node[node]['centroid'],
                description=""))
        centroids.add(G.node[node]['centroid'])

    with open(output_dir + "pan_genome_reference.fa", 'w') as outfile:
        SeqIO.write(records, outfile, "fasta")

    return


# TODO: come with a nice weighting to account for the number of observations
def generate_gene_mobility(G, output_dir):

    with open(output_dir + "gene_mobility.csv", 'w') as outfile:
        outfile.write("gene_id,annotation,count,degree,entropy\n")
        for node in G.nodes():
            entropy = 0
            for edge in G.edges(node):
                p = G[edge[0]][edge[1]]['weight'] / (
                    1.0 * G.node[node]['size'])
                entropy -= p * np.log(p)
            outfile.write(",".join([
                G.node[node]['name'], G.node[node]['annotation'],
                str(G.node[node]['size']),
                str(G.degree[node]), "{:.5f}".format(entropy)
            ]) + "\n")


def generate_common_struct_presence_absence(G,
                                            output_dir,
                                            n_members,
                                            min_variant_support=2):

    struct_variants = {}
    for node in G.nodes():
        if G.degree[node] < 3: continue  #skip as linear
        for path in iter.combinations(G.edges(node), 2):
            in_both = (set(G[path[0][0]][path[0][1]]['members']) & set(
                G[path[1][0]][path[1][1]]['members']))
            if len(in_both) >= min_variant_support:
                struct_variants[(path[0][0], path[0][1], path[1][1])] = in_both

    header = []
    for variant in struct_variants:
        header.append("-".join([
            G.node[variant[1]]['name'], G.node[variant[0]]['name'],
            G.node[variant[2]]['name']
        ]))

    with open(output_dir + "struct_presence_absence.csv", 'w') as outfile:
        outfile.write(",".join(header) + "\n")
        for member in range(n_members):
            variant_calls = []
            for variant in struct_variants:
                if str(member) in struct_variants[variant]:
                    variant_calls.append("1")
                else:
                    variant_calls.append("0")
            outfile.write(",".join(variant_calls) + "\n")

    return
