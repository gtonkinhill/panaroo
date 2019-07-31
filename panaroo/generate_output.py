from joblib import Parallel, delayed
import networkx as nx
from collections import defaultdict
import numpy as np
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import itertools as iter
from tqdm import tqdm

from panaroo.generate_alignments import *


def generate_roary_gene_presence_absence(G, mems_to_isolates, output_dir):

    # arange isolates
    isolates = []
    mems_to_index = {}
    for i, mem in enumerate(mems_to_isolates):
        isolates.append(mems_to_isolates[mem])
        mems_to_index[str(mem)] = i

    # generate file
    with open(output_dir + "gene_presence_absence.csv", 'w') as csv_outfile, \
         open(output_dir + "gene_presence_absence.Rtab", 'w') as Rtab_outfile:
        header = [
            "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
            "No. sequences", "Avg sequences per isolate", "Genome Fragment",
            "Order within Fragment", "Accessory Fragment",
            "Accessory Order with Fragment", "QC", "Min group size nuc",
            "Max group size nuc", "Avg group size nuc"
        ] + isolates
        csv_outfile.write(",".join(header) + "\n")
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")

        # Iterate through coponents writing out to file
        used_gene_names = set([""])
        unique_id_count = 0
        frag = 0
        for component in nx.connected_components(G):
            frag += 1
            count = 0
            for node in component:
                count += 1
                name = '_'.join(
                    G.node[node]['annotation'].strip().strip(';').split(';'))
                name = ''.join(e for e in name if e.isalnum() or e == "_")
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
                entry.append((1.0 * len(G.node[node]['seqIDs'])) /
                             len(set(G.node[node]['members'])))
                entry.append(frag)
                entry.append(count)
                entry += ["", "", ""]
                entry.append(np.min(G.node[node]['lengths']))
                entry.append(np.max(G.node[node]['lengths']))
                entry.append(np.mean(G.node[node]['lengths']))
                pres_abs = [""] * len(isolates)
                for seq in G.node[node]['seqIDs']:
                    sample_id = mems_to_index["_".join(seq.split("_")[:-2])]
                    if pres_abs[sample_id]=="": #ensures we only take the first one
                        pres_abs[sample_id] = seq
                entry += pres_abs
                csv_outfile.write(",".join([str(e) for e in entry]) + "\n")
                Rtab_outfile.write(entry[0] + "\t")
                Rtab_outfile.write("\t".join(
                    (["0" if e == "" else "1" for e in pres_abs])) + "\n")

    return G


def generate_pan_genome_reference(G, output_dir, split_paralogs=False):

    # need to treat paralogs differently?
    centroids = set()
    records = []

    for node in G.nodes():
        if not split_paralogs and G.node[node]['centroid'].split(";")[0] in centroids:
            continue
        records.append(
            SeqRecord(Seq(G.node[node]['dna'].split(";")[0], generic_dna),
                      id=G.node[node]['centroid'].split(";")[0],
                      description=""))
        centroids.add(G.node[node]['centroid'].split(";")[0])

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
                p = G[edge[0]][edge[1]]['weight'] / (1.0 *
                                                     G.node[node]['size'])
                entropy -= p * np.log(p)
            outfile.write(",".join([
                G.node[node]['name'], G.node[node]['annotation'],
                str(G.node[node]['size']),
                str(G.degree[node]), "{:.5f}".format(entropy)
            ]) + "\n")


def generate_common_struct_presence_absence(G,
                                            output_dir,
                                            mems_to_isolates,
                                            min_variant_support=2):

    # arange isolates
    isolates = []
    members = []
    for mem in mems_to_isolates:
        isolates.append(mems_to_isolates[mem])
        members.append(mem)

    struct_variants = {}
    for node in G.nodes():
        if G.degree[node] < 3: continue  #skip as linear
        for path in iter.combinations(G.edges(node), 2):
            in_both = (set(G[path[0][0]][path[0][1]]['members'])
                       & set(G[path[1][0]][path[1][1]]['members']))
            if len(in_both) >= min_variant_support:
                struct_variants[(path[0][0], path[0][1], path[1][1])] = in_both

    header = []
    for variant in struct_variants:
        header.append("-".join([
            G.node[variant[1]]['name'], G.node[variant[0]]['name'],
            G.node[variant[2]]['name']
        ]))

    with open(output_dir + "struct_presence_absence.Rtab", 'w') as Rtab_outfile:
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")
        for h, variant in zip(header, struct_variants):
            variant_calls = [h]
            for member in members:
                if member in struct_variants[variant]:
                    variant_calls.append("1")
                else:
                    variant_calls.append("0")
            Rtab_outfile.write("\t".join(variant_calls) + "\n")

    return


def generate_pan_genome_alignment(G, temp_dir, output_dir, threads, aligner,
                                  isolates):
    #Make a folder for the output alignments
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        None
    #Multithread writing gene sequences to disk (temp directory) so aligners can find them
    unaligned_sequence_files = Parallel(n_jobs=threads)(
        delayed(output_sequence)(G.node[x], isolates, temp_dir, output_dir)
        for x in tqdm(G.nodes()))
    #Get Biopython command calls for each output gene sequences
    commands = [
        get_alignment_commands(fastafile, output_dir, aligner, threads)
        for fastafile in unaligned_sequence_files
    ]
    #Run these commands in a multi-threaded way
    multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                          threads, aligner)
    return


def get_core_gene_nodes(G, threshold, num_isolates):
    #Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes():
        if float(G.node[node]["size"]) / float(num_isolates) > threshold:
            core_nodes.append(node)
    return core_nodes


def concatenate_core_genome_alignments(core_names, output_dir):

    alignments_dir = output_dir + "/aligned_gene_sequences/"
    #Open up each alignment that is assosciated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]

    #Read in all these alginemnts
    gene_alignments = []
    gene_names = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, 'fasta')
        gene_dict = {}
        for record in alignment:
            genome_id = record.id.split(";")[0]
            gene_dict[genome_id] = (record.id, record.seq)
            gene_length = len(record.seq)
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))

    #Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    #Write out the two output files
    SeqIO.write(isolate_aln, output_dir + 'core_gene_alignment.aln', 'fasta')

    write_alignment_header(gene_alignments, output_dir)
    return core_filenames


def generate_core_genome_alignment(G, temp_dir, output_dir, threads, aligner,
                                   isolates, threshold, num_isolates):
    #Make a folder for the output alignments TODO: decide whether or not to keep these
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        None
    #Get core nodes
    core_genes = get_core_gene_nodes(G, threshold, num_isolates)
    core_gene_names = [G.node[x]["name"] for x in core_genes]
    #Output core node sequences
    unaligned_sequence_files = Parallel(n_jobs=threads)(
        delayed(output_sequence)(G.node[x], isolates, temp_dir, output_dir)
        for x in tqdm(core_genes))
    #Get alignment commands
    commands = [
        get_alignment_commands(fastafile, output_dir, aligner, threads)
        for fastafile in unaligned_sequence_files
    ]
    #Run alignment commands
    multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                          threads, aligner)
    #Concatenate them together to produce the two output files
    concatenate_core_genome_alignments(core_gene_names, output_dir)
    return
