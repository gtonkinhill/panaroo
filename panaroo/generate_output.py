from joblib import Parallel, delayed
import networkx as nx
from collections import defaultdict
import numpy as np
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools as iter
from tqdm import tqdm

from .generate_alignments import *


def generate_roary_gene_presence_absence(
    G, mems_to_isolates, orig_ids, ids_len_stop, output_dir
):

    # arange isolates
    isolates = []
    mems_to_index = {}
    for i, mem in enumerate(mems_to_isolates):
        isolates.append(mems_to_isolates[mem])
        mems_to_index[str(mem)] = i

    # generate file
    with open(
        output_dir + "gene_presence_absence_roary.csv", "w"
    ) as roary_csv_outfile, open(
        output_dir + "gene_presence_absence.csv", "w"
    ) as csv_outfile, open(
        output_dir + "gene_presence_absence.Rtab", "w"
    ) as Rtab_outfile:
        header = [
            "Gene",
            "Non-unique Gene name",
            "Annotation",
            "No. isolates",
            "No. sequences",
            "Avg sequences per isolate",
            "Genome Fragment",
            "Order within Fragment",
            "Accessory Fragment",
            "Accessory Order with Fragment",
            "QC",
            "Min group size nuc",
            "Max group size nuc",
            "Avg group size nuc",
        ] + isolates
        roary_csv_outfile.write(",".join(header) + "\n")
        csv_outfile.write(",".join(header[:3] + isolates) + "\n")
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")

        # Iterate through coponents writing out to file
        used_gene_names = set([""])
        unique_id_count = 0
        frag = 0
        entry_list = []
        entry_ext_list = []
        pres_abs_list = []
        entry_sizes = []
        entry_count = 0
        for component in nx.connected_components(G):
            frag += 1
            count = 0
            for node in component:
                count += 1
                len_mode = max(
                    G.nodes[node]["lengths"], key=G.nodes[node]["lengths"].count
                )
                name = "~~~".join(
                    [
                        gn
                        for gn in G.nodes[node]["annotation"]
                        .strip()
                        .strip(";")
                        .split(";")
                        if gn != ""
                    ]
                )
                name = "".join(e for e in name if e.isalnum() or e in ["_", "~"])
                if name.lower() not in used_gene_names:
                    entry = [name]
                    used_gene_names.add(name.lower())
                    G.nodes[node]["name"] = name
                else:
                    G.nodes[node]["name"] = "group_" + str(unique_id_count)
                    entry = [G.nodes[node]["name"]]
                    unique_id_count += 1
                entry.append(G.nodes[node]["annotation"])
                entry.append(G.nodes[node]["description"])
                entry.append(G.nodes[node]["size"])
                entry.append(len(G.nodes[node]["seqIDs"]))
                entry.append(
                    (1.0 * len(G.nodes[node]["seqIDs"])) / G.nodes[node]["size"]
                )
                entry.append(frag)
                entry.append(count)
                entry += ["", "", ""]
                entry.append(np.min(G.nodes[node]["lengths"]))
                entry.append(np.max(G.nodes[node]["lengths"]))
                entry.append(np.mean(G.nodes[node]["lengths"]))
                pres_abs = [""] * len(isolates)
                pres_abs_ext = [""] * len(isolates)
                entry_size = 0
                for seq in G.nodes[node]["seqIDs"]:
                    sample_id = mems_to_index["_".join(seq.split("_")[:-2])]
                    if pres_abs[sample_id] == "":  # ensures we only take the first one
                        if seq in orig_ids:
                            pres_abs[sample_id] = orig_ids[seq]
                            pres_abs_ext[sample_id] = orig_ids[seq]
                        else:
                            pres_abs[sample_id] = seq
                            pres_abs_ext[sample_id] = seq
                        entry_size += 1
                    else:
                        # this is similar to PIRATE output
                        if seq in orig_ids:
                            pres_abs[sample_id] += ";" + orig_ids[seq]
                            pres_abs_ext[sample_id] += ";" + orig_ids[seq]
                        else:
                            pres_abs[sample_id] += ";" + seq
                            pres_abs_ext[sample_id] += ";" + seq
                    if (abs(ids_len_stop[seq][0] - len_mode) / len_mode) > (
                        0.05 * len_mode
                    ):
                        pres_abs_ext[sample_id] += "_len"
                    if ids_len_stop[seq][1]:
                        pres_abs_ext[sample_id] += "_stop"

                entry += pres_abs
                entry_list.append(entry)
                entry_ext_list.append(entry[:3] + pres_abs_ext)
                pres_abs_list.append(pres_abs)
                entry_sizes.append((entry_size, entry_count))
                entry_count += 1

        # sort so that the most common genes are first (as in roary)
        entry_sizes = sorted(entry_sizes, reverse=True)
        for s, i in entry_sizes:
            roary_csv_outfile.write(",".join([str(e) for e in entry_list[i]]) + "\n")
            csv_outfile.write(",".join([str(e) for e in entry_ext_list[i]]) + "\n")
            Rtab_outfile.write(entry_list[i][0] + "\t")
            Rtab_outfile.write(
                "\t".join((["0" if e == "" else "1" for e in pres_abs_list[i]])) + "\n"
            )

    return G


def generate_pan_genome_reference(G, output_dir, split_paralogs=False):

    # need to treat paralogs differently?
    centroids = set()
    records = []

    for node in G.nodes():
        if not split_paralogs and G.nodes[node]["centroid"][0] in centroids:
            continue
        records.append(
            SeqRecord(
                Seq(max(G.nodes[node]["dna"], key=lambda x: len(x))),
                id=G.nodes[node]["name"],
                description="",
            )
        )
        for centroid in G.nodes[node]["centroid"]:
            centroids.add(centroid)

    with open(output_dir + "pan_genome_reference.fa", "w") as outfile:
        SeqIO.write(records, outfile, "fasta")

    return


def generate_common_struct_presence_absence(
    G, output_dir, mems_to_isolates, min_variant_support=2
):

    # arange isolates
    isolates = []
    members = []
    for mem in mems_to_isolates:
        isolates.append(mems_to_isolates[mem])
        members.append(mem)

    struct_variants = {}
    for node in G.nodes():
        if G.degree[node] < 3:
            continue  # skip as linear
        for path in iter.combinations(G.edges(node), 2):
            in_both = (
                G[path[0][0]][path[0][1]]["members"]
                & G[path[1][0]][path[1][1]]["members"]
            )
            if len(in_both) >= min_variant_support:
                struct_variants[(path[0][0], path[0][1], path[1][1])] = in_both

    header = []
    for variant in struct_variants:
        header.append(
            "-".join(
                [
                    G.nodes[variant[1]]["name"],
                    G.nodes[variant[0]]["name"],
                    G.nodes[variant[2]]["name"],
                ]
            )
        )

    with open(output_dir + "struct_presence_absence.Rtab", "w") as Rtab_outfile:
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


def generate_pan_genome_alignment(G, temp_dir, output_dir, threads, aligner, isolates):
    # Make a folder for the output alignments
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        None
    # Multithread writing gene sequences to disk (temp directory) so aligners can find them
    unaligned_sequence_files = Parallel(n_jobs=threads)(
        delayed(output_sequence)(G.nodes[x], isolates, temp_dir, output_dir)
        for x in tqdm(G.nodes())
    )
    # remove single sequence files
    unaligned_sequence_files = filter(None, unaligned_sequence_files)
    # Get Biopython command calls for each output gene sequences
    commands = [
        get_alignment_commands(fastafile, output_dir, aligner, threads)
        for fastafile in unaligned_sequence_files
    ]
    # Run these commands in a multi-threaded way
    multi_align_sequences(
        commands, output_dir + "aligned_gene_sequences/", threads, aligner
    )
    return


def get_core_gene_nodes(G, threshold, num_isolates):
    # Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes():
        # if G.nodes[node]['paralog']: continue
        if float(G.nodes[node]["size"]) / float(num_isolates) > threshold:
            core_nodes.append(node)
    return core_nodes


def concatenate_core_genome_alignments(core_names, output_dir):

    alignments_dir = output_dir + "/aligned_gene_sequences/"
    # Open up each alignment that is assosciated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [x for x in alignment_filenames if x.split(".")[0] in core_names]

    # Read in all these alignments
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, "fasta")
        gene_dict = {}
        for record in alignment:
            if record.id[:3] == "_R_":
                record.id = record.id[3:]
            genome_id = record.id.split(";")[0]
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)
            gene_length = len(record.seq)
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))
    # Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    # Write out the two output files
    SeqIO.write(isolate_aln, output_dir + "core_gene_alignment.aln", "fasta")
    write_alignment_header(gene_alignments, output_dir, "core_alignment_header.embl")

    # Write out filtered files
    conserved_alignments = set()
    for gene in gene_alignments:
        if is_conserved(list(gene[1].values()), 0.9):
            conserved_alignments.add(gene[0])

    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if gene[0] in conserved_alignments:
                if iso in gene[1]:
                    seq += gene[1][iso][1]
                else:
                    seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    # Write out the two output files
    SeqIO.write(isolate_aln, output_dir + "core_gene_alignment_filtered.aln", "fasta")
    write_alignment_header(
        [g for g in gene_alignments if g[0] in conserved_alignments],
        output_dir,
        "core_alignment_filtered_header.embl",
    )

    return core_filenames


def generate_core_genome_alignment(
    G, temp_dir, output_dir, threads, aligner, isolates, threshold, num_isolates
):
    # Make a folder for the output alignments TODO: decide whether or not to keep these
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        None
    # Get core nodes
    core_genes = get_core_gene_nodes(G, threshold, num_isolates)
    core_gene_names = [G.nodes[x]["name"] for x in core_genes]
    # Output core node sequences
    unaligned_sequence_files = Parallel(n_jobs=threads)(
        delayed(output_sequence)(G.nodes[x], isolates, temp_dir, output_dir)
        for x in tqdm(core_genes)
    )
    # remove single sequence files
    unaligned_sequence_files = filter(None, unaligned_sequence_files)
    # Get alignment commands
    commands = [
        get_alignment_commands(fastafile, output_dir, aligner, threads)
        for fastafile in unaligned_sequence_files
    ]
    # Run alignment commands
    multi_align_sequences(
        commands, output_dir + "aligned_gene_sequences/", threads, aligner
    )
    # Concatenate them together to produce the two output files
    concatenate_core_genome_alignments(core_gene_names, output_dir)
    return


def generate_summary_stats(output_dir):
    with open(output_dir + "gene_presence_absence_roary.csv", "r") as inhandle:
        gene_presence_absence = inhandle.read().splitlines()[1:]
    noSamples = len(gene_presence_absence[0].split(",")) - 14
    # Layout categories
    noCore = 0
    noSoftCore = 0
    noShell = 0
    noCloud = 0
    total_genes = 0
    # Iterate through GPA and summarise
    for gene in gene_presence_absence:
        proportion_present = float(gene.split(",")[3]) / noSamples * 100.0
        if proportion_present >= 99:
            noCore += 1
        elif proportion_present >= 95:
            noSoftCore += 1
        elif proportion_present >= 15:
            noShell += 1
        else:
            noCloud += 1
        total_genes += 1

    # write output
    with open(output_dir + "summary_statistics.txt", "w") as outfile:
        output = (
            "Core genes\t(99% <= strains <= 100%)\t"
            + str(noCore)
            + "\n"
            + "Soft core genes\t(95% <= strains < 99%)\t"
            + str(noSoftCore)
            + "\n"
            + "Shell genes\t(15% <= strains < 95%)\t"
            + str(noShell)
            + "\n"
            + "Cloud genes\t(0% <= strains < 15%)\t"
            + str(noCloud)
            + "\n"
            + "Total genes\t(0% <= strains <= 100%)\t"
            + str(total_genes)
        )
        outfile.write(output)

    return True
