import networkx as nx
from collections import defaultdict
import numpy as np

def generate_roary_gene_presence_absence(G, file_names, dna_file,
    output_dir):

    # clean file names
    file_names = [f.replace(",", "") for f in file_names]

    # calculate node average, min and max lengths
    seq_lengths = {}
    for rec in SeqIO.parse(dna_file, "fasta"):
        seq_lengths[str(rec.id))] = len(rec.seq)
    node_lengths = defaultdict(list)
    for node in G.nodes():
        for seq in G.node[node][seqIDs]:
            node_lengths[node].append(seq_lengths[seq])


    # generate file
    with open(output_dir + "gene_presence_absence.csv", 'w') as outfile:
        header = ["Gene","Non-unique Gene name","Annotation","No. isolates",
            "No. sequences","Avg sequences per isolate","Genome Fragment",
            "Order within Fragment","Accessory Fragment",
            "Accessory Order with Fragment","QC","Min group size nuc",
            "Max group size nuc","Avg group size nuc"
        ] + file_names
        outfile.write(",".join(header) + "\n")

        # Iterate through coponents writing out to file
        used_gene_names = set()
        unique_id_count = 0
        frag = 0
        for component in nx.connected_components(G):
            frag+=1
            count=0
            for node in component:
                count+=1
                if G.node[node]['annotation'] not in used_gene_names:
                    entry = [G.node[node]['annotation']]
                    used_gene_names.add(G.node[node]['annotation'])
                else:
                    entry = ["group_" + str(unique_id_count)]
                    unique_id_count += 1
                entry.append(G.node[node]['annotation'])
                entry.append(G.node[node]['description'])
                entry.append(len(G.node[node]['seqIDs']))
                entry.append(len(set(G.node[node]['members'])))
                entry.append((1.0*len(G.node[node]['seqIDs']))/
                    len(set(G.node[node]['members'])))
                entry.append(frag)
                entry.append(count)
                entry.append(["", "", ""])
                entry.append(np.min(node_lengths[node]))
                entry.append(np.max(node_lengths[node]))
                entry.append(np.mean(node_lengths[node]))
                pres_abs = [""]*len(file_names)
                for seq in G.node[node]['seqIDs']:
                    sample_id = int(seq.split("_")[0])
                    pres_abs[sample_id] = seq
                entry.append(pres_abs)
                outfile.write(",".join(entry) + "\n")

    return
