# Test: nodes whose centroids are different sequence IDs but
# identical sequences must not lose a centroid when they are merged.
#
# Multithreaded cd-hit (-T > 1) does not guarantee that identical sequences end up in
# the same cluster, so two clusters can have different centroid IDs with identical
# protein/DNA sequences.
#
# Calls the code in the order panaroo's main() uses: collapse mistranslations
# (DNA), then collapse gene families (protein) to try to catch this problem on a
# minimal example.
#
import random
import tempfile

import networkx as nx
from intbitset import intbitset

from panaroo.clean_network import collapse_families
from panaroo.merge_nodes import merge_node_cluster

CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in ("TAA", "TAG", "TGA")]


def random_gene(rng, n_codons=100):
    from Bio.Seq import Seq
    dna = "ATG" + "".join(rng.choice(CODONS) for _ in range(n_codons)) + "TAA"
    return dna, str(Seq(dna).translate(to_stop=True))


def add_gene_node(G, node_id, centroid, seq_ids, dna, protein):
    members = intbitset([int(s.split("_")[0]) for s in seq_ids])
    G.add_node(node_id,
               size=len(members),
               centroid=[centroid],
               maxLenId=0,
               members=members,
               seqIDs=set(seq_ids),
               hasEnd=False,
               protein=[protein],
               dna=[dna],
               annotation="",
               description="",
               lengths=[len(dna)] * len(seq_ids),
               longCentroidID=(len(dna), centroid),
               paralog=False,
               mergedDNA=False)


def add_edge(G, a, b):
    members = G.nodes[a]["members"] & G.nodes[b]["members"]
    G.add_edge(a, b, size=len(members), members=members)


def build_graph():
    """Hub gene H (genomes 0,1) with three neighbours.

    A and B have identical DNA and protein sequences in genome 0 and genome 1
    but are separate nodes with different centroid IDs, as happens when
    multithreaded cd-hit splits identical sequences into two clusters. C and D
    are differing genes that give H a degree > 2 to run collapse_families.
    """
    rng = random.Random(1)
    hub_dna, hub_prot = random_gene(rng)
    dup_dna, dup_prot = random_gene(rng)
    c_dna, c_prot = random_gene(rng)
    d_dna, d_prot = random_gene(rng)

    G = nx.Graph()
    add_gene_node(G, 1, "0_0_0", ["0_0_0", "1_0_0"], hub_dna, hub_prot)
    add_gene_node(G, 2, "0_0_1", ["0_0_1"], dup_dna, dup_prot)   # A
    add_gene_node(G, 3, "1_0_1", ["1_0_1"], dup_dna, dup_prot)   # B: same seq as A
    add_gene_node(G, 4, "0_0_2", ["0_0_2", "1_0_2"], c_dna, c_prot)
    add_gene_node(G, 5, "0_0_3", ["0_0_3", "1_0_3"], d_dna, d_prot)
    for n in (2, 3, 4, 5):
        add_edge(G, 1, n)

    seqid_to_centroid = {}
    for n in G.nodes():
        for sid in G.nodes[n]["seqIDs"]:
            seqid_to_centroid[sid] = G.nodes[n]["centroid"][0]
    return G, seqid_to_centroid


def assert_parallel(G):
    for n in G.nodes():
        nd = G.nodes[n]
        assert len(nd["centroid"]) == len(nd["protein"]) == len(nd["dna"]), (
            n, nd["centroid"], nd["protein"], nd["dna"])


def test_merge_node_cluster_keeps_centroids_and_sequences_parallel():
    G, _ = build_graph()
    dup_protein = G.nodes[2]["protein"][0]
    dup_dna = G.nodes[2]["dna"][0]
    G = merge_node_cluster(G, [2, 3], 99, multi_centroid=False)

    nd = G.nodes[99]
    assert sorted(nd["centroid"]) == ["0_0_1", "1_0_1"]
    assert_parallel(G)
    assert nd["protein"] == [dup_protein, dup_protein]
    assert nd["dna"] == [dup_dna, dup_dna]


def test_collapse_families_with_identical_centroid_sequences():
    G, seqid_to_centroid = build_graph()
    with tempfile.TemporaryDirectory() as tmp:
        G = collapse_families(G,
                              seqid_to_centroid=seqid_to_centroid,
                              outdir=tmp + "/",
                              dna_error_threshold=0.98,
                              correct_mistranslations=True,
                              quiet=True)[0]
        merged = [n for n in G.nodes() if "1_0_1" in G.nodes[n]["seqIDs"]]
        assert len(merged) == 1
        assert {"0_0_1", "1_0_1"} <= set(G.nodes[merged[0]]["centroid"])
        assert_parallel(G)

        G, distances, centroid_to_index = collapse_families(
            G,
            seqid_to_centroid=seqid_to_centroid,
            outdir=tmp + "/",
            family_threshold=0.7,
            correct_mistranslations=False,
            quiet=True)
    assert set(seqid_to_centroid.values()) <= set(centroid_to_index)
    assert_parallel(G)


if __name__ == "__main__":
    test_collapse_families_with_identical_centroid_sequences()
    print("OK")
