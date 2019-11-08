# simple test to see if the refind alignment is returning the correct sequence forward/reverse
from panaroo.find_missing import search_dna
import sys
from Bio.Seq import translate, reverse_complement, Seq

def test_aln(datafolder):

    # set up search
    db_seq = ("CAGATCACCGGGTGCCCGGCCGCGCGCACCGCCTCCACCAGCGGCGGCAGGCGCTCGCCG" +
        "ACCTTCTGCGCGCCCATCCGCGCGATCAGCGTCAGGCGCCCCGGCTCGCGGCGCGGGTCG" +
        "AGGCGCTCGCAGAGCGCCAGCAACTGGTCGCGGCCGATCTCCGGACCGACCTTGCAGGCC" +
        "ACCGGGTTGAGCACCTCGGCCAGCAGCGCCACATGGGCGCCGTCGACCTGGCGGGTGCGC" +
        "TCGCCGATCCACGGCCAGTGGGTCGAACCGAGATAGACCCGGCGCTGCTCGTCCTCGCGC" +
        "AGCATCGACAGCTCGTAGTCGAGCAGCAGCATCTCGTGGCTGGTCCAGACCGGCGAGGCA" +
        "TTCGCCTCCTGCCCGGACGCGGCGTCCCAGCCCAGGTGGCGCATGATGTTGCGCGCCGCC" +
        "GCGTAGCCCTTGAGGATCCGCTGCGGATCGGCCCGGCGCTGCTCGGCATGGGCCTCGCGG" +
        "CCGTTGACCATGTCGCCGCGATAGACCGGCAGGGTCTGCTC")

    search_sequence = "CCGGCCGCGCGCACCGCCTCCACCAGCGGCGGCAG"
    prop_match = 0.5
    pairwise_id_thresh = 0.99

    # Test forward

    seq, loc = search_dna(db_seq, search_sequence, prop_match, pairwise_id_thresh, False)

    assert(search_sequence==seq)
    assert(search_sequence==db_seq[loc[0]:loc[1]])

    # Test reverse

    search_sequence_rev = str(Seq(search_sequence).reverse_complement())

    seq, loc = search_dna(db_seq, search_sequence_rev, prop_match, pairwise_id_thresh, False)

    assert(search_sequence_rev==seq)
    assert(search_sequence==db_seq[loc[0]:loc[1]])

    # Test at ends

    search_sequence = "TTTTCAGATCACCGGGTGC"

    seq, loc = search_dna(db_seq, search_sequence, prop_match, pairwise_id_thresh, False)

    assert("CAGATCACCGGGTGC"==seq)
    assert("CAGATCACCGGGTGC"==db_seq[loc[0]:loc[1]])

    return