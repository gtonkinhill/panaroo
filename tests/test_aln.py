# test if identical GFF3 files give the correct fasta alignment output
from panaroo.__main__ import main
import sys
import glob
from Bio import SeqIO


def test_aln(datafolder):

    # run panaroo
    sys.argv = ["", "-i", 
        datafolder + "aln.gff", 
        datafolder + "aln.gff", 
        "-o", datafolder,
        "-a", "pan"]
    main()

    # check alignments
    fasta_files = glob.glob(datafolder + "aligned_gene_sequences/*.fas")

    assert len(fasta_files)==9

    for f in fasta_files:
        assert len(list(SeqIO.parse(f, "fasta")))==2

    return