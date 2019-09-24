# test if messy GFF3 files give the correct presence/absence matrix
from panaroo.__main__ import main
import numpy as np
import sys


def test_messy(datafolder):

    # run panaroo
    sys.argv = ["", "-i", 
        datafolder + "aa1.gff", 
        datafolder + "aa2.gff", 
        datafolder + "aa3.gff", 
        datafolder + "aa4.gff", 
        "-o", datafolder]
    main()

    # read p/a file
    pa = np.genfromtxt(datafolder + "gene_presence_absence.Rtab",
        delimiter="\t", skip_header=1)

    assert pa.shape == (5108, 5)
    assert np.sum(pa[:,1:])==20421

    return