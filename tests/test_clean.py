# test if identical GFF3 files give the correct presence/absence matrix
from panaroo.__main__ import main
import numpy as np
import sys


def test_clean(datafolder):

    # run panaroo
    sys.argv = ["", "-i", 
        datafolder + "paralog.gff", 
        datafolder + "paralog.gff", 
        datafolder + "paralog.gff", 
        "-o", datafolder]
    main()

    # read p/a file
    pa = np.genfromtxt(datafolder + "gene_presence_absence.Rtab",
        delimiter="\t", skip_header=1)

    assert pa.shape == (5782, 4)
    assert (pa.shape[0]*(pa.shape[1]-1))==np.sum(pa[:,1:])

    return