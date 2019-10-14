# test if messy GFF3 files give the correct presence/absence matrix
from panaroo.__main__ import main
import numpy as np
import sys
import tempfile


def test_messy(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        # run panaroo
        sys.argv = ["", "-i", 
            datafolder + "aa1.gff", 
            datafolder + "aa2.gff", 
            datafolder + "aa3.gff", 
            datafolder + "aa4.gff", 
            "-o", tmpoutdir]
        main()

        # read p/a file
        pa = np.genfromtxt(tmpoutdir + "/gene_presence_absence.Rtab",
            delimiter="\t", skip_header=1)

        assert pa.shape == (5110, 5)
        assert np.sum(pa[:,1:])==20424

    return