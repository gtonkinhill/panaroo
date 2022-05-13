# test if the panaroo merge function is working correctly
from panaroo.__main__ import main
from panaroo.integrate import main as integrate_main
import numpy as np
import sys
import os
import tempfile

def test_integrate(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir, \
        tempfile.TemporaryDirectory() as tmpdirA, \
        tempfile.TemporaryDirectory() as tmpdirB:
        # run panaroo on pairs of gff3 files
        sys.argv = ["", "-i",
            datafolder + "aa1.gff",
            datafolder + "aa2.gff",
            "--no_clean_edges",
            "--clean-mode", "sensitive",
            "-o", tmpdirA]
        main()

        # merge the result
        sys.argv = ["", 
            "-d", tmpdirA,
            "-i", datafolder + "aa3.gff",
            "-o", tmpoutdir]
        try:
            integrate_main()
        except SystemExit:
            pass

        assert os.path.isfile(tmpoutdir + "/final_graph.gml")

        # read gene p/a file
        pa = np.genfromtxt(tmpoutdir + "/gene_presence_absence.Rtab",
            delimiter="\t", skip_header=1)
        assert pa.shape[0]==5151
        assert np.sum(pa[:,1:])==15339.0

        # read struct p/a file
        pa = np.genfromtxt(tmpoutdir + "/struct_presence_absence.Rtab",
            delimiter="\t", skip_header=1)
        assert pa.shape == (106, 4)
        assert np.sum(pa[:,1:])==212

    return
