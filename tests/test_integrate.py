# test if the panaroo merge function is working correctly
from panaroo.__main__ import main
from panaroo.integrate import main as integrate_main
import numpy as np
import sys
import os
import tempfile

def test_integrate(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        with tempfile.TemporaryDirectory() as tmpoutdirB:
            # run panaroo on pairs of gff3 files
            sys.argv = ["", "-i",
                datafolder + "aa1.gff",
                datafolder + "aa2.gff",
                "--no_clean_edges",
                "--clean-mode", "sensitive",
                "-o", tmpoutdir]
            main()

            # merge the result
            sys.argv = ["", 
                "-d", tmpoutdir,
                "-i", datafolder + "aa3.gff",
                "-o", tmpoutdirB]
            try:
                integrate_main()
            except SystemExit:
                pass

            assert os.path.isfile(tmpoutdirB + "/final_graph.gml")

            # read gene p/a file
            pa = np.genfromtxt(tmpoutdirB + "/gene_presence_absence.Rtab",
                delimiter="\t", skip_header=1)
            assert pa.shape[0]==5151
            assert np.sum(pa[:,1:])==15339.0

            # read struct p/a file
            pa = np.genfromtxt(tmpoutdirB + "/struct_presence_absence.Rtab",
                delimiter="\t", skip_header=1)
            assert pa.shape == (106, 4)
            assert np.sum(pa[:,1:])==212

    return
