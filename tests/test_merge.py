# test if the panaroo merge function is working correctly
from panaroo.__main__ import main
from panaroo.merge_graphs import main as merge_main
import numpy as np
import sys
import os
import tempfile


def test_merge(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        with tempfile.TemporaryDirectory() as tmpdirA:
            # run panaroo on pairs of gff3 files
            sys.argv = ["", "-i", 
                datafolder + "aa1.gff", 
                datafolder + "aa2.gff", 
                datafolder + "aa3.gff",
                datafolder + "aa4.gff",
                "-o", tmpdirA]
            main()
            # read gene p/a file
            pa = np.genfromtxt(tmpdirA + "/gene_presence_absence.Rtab",
                delimiter="\t", skip_header=1)

            assert pa.shape == (5109, 5)
            assert np.sum(pa[:,1:])==20424

            # merge the result
            sys.argv = ["", "-d", 
                tmpdirA, 
                tmpdirA, 
                "-o", tmpoutdir]
            merge_main()

        assert os.path.isfile(tmpoutdir + "/merged_final_graph.gml")

        # read gene p/a file
        pa = np.genfromtxt(tmpoutdir + "/gene_presence_absence.Rtab",
            delimiter="\t", skip_header=1)

        assert pa.shape == (5109, 9)
        assert np.sum(pa[:,1:])==40848

        # read struct p/a file
        pa = np.genfromtxt(tmpoutdir + "/struct_presence_absence.Rtab",
            delimiter="\t", skip_header=1)
        
        assert pa.shape == (134 , 9)
        assert np.sum(pa[:,1:])==536

    return