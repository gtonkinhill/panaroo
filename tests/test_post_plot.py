# test if the correct post analysis plots can be generated
from panaroo.__main__ import main
from panaroo.generate_abundance_plots import main as main_plots
import sys
import os.path
import tempfile


def test_post_plot(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        with tempfile.TemporaryDirectory() as tmpdir:
            # run panaroo
            sys.argv = ["", "-i", 
                datafolder + "aa1.gff", 
                datafolder + "aa2.gff", 
                "-o", tmpdir]
            main()

            # generate plots
            sys.argv = ["", "-i", 
                tmpdir + "/gene_presence_absence_roary.csv", 
                "-o", tmpoutdir,
                "--graph_type", "all"]
            main_plots()

        # check output
        assert os.path.isfile(tmpoutdir + "/jack1.png")
        assert os.path.isfile(tmpoutdir + "/jack2.png")
        assert os.path.isfile(tmpoutdir + "/ICE.png")
        assert os.path.isfile(tmpoutdir + "/chao2.png")
        assert os.path.isfile(tmpoutdir + "/acc.png")

    return


