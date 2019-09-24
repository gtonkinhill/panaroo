# test if the correct post analysis plots can be generated
from panaroo.__main__ import main
from panaroo.generate_abundance_plots import main as main_plots
import sys
import os.path
import tempfile


def test_post_plot(datafolder):

    with tempfile.TemporaryDirectory() as tmpdir:
        # run panaroo
        sys.argv = ["", "-i", 
            datafolder + "aa1.gff", 
            datafolder + "aa2.gff", 
            "-o", tmpdir]
        main()

        # generate plots
        sys.argv = ["", "-i", 
            tmpdir + "/gene_presence_absence.csv", 
            "-o", datafolder,
            "--graph_type", "all"]
        main_plots()

    # check output
    assert os.path.isfile(datafolder + "jack1.png")
    assert os.path.isfile(datafolder + "jack2.png")
    assert os.path.isfile(datafolder + "ICE.png")
    assert os.path.isfile(datafolder + "chao2.png")
    assert os.path.isfile(datafolder + "acc.png")

    return


