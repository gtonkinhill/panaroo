# test if the correct qc plots can be generated
from panaroo.generate_qc_plots import main
import sys
import os.path
import tempfile


def test_qc(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        # run panaroo
        sys.argv = [
            "", "-i", datafolder + "aa1.gff", datafolder + "aa2.gff", "-o",
            tmpoutdir, "--graph_type", "all", "--ref_db",
            datafolder + "test_mash_ref.msh"
        ]
        main()

        assert os.path.isfile(tmpoutdir + "/MDS_mash_plot.html")
        assert os.path.isfile(tmpoutdir + "/MDS_mash_plot.png")
        assert os.path.isfile(tmpoutdir + "/ncontigs_barplot.png")
        assert os.path.isfile(tmpoutdir + "/ncontigs_boxplot.html")
        assert os.path.isfile(tmpoutdir + "/ngenes_barplot.png")
        assert os.path.isfile(tmpoutdir + "/ngenes_boxplot.html")
        assert os.path.isfile(tmpoutdir + "/mash_contamination_barplot.html")

        with open(tmpoutdir + "/mash_contamination_hits.tab", 'r') as infile:
            assert len(infile.readlines()) == 2

        with open(tmpoutdir + "/mash_dist.txt", 'r') as infile:
            assert len(infile.readlines()) == 3

    return