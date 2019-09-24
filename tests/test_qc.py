# test if the correct qc plots can be generated
from panaroo.generate_qc_plots import main
import sys
import os.path


def test_qc(datafolder):

    # run panaroo
    sys.argv = ["", "-i", 
        datafolder + "aa1.gff", 
        datafolder + "aa2.gff", 
        "-o", datafolder,
        "--graph_type", "all",
        "--ref_db", datafolder + "test_mash_ref.msh"]
    main()

    assert os.path.isfile(datafolder + "MDS_mash_plot.html")
    assert os.path.isfile(datafolder + "MDS_mash_plot.png")
    assert os.path.isfile(datafolder + "ncontigs_barplot.png")
    assert os.path.isfile(datafolder + "ncontigs_boxplot.html")
    assert os.path.isfile(datafolder + "ngenes_barplot.png")
    assert os.path.isfile(datafolder + "ngenes_boxplot.html")
    assert os.path.isfile(datafolder + "mash_contamination_barplot.html")

    with open(datafolder + "mash_contamination_hits.tab", 'r') as infile:
        assert len(infile.readlines())==2

    with open(datafolder + "mash_dist.txt", 'r') as infile:
        assert len(infile.readlines())==3

    return