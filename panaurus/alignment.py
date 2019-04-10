import networkx as nx
import numpy as np
from Bio.Align.Applications import PrankCommandLine
from Bio.Align.Applications import MafftCommandline

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

def align_gene(G)
