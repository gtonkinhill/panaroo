from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

def output_sequence(node, temp_directory):
    sequences = []
    name = node["centroid"]
    raw_sequences =  ";".split(node["dna"])
    isolates = ";".split(node["members"])
    for i in range(len(raw_sequences)):
        sequences.append(SeqRecord(Seq(raw_sequences[i], generic_dna), id=isolates[i]))
    outname = temp_directory + '/' + name + ".fasta"
    SeqIO.write(sequences, outname, 'fasta')
    return outname

def get_alignment_commands(fastafile_name, aligner, no_threads):
    geneName = fastafile_name.split('.')[0]
    if aligner == "mafft":
        command = MafftCommandline(input=fastafile_name, auto=True, thread=no_threads, nuc=True)
    elif aligner == "prank":
        command = PrankCommandline(d=fastafile_name, o=geneName, f=8, codon=True)
    elif aligner == "clustal":
        command = ClustalOmegaCommandline(infile=fastafile_name, outfile=geneName+".aln.fas", seqtype="DNA", threads=no_threads)
    return command

def align_sequences(command, aligner):
    if aligner == "mafft":
        name = str(command).split("input=")[1].split().split('.')[0]
        stdout, stderr = command()
        with open(name+'.aln.fas', 'w+') as handle:
            handle.write(stdout)
    else:
        stdout, stderr = command()
    return True
