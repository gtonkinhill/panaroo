from hmmer import run_jackhmmer
from Bio import SeqIO


def cluster_jackhmmer(input_file, output_dir, evalue, n_cpu):

    # temp file name
    temp_prot_file = output_dir + "/temp_protein_sequence.fasta"
    temp_query_file = output_dir + "/temp_query_sequence.fasta"
    temp_hmmer_file = output_dir + "/temp_hmmer_tblout.txt"
    temp_align_file = output_dir + "/temp_hmmer_align.txt"

    # read in input to dictionary and calculate sequence lengths
    seq_dict = {rec.id : (str(rec.seq), len(rec.seq)) for rec in SeqIO.parse(
        input_file, "fasta")}

    # iterative run jackhmmer using the longest sequence as a query and removing
    # matching  sequences
    # TODO: add in ability to take initial hmmer models (which are run prior)
    while(len(seq_dict)>1):
        # find longest sequence and write remaining to temp file
        with open(temp_prot_file, 'w') as outfile:
            max_length = -1
            for id in seq_dict:
                if seq_dict[id][1] > max_length:
                    if max_length>0:
                        outfile.write(">" + long_seq[0] + "\n" +
                            long_seq[1]+ "\n")
                    long_seq = (id, seq_dict[id][0])
                    max_length = seq_dict[id][1]

                else:
                    outfile.write(">" + id + "\n" +
                        seq_dict[id][0]+ "\n")

        # write query file (longest)
        with open(temp_query_file, 'w')  as outfile:
            outfile.write(">" + long_seq[0] + "\n" +
                long_seq[1]+ "\n")

        # run jackhmmer
        run_jackhmmer(temp_query_file, temp_prot_file,
            output_file=output_dir+"/temp_jackhmmer_out.txt",
            align_file=temp_align_file,
            E=evalue, incE=evalue,
            tblout=temp_hmmer_file,
            domtblout=output_dir+"/temp_jackhmmer_domtblout.txt",
            n_cpu=n_cpu)

        break
        # remove found sequences from seq_dict
        # with open(temp_hmmer_file, 'rU') as infile:
        #     for line in infile:
        #

    # return a list of hmmer models
    return []
