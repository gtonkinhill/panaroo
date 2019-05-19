#Takes .gff output from prokka and outputs combined gene/protien sequences from each isolate
import os
from pybedtools import BedTool
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate, Seq
from io import StringIO
import tempfile

@profile
def get_gene_sequences(gff_file, file_number, temp_dir,
    prot_handle, dna_handle, csv_handle):
    #Get name and separate the prokka GFF into separate GFF and FASTA files

    # read in and write output for pybedtools
    lines = gff_file.read()
    split = lines.split('##FASTA')

    temp_gff_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)
    temp_gff_file.close()
    temp_fasta_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)
    temp_fasta_file.close()
    with open(temp_gff_file.name, 'w') as outfile:
        outfile.write(split[0].strip())
    with open(temp_fasta_file.name, 'w') as outfile:
        outfile.write(split[1].strip())

    # run pybedtools
    annotations = BedTool(temp_gff_file.name)
    annotations.sort()
    temp = annotations.sequence(fi=temp_fasta_file.name, s=True)

    # run through annotations/fasta and create output
    gene_count = 0
    contig_count = 0
    prev_contig = None
    dna_out = []
    protein_out = []
    annotation_out = []

    for ann, seq in zip(annotations, SeqIO.parse(annotations.seqfn, 'fasta')):
        
        # only consider CDS for now
        if "CDS" not in ann.fields[2]: continue

        # generate internal gene id
        if (ann.chrom!=prev_contig):
            if prev_contig is not None:
                gene_count = 0
                contig_count += 1
        clustering_id = "_".join([str(file_number), str(contig_count),
            str(gene_count)])

        # dna
        seq.id = clustering_id
        seq.description = ""
        dna_out.append(seq)

        # protein
        # prot = seq.translate()
        # prot.id = seq.id
        # prot.description = ""
        prot = SeqRecord(
            translate(seq.seq),
            id=seq.id,
            description = ""
        )
        if prot.seq[-1] == "*":
            prot.seq = prot.seq[0:-1]
        if "*" in prot.seq:
            print(prot)
            raise ValueError("Premature stop codon in a gene!")
        protein_out.append(prot)

        # annotations
        try:
            gene_name = ann.attrs['gene']
        except KeyError:
            gene_name = ""

        if gene_name == "":
            try:
                gene_name = ann.attrs['Name']
            except KeyError:
                gene_name = ""
        gene_description = " ".join(ann.fields).split("product=")[-1]

        annotation_out.append((
            os.path.splitext(os.path.basename(gff_file.name))[0], 
            ann.chrom,
            clustering_id, 
            ann.name,
            str(prot.seq),
            str(seq.seq),
            gene_name,
            gene_description
        ))

        prev_contig = ann.chrom
        gene_count += 1

    # write output
    SeqIO.write(dna_out, dna_handle, "fasta")
    SeqIO.write(protein_out, prot_handle, "fasta")
    for line in annotation_out:
        csv_handle.write(",".join(line) + '\n')
    
    # remove temporary files
    os.remove(temp_gff_file.name)
    os.remove(temp_fasta_file.name)
    annotations.delete_temporary_history(ask=False)

    return True


def process_prokka_input(gff_list, output_dir, temp_dir):
    try:
        protienHandle = open(output_dir + "combined_protein_CDS.fasta", 'w+')
        DNAhandle = open(output_dir + "combined_DNA_CDS.fasta", 'w+')
        csvHandle = open(output_dir + "gene_data.csv", 'w+')
        csvHandle.write(
            "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description\n"
        )
        for gff_no, gff in enumerate(gff_list):
            get_gene_sequences(gff, gff_no, temp_dir,
                protienHandle, DNAhandle, csvHandle)
        protienHandle.close()
        DNAhandle.close()
        csvHandle.close()
        return True
    except:
        print("Error reading prokka input!")
        raise RuntimeError("Error reading prokka input!")


if __name__ == "__main__":
    #used for debugging purpopses
    import sys
    thing = process_prokka_input([open(f, 'rU') for f in sys.argv[1:]], "", "")
