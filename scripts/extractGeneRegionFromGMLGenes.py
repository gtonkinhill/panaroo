#Takes a GML file from panaroo and an assembly and extracts the region surrounding a given pair of genes
#Will only extract the region in samples with a path between the genes
#To run: python extractGeneRegionFromGMLGenes.py -n .gmlFile -d gene_data.csv -g 2xGenesOfInterest -l LengthOfRegionToExtract -t AssemblyFileTranslation -o OutFile

import networkx as nx
import gffutils
from Bio import SeqIO
from io import StringIO
import argparse


def getGFFTranslation(
        translationFile
):  #Takes a file containing GFF file name and sample name and returns a dictionary with sample name keys and GFF file name objects
    gffDict = {}
    for line in translationFile:  #Iterate through the files
        gffDict[line.strip().split("\t")[1]] = line.strip().split("\t")[0]
    return gffDict


def clean_gff_string(
        gff_string):  #Removes other "##" starting lines from a gff file
    splitlines = gff_string.splitlines()
    lines_to_delete = []
    for index in range(len(splitlines)):
        if '##sequence-region' in splitlines[index]:
            lines_to_delete.append(index)
    for index in sorted(lines_to_delete, reverse=True):
        del splitlines[index]
    cleaned_gff = "\n".join(splitlines)
    return cleaned_gff


def get_gene_coordinates(
        gff_file, contig, genes, gene1
):  #Takes a GFF file, splits into CDS and FASTA and identifies the coordinates of the 2 genes
    gff = open(gff_file).read()  #Import the GFF file
    split = gff.split("##FASTA")

    geneCoordinates = []  #Will be filled with the coordinates of the 2 genes

    with StringIO(split[1]
                  ) as temp_fasta:  #Import the sequences from the GFF as fasta
        sequences = list(SeqIO.parse(temp_fasta, "fasta"))

    parsed_gff = gffutils.create_db(clean_gff_string(split[0]),
                                    dbfn=":memory:",
                                    force=True,
                                    keep_order=True,
                                    from_string=True)

    for geneSample in parsed_gff.all_features(
            featuretype=()):  #Iterate through the genes
        if geneSample.id == genes[0] or geneSample.id == genes[
                1]:  #Check if the gene is one of the genes of interest
            geneCoordinates.append(geneSample.start)
            geneCoordinates.append(geneSample.stop)
        if geneSample.id == gene1:  #Check if the gene is gene 1
            strand = geneSample.strand

    firstPosition = min(
        geneCoordinates
    ) - 1001  #The first position to be extracted from the contig
    endPosition = max(
        geneCoordinates
    ) + 1000  #The end position to be extracted from the contig
    if firstPosition < 0:
        firstPosition = 0
    if endPosition > len(sequences[int(contig)].seq):
        endPosition = len(sequences[int(contig)].seq)

    if strand == "+":  #Check if the sequence needs to be reverse complemented
        extractSequence = sequences[int(contig)].seq[
            firstPosition:
            endPosition]  #Assign to the sequence in the region of interest
    else:
        extractSequence = sequences[int(
            contig)].seq[firstPosition:endPosition].reverse_complement()

    return extractSequence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", help=".gml output file from panaroo")
    parser.add_argument("-d", help="gene_data.csv file from panaroo")
    parser.add_argument(
        "-g",
        nargs=2,
        help=
        "Genes of interest, the region surrounding these genes will be extracted"
    )
    parser.add_argument(
        "-l",
        help=
        "The length of the region on either side of the genes of interest that will be extracted, default = 1000 nucleotides",
        default="1000")
    parser.add_argument(
        "-t",
        help=
        "Translation file that maps the names of the GFF files onto the sample names in panaroo, tab delimited, e.g. sample1_assembly.fasta TAB 0"
    )
    parser.add_argument("-o", help="Output file")
    args = parser.parse_args()

    geneGraph = nx.read_gml(args.n)  #Import the gene graph
    geneData = open(args.d).readlines()
    gene1 = args.g[0]  #The name of the first gene
    gene2 = args.g[1]  #The name of the second gene

    gffTranslation = open(args.t).readlines(
    )  #Import the GFF files and their translated sample names
    gffFiles = getGFFTranslation(gffTranslation)

    outFile = open(args.o, "w")

    geneSequences = []  #Will be filled with the gene identifiers of interest
    gene1Sequences = [
    ]  #Will be filled with the ids for gene 1 for strand determination

    for node in geneGraph.nodes():  #Iterate through the nodes in the graph
        if geneGraph.node[node]["name"] == gene1 or geneGraph.node[node][
                "name"] == gene2:
            geneSequences.append(geneGraph.node[node]["seqIDs"])
        if geneGraph.node[node]["name"] == gene1:
            gene1Sequences.append(geneGraph.node[node]["seqIDs"])

    geneSeqIDs = [item for sublist in geneSequences
                  for item in sublist]  #Flatten the seqIDs into a single list

    geneNames = {
    }  #Will be filled with the names of the genes used in the GFF file in each sample
    contigNames = {
    }  #Will be filled with the contig the genes of interest are on in the GFF file in each sample
    gene1Names = {
    }  #Will be filled with the name of gene 1 in the GFF file in each sample

    for gene in geneData:  #Iterate through the genes
        for geneName in geneSeqIDs:  #Iterate through the gene names
            if "," + geneName + "," in gene:  #Check if the gene name is in the gene
                if geneName.split(
                        "_"
                )[0] in geneNames:  #Check if the sample already has a gene
                    geneNames[geneName.split("_")[0]].append(
                        gene.strip().split(",")[3])
                    contigNames[geneName.split("_")[0]].append(
                        geneName.split("_")[1])
                else:
                    geneNames[geneName.split("_")[0]] = [
                        gene.strip().split(",")[3]
                    ]
                    contigNames[geneName.split("_")[0]] = [
                        geneName.split("_")[1]
                    ]
                if geneName in gene1Sequences[0]:  #Check if the gene is gene 1
                    gene1Names[geneName.split("_")[0]] = gene.strip().split(
                        ",")[3]

    for sample in geneNames:  #Iterate through the samples
        if len(geneNames[sample]) == 2:  #Check if the sample has the 2 genes
            if contigNames[sample][0] == contigNames[sample][
                    1]:  #Check if there is a path between the 2 genes
                outFile.write(">" + sample + "\n")
                outFile.write(
                    str(
                        get_gene_coordinates(
                            gffFiles[sample], contigNames[sample][0],
                            geneNames[sample], gene1Names[sample])) + "\n")
            else:
                print("The genes in sample " + str(sample) +
                      " are on different contigs")
        else:
            print("Sample " + str(sample) + " does not contain the 2 genes")

    outFile.close()

#    for gene in geneData: #Iterate through the genes
#        for geneName1 in gene1SeqIDs: #Iterate through the names of gene 1
#            if "," + geneName1 + "," in gene: #Check if the gene name is in the gene
#                gene1Names[geneName1.split("_")[0]] = gene.strip().split(",")[3]
#        for geneName2 in gene2SeqIDs: #Iterate through the names of gene 2
#            if "," + geneName2 + "," in gene: #Check if the gene name is in the gene
#                gene2Names[geneName2.split("_")[0]] = gene.strip().split(",")[3]