#import shutil
#import tempfile
import os
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq

#debugging purposes
import numpy as np
from time import time


import isvalid
#from .__init__ import __version__

def get_options():
    import argparse

    description = 'Generate GFF annotation files for each isolate from the pan-genome graph'
    parser = argparse.ArgumentParser(description=description,
                                     prog='generate_gff_files')

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory",
                         type=lambda x: isvalid.is_valid_folder(parser, x))
    
    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help=("input GFF3 files used to run panaroo. These are required " +
              "Can also take a file listing each gff file line by line."),
        type=str,
        nargs='+')

    # GFF Options
    core = parser.add_argument_group('GFFs')
    core.add_argument(
        "-f",
        "--f",
        dest="format",
        help=("Output format for the GFF files. 'prokka' with included" +
              " FASTA, or strict 'gff3'. Default: 'prokka'"),
        type=str,
        choices={'prokka', 'gff3'},
        default='prokka')


    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--verbose",
                        dest="verbose",
                        help="print additional output",
                        action='store_true',
                        default=False)
#    parser.add_argument('--version',
#                        action='version',
#                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)

def parse_all_gffs(list_of_isolate_names, input_list, verbose):
    ordered_parsed_gffs = []
    if verbose == True:
        list_of_isolate_names = tqdm(list_of_isolate_names)
    for isolate_id in list_of_isolate_names:
        parsed_gff = {}
        #get the right file in the original processing order
        input_file = next((s for s in input_list if isolate_id in s))
        raw_file = open(input_file, 'r')
        #read in and split off FASTA portion
        lines = raw_file.read().replace(',', '')
        split = lines.split('##FASTA')
        if len(split) > 1:
            parsed_gff["fasta"] = split[1]
        else:
            parsed_gff["fasta"] = None
        #get the header portion, remaining GFF
        header = []
        body = []
        for line in split[0].splitlines():
            if "##" in line:
                header.append(line)
            else:
                body.append(line)
        parsed_gff["header"] = "\n".join(header)
        parsed_gff["body"] = body
        ordered_parsed_gffs.append(parsed_gff)
    return ordered_parsed_gffs
        
def parse_gff_body(gff_body_list):
    parsed_gff = []
    for gff_line in gff_body_list:
        parsed_gff_line = {}
        initial_split = gff_line.split()
        parsed_gff_line["seqid"] = initial_split[0]
        parsed_gff_line["type"] = initial_split[2]
        parsed_gff_line["start"] = initial_split [3]
        parsed_gff_line["end"] = initial_split [4]
        parsed_gff_line["score"] = initial_split [5]
        parsed_gff_line["strand"] = initial_split [6]
        parsed_gff_line["phase"] = initial_split [7]
        attribute_split = " ".join(initial_split[8:]).split(";")
        for attribute in attribute_split:
            if "ID=" in attribute:
                parsed_gff_line["ID"] = attribute.split("=")[1]
            elif "eC_number=" in attribute:
                parsed_gff_line["eC_number"] = attribute.split("=")[1]
            elif "inference=" in attribute:
                parsed_gff_line["inference"] = attribute.split("=")[1]
        parsed_gff.append(parsed_gff_line)

    return parsed_gff
    
def process_refound_gene(refound_id, pangenome_id, parsed_gff, output_dir, G):
    #Get the refound sequence
    with open(output_dir + "gene_data.csv", 'r') as infile:
       next(infile)
       for line in infile:
           splitline = line.split(",")
           geneid = splitline[2]
           if refound_id == geneid:
               refound_seq = Seq(splitline[5])
    #Parse it so that there is easy reverse-translation
    isolate_raw_fasta =  StringIO(parsed_gff["fasta"]) 
    isolate_fasta = SeqIO.parse(isolate_raw_fasta, "fasta")
    #Set up checks to ensure it is found
    start = 0
    found = False
    #Find the refound sequence's start position, strand
    for seq_record in isolate_fasta:
        if refound_seq in seq_record.seq:
            start = seq_record.seq.index(refound_seq)
            strand = "+"
            scaffold_id = seq_record.id
            found = True
        elif refound_seq.reverse_complement() in seq_record.seq:
            start = seq_record.seq.index(refound_seq.reverse_complement())
            strand = "-"
            scaffold_id = seq_record.id
            found = True
    #Catch cases where it can't be found? Theoretically impossible
    if (start == 0) and (found ==False):
        print("Refound gene not found in  original scaffold!")
        print(refound_seq)
        print(parsed_gff["header"].split.readlines()[1])
        return None
    #Get additional data required for GFF annotation
    stop = start + len(refound_seq)
    gene_name = G.nodes[pangenome_id]["annotation"]
    if G.nodes[pangenome_id]["paralog"] == 1:
        has_paralog = "True"
    else:
        has_paralog = "False"
    gene_description = G.nodes[pangenome_id]["description"]    
    #Put it all together and return a line of the GFF body
    gff_attributes = ";".join(["ID="+refound_id,
                               "name="+gene_name,
                               "description="+gene_description,
                               "inference=Panaroo absent gene DNA re-finding",
                               "has_pangenome_paralog="+has_paralog])
    gff_line = [scaffold_id, "Panaroo_refound", ".", str(start+1), str(stop+1), ".", 
                strand, ".", gff_attributes]
    return "\t".join(gff_line)

def create_new_gffs(isolate_index, parsed_gffs, pp_isolate_genes,
                    gene_name_dic, outdir, gff_format, G):

    #set up variables, add original GFF3 header to output
    new_gff_body_lines = []
    new_gff_body_lines.append(parsed_gffs[isolate_index]["header"])
    parsed_original_gffbody = parse_gff_body(parsed_gffs[isolate_index]["body"])
    pangenome_isolate_genes = 0
    global_loop_start = time()
    #Need to go through all the pangenome genes
    for pangenome_gene in pp_isolate_genes[str(isolate_index)]:
        pangenome_isolate_genes += 1
        gene_times = []
        refound_times = []
        nonrefound_times = []
        #And all the genes from this isolate with the pan-genome gene
        for gene in pp_isolate_genes[str(isolate_index)][pangenome_gene]:
            gene_loop_start = time()
            if "refound" in gene:
                refound_start = time()
                #Deal with refound genes seperately, don't need original gff
                refound_line = process_refound_gene(gene, pangenome_gene, parsed_gffs[isolate_index], outdir, G)
                new_gff_body_lines.append(refound_line)
                refound_end = time()
                refound_times.append(refound_end-refound_start)
            else:
                #Identify original annotation
                nonrefound_start = time()
                original_gene_data = list(filter(lambda isolategff: 
                                  isolategff['ID'] == gene_name_dic[gene],
                                  parsed_original_gffbody))
                #Check that this is 1:1
                if len(original_gene_data) > 1:
                    print("Uh Oh, too much data")
                    print(original_gene_data)

                else:
                    original_gene_data = original_gene_data[0]
                #Get various other metadata for gene required for GFF3
                gene_name = G.nodes[pangenome_gene]["annotation"]
                if gene_name == "":
                    gene_name = "No_name"
                if G.nodes[pangenome_gene]["paralog"] == 1:
                    has_paralog = "True"
                else:
                    has_paralog = "False"
                gene_description = G.nodes[pangenome_gene]["description"]

                new_attributes = ";".join(["ID="+gene,
                                          "name="+gene_name,
                                          "description="+gene_description,
                                          "prepanaroo_ID="+original_gene_data["ID"],
                                          "eC_number="+original_gene_data.get("eC_number", str(None)),
                                          "prepanaroo_inference="+original_gene_data["inference"],
                                          "has_pangenome_paralog="+has_paralog])
                new_gene_line = "\t".join([original_gene_data["seqid"], 
                                      "Panaroo", 
                                      original_gene_data["type"],
                                      original_gene_data["start"],
                                      original_gene_data["end"],
                                      original_gene_data["score"],
                                      original_gene_data["strand"],
                                      original_gene_data["phase"],
                                      new_attributes])
                new_gff_body_lines.append(new_gene_line)
                nonrefound_end = time()
                nonrefound_times.append(nonrefound_end-nonrefound_start)
        gene_loop_end = time()
        gene_times.append(gene_loop_end-gene_loop_start)
                
    global_loop_end = time()
    global_iteration_length = {global_loop_end - global_loop_start}
    print("One iteration of an isolate took: " + str(global_iteration_length) + " seconds")
    print("Average per-gene time: " + str(np.mean(gene_times)))
    print("Average non-refound time: " + str(np.mean(nonrefound_times)))
    print("Average refound time: " + str(np.mean(refound_times)))
    if gff_format == "prokka":
        new_gff_body_lines.append("##FASTA")
        new_gff_body_lines.append(parsed_gffs[isolate_index]["fasta"])

    return new_gff_body_lines

def output_gff(isolate_name, gff_lines, outdir):
    gff = "\n".join(gff_lines)
    with open(outdir+"postpanaroo_gffs/"+isolate_name+"_panaroo.gff", "w") as outhandle:
        outhandle.write(gff)
    return True               
                
            

def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # Create temporary directory
    #temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # Load isolate names, initial ID to clustering ID mapping
    if args.verbose:
        print("Loading panaroo output and input gff files...")
    seen = set()
    isolate_names = []
    gene_names = {}
    with open(args.output_dir + "gene_data.csv", 'r') as infile:
        next(infile)
        for line in infile:
            splitinfo = line.split(",")
            iso = splitinfo[0]
            gene_names[splitinfo[2]] = splitinfo[3]
            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)

    # Load graph
    G = nx.read_gml(args.output_dir + "final_graph.gml")

    #parse input GFFs for headers, start/stop positions, FASTA
    parsed_gffs = parse_all_gffs(isolate_names, args.input_files, args.verbose)

    #Transform isolates-per-gene to localgeneids per pangenomeid per isolate
    if args.verbose:
        print("Cross-referencing input gffs and the pangenome...")

    isolate_genes = {}
    for pangenome_gene_id in tqdm(G.nodes):
        for isolate_id in G.nodes[pangenome_gene_id]["genomeIDs"].split(";"):
            isolate_genes[isolate_id] = isolate_genes.get(isolate_id, {})
            for isolate_geneid in G.nodes[pangenome_gene_id]["seqIDs"]:
                if isolate_geneid.split("_")[0] == isolate_id:
                    isolate_genes[isolate_id][pangenome_gene_id] = isolate_genes[isolate_id].get(
                        pangenome_gene_id, []) + [isolate_geneid]
            
    
    #create and output new GFF files, multithreaded
    if args.verbose:
        print("Creating new gff files...")
    new_gffs = Parallel(n_jobs=args.n_cpu, prefer="threads")(
        delayed(create_new_gffs)(x, parsed_gffs, isolate_genes,
                    gene_names, args.output_dir, args.format, G) for x in 
                    tqdm(range(len(isolate_names))))
    
    #temporarily single threaded, need to refactor the gff function
    #new_gffs = create_new_gffs(isolate_names, parsed_gffs, isolate_genes,
    #                gene_names, args.output_dir, args.format, G) 
    
    if args.verbose:
        print("Writing output...")
    if not os.path.exists(args.output_dir+"postpanaroo_gffs"):
        os.mkdir(args.output_dir+"postpanaroo_gffs")
    for index in range(len(new_gffs)):   
        output_gff(isolate_names[index], new_gffs[index], args.output_dir)

    # remove temporary directory
    #shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
