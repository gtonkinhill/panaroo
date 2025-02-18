#import shutil
#import tempfile
import os
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq

#janky workaround to run this as a script 
try:
    from .isvalid import is_valid_folder
except ImportError as e: 
    from isvalid import is_valid_folder
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
                         type=lambda x: is_valid_folder(parser, x))
    
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
    if len(input_list) == 1:
        with open(input_list[0], 'r') as inhandle:
            input_list = inhandle.read().splitlines()
    if verbose == True:
        list_of_isolate_names = tqdm(list_of_isolate_names)
    for isolate_id in list_of_isolate_names:
        parsed_gff = {}
        #get the right file in the original processing order
        base_input = [os.path.basename(input_id) for input_id in input_list]
        base_input = [input_id.replace(".gff", "") for input_id in base_input]
        input_file = next((s[1] for s in zip(base_input, input_list) if isolate_id == s[0]))
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
        initial_split = gff_line.split("\t")
        if initial_split[2] == "CDS":
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
                elif "locus_tag=" in attribute:
                    parsed_gff_line["Locus_tag"] = attribute.split("=")[1]
            parsed_gff.append(parsed_gff_line)

    return parsed_gff
    
def process_refound_gene(refound_id, pangenome_id, parsed_gff, refound_seqs, 
                         output_dir, G):
    #Get the refound sequence
    refound_seq = refound_seqs[refound_id]

    #Set up checks to ensure it is found
    scaffold_id = refound_seqs[refound_id][0]
    start = refound_seqs[refound_id][1]
    stop = refound_seqs[refound_id][2]
    strand = refound_seqs[refound_id][3]

    #Get additional data required for GFF annotation
    combined_gene_name = G.nodes[pangenome_id]["annotation"]
    gene_name = "_".join(combined_gene_name.strip(";").split(";")) 
    panaroo_name = G.nodes[pangenome_id]["name"]
    
    if G.nodes[pangenome_id]["paralog"] == 1:
        has_paralog = "True"
    else:
        has_paralog = "False"
    gene_description = G.nodes[pangenome_id]["description"]    
    #Put it all together and return a line of the GFF body
    gff_attributes = ";".join(["ID="+refound_id,
                               "locus_tag="+refound_id,
                               "name="+gene_name,
                               "description="+gene_description,
                               "panaroo_gene_cluster="+panaroo_name,
                               "inference=panaroo refound gene",
                               "has_pangenome_paralog="+has_paralog])
    gff_line = [scaffold_id, "Panaroo_refound", "candidate_gene", str(start+1), str(stop), ".", 
                strand, ".", gff_attributes]
    return "\t".join(gff_line)

def create_new_gffs(isolate_index, parsed_gffs, pp_isolate_genes,
                    gene_name_dic, refound_seqs, outdir, gff_format, G):

    #set up variables, add original GFF3 header to output
    new_gff_body_lines = []
    new_gff_body_lines.append(parsed_gffs[isolate_index]["header"])
    parsed_original_gffbody = parse_gff_body(parsed_gffs[isolate_index]["body"])
    pangenome_isolate_genes = 0

    #Need to go through all the pangenome genes
    #If isolates have been entirely removed due to being contaminats
    if str(isolate_index) not in pp_isolate_genes.keys():
        return []
    for pangenome_gene in pp_isolate_genes[str(isolate_index)]:
        pangenome_isolate_genes += 1
        #And all the genes from this isolate with the pan-genome gene
        for gene in pp_isolate_genes[str(isolate_index)][pangenome_gene]:
            if "refound" in gene:
                #Deal with refound genes seperately, don't need original gff
                refound_line = process_refound_gene(gene, pangenome_gene, 
                                                    parsed_gffs[isolate_index], 
                                                    refound_seqs, outdir, G)
                new_gff_body_lines.append(refound_line)
            else:
                #Identify original annotation
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
                combined_gene_name = G.nodes[pangenome_gene]["annotation"]
                gene_name = "_".join(combined_gene_name.strip(";").split(";"))                
                panaroo_name = G.nodes[pangenome_gene]["name"]

                
                if gene_name == "":
                    gene_name = "No_name"
                if G.nodes[pangenome_gene]["paralog"] == 1:
                    has_paralog = "True"
                else:
                    has_paralog = "False"
                gene_description = G.nodes[pangenome_gene]["description"]

                new_attributes = ";".join(["ID="+original_gene_data["ID"],
                                          "name="+gene_name,
                                          "locus_tag="+original_gene_data.get("Locus_tag", str("Panaroo_refound")),
                                          "description="+gene_description,
                                          "pangenome_id="+str(pangenome_gene),
                                          "panaroo_ID="+gene,
                                          "panaroo_gene_cluster="+panaroo_name,
                                          "eC_number="+original_gene_data.get("eC_number", str(None)),
                                           "prepanaroo_inference="+original_gene_data.get("inference", 
                                                                                          "Unknown_inference"),
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

    # Sort the annotation body of the gff file based on the first coordinate
    new_gff_body_lines[1:] = sorted(new_gff_body_lines[1:], key = lambda x: (x.split('\t')[0], int(x.split('\t')[3])))
                
    if gff_format == "prokka":
        new_gff_body_lines.append("##FASTA")
        new_gff_body_lines.append(parsed_gffs[isolate_index]["fasta"].strip())

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

    # Load isolate names, initial ID to clustering ID mapping
    if args.verbose:
        print("Loading panaroo output and input gff files...")
    seen = set()
    isolate_names = []
    gene_names = {}
    refound_seqs = {}
    with open(args.output_dir + "gene_data.csv", 'r') as infile:
        next(infile)
        for line in infile:
            splitinfo = line.strip().split(",")
            iso = splitinfo[0]
            gene_names[splitinfo[2]] = splitinfo[3]
            clusterid = splitinfo[2]
            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)
            if "refound" in splitinfo[2]:
                loc, strand = splitinfo[-1].split(';')
                loc = loc.split(':')[1].split('-')
                strand = strand.split(':')[1]
                refound_seqs[clusterid] = (splitinfo[1], int(loc[0]), int(loc[1]), strand)

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
            if type(G.nodes[pangenome_gene_id]["seqIDs"]) == list:                
                for isolate_geneid in G.nodes[pangenome_gene_id]["seqIDs"]:
                    if isolate_geneid.split("_")[0] == isolate_id:
                        isolate_genes[isolate_id][pangenome_gene_id] = isolate_genes[isolate_id].get(
                            pangenome_gene_id, []) + [isolate_geneid]
            elif type(G.nodes[pangenome_gene_id]["seqIDs"]) == str:
                isolate_geneid = G.nodes[pangenome_gene_id]["seqIDs"]
                if isolate_geneid.split("_")[0] == isolate_id:
                        isolate_genes[isolate_id][pangenome_gene_id] = isolate_genes[isolate_id].get(
                            pangenome_gene_id, []) + [isolate_geneid]
            
    
    #create and output new GFF files, multithreaded
    if args.verbose:
        print("Creating new gff files...")
    new_gffs = Parallel(n_jobs=args.n_cpu, prefer="threads")(
        delayed(create_new_gffs)(x, parsed_gffs, isolate_genes,
                    gene_names, refound_seqs, args.output_dir, args.format, G) for x in 
                    tqdm(range(len(isolate_names))))
    
    if args.verbose:
        print("Writing output...")
    if not os.path.exists(args.output_dir+"postpanaroo_gffs"):
        os.mkdir(args.output_dir+"postpanaroo_gffs")
    for index in range(len(new_gffs)):   
        output_gff(isolate_names[index], new_gffs[index], args.output_dir)


    return


if __name__ == '__main__':
    main()
