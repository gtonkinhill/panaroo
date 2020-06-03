#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script downloads fasta and gff files in fasta format for accessions of interest. These are then joined and cleaned for direct input into panaroo.
"""

import os

from .Entrez import fasta_downloader, gff_downloader
from .clean_gffs import clean_gffs

def main():
    
    args = get_options()
    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    current_dir = os.getcwd() + '/'
    os.makedirs(str(current_dir + 'genomes'))
    os.makedirs(str(current_dir + 'gffs_raw'))
    os.makedirs(str(current_dir + args.output_dir))
    accessions = args.accessions.split(',')
    headers = fasta_downloader(args.email, list(accessions))#, args.complete)
        
    headers['gff_filenames'] = gff_downloader(args.email, accessions)
    
    headers.apply(lambda row: clean_gffs(row["genome_ids"], row["gff_filenames"], 'genomes', 'gffs_raw', args.output_dir), axis = 1) #Count the number of guide RNA off targets
    
    return

def get_options(): #options for downloading and cleaning
    
    import argparse

    description = 'Download and process assemblies for panaroo input'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo-retrieve')

    io_opts = parser.add_argument_group('Entrez')
    
    io_opts.add_argument("-s",
                        "--search_term",
                        dest="accessions",
                        required=True,
                        help='sequences to download, specify a species for all accessions or specific accessions of interest. Separate accessions by "," with no spaces',
                        type=str) #Specify the search term for entrez
     
   # io_opts.add_argument("--complete",
                        #  dest="complete",
                         # help="only retrieve complete genomes",
                         # action='complete',
                         # default=None) #Option to only retrieve complete genomes  
    
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory for cleaned gffs",
                        type=str) #Specify the search term for entrez

    io_opts.add_argument("-e",
                        "--email",
                        dest="email",
                        required=True,
                        help="specify email for entrez access",
                        type=str) #Specify the search term for entrez
    
    #io_opts.add_argument("--store",
                          #dest="store",
                         # help="store fastas in output directory",
                          #action='store',
                         # default=None) #Option to only retrieve complete genomes  
    
    args = parser.parse_args()
    
    return (args)

if __name__ == '__main__':
    
    main()
