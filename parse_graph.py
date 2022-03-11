#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sam Lipworth
"""
import sys
import igraph as ig
import argparse
from Bio import SeqIO

def get_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    
    required.add_argument('-g', '--genes', action='store',
                          required=True,
                          help='List of genes to export')
    
    required.add_argument('-gml', '--graph', action='store',
                          required=True,
                          help='Panaroo final.gml graph'),
    required.add_argument('-f', '--fasta', action='store',
                          required=True,
                          help='Panaroo combined CDS file')
    required.add_argument('-o', '--outfile', action='store',
                          required=True,
                          help='Output file prefix')
    
    args = parser.parse_args(None if sys.argv[1:] else ('-h'))
    return args
    

def get_gene_ids(gene):

    args=get_arguments()
    print('Loading graph')
    g = ig.read(args.graph)
    print('Finding genes')
    query=g.vs.find(name=gene)
    query=query["geneIDs"]
    query=query.split(';')
    
    return(query)


def parse_fasta(gene,query):
    
    args=get_arguments()
    out_fasta = open(str(args.outfile +'_' + gene + '.'+ 'fasta'),'w')
    print('Loading panaroo fasta file')
    gene_fasta = SeqIO.parse(args.fasta,'fasta')
    print('Writing output to fasta files')
    for seq in gene_fasta:
    
        if seq.id in query: 
            
        
            SeqIO.write(seq,out_fasta, "fasta")
        
def main():
    args=get_arguments()
    
    
    with open(args.genes) as gene_file:
        for g in gene_file:
            g = g.strip('\n')
            print(g)
            queries=get_gene_ids(g)
            parse_fasta(g,queries)
    
    
if __name__ == "__main__":
    main()
    