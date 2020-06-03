#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to remove non-gene annotations from NCBI GFF files  
"""

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


def clean_gffs(genome_id, gff_filename, input_fastas, input_gffs, output_dir):
    import re 
    
    gff = input_gffs + '/' + gff_filename + '.gff'
    o = open(gff, 'r')
    stored = str(o.read()) 
    
    d = {"Dbxref=": "dbxref=", "Notes=": "note="}
    stored = replace_all(stored, d) #Ensure GFF annotation format is consistent with prokka output
    
    gene_list = stored.splitlines() 
    title = gene_list[7]
    gene_list = gene_list[9:-1] #remove non-essential content
    
    genes = []
    
    for y in range(len(gene_list)):
        if re.split(r'\t+', gene_list[y])[2] == 'gene':
            genes.append(gene_list[y])
    
    cleaned_gffs = "\n".join(str(z) for z in genes)
    cleaned_gffs = title + '\n' + cleaned_gffs

    #Concatenate downloaded fasta and reformatted GFF
    fasta_filename = genome_id
    fasta_file = input_fastas +'/' + fasta_filename + '.fasta'
    cleaned_gffs = cleaned_gffs + '\n' + '##FASTA' + '\n' + str((open(fasta_file,'r').read())) 
    filename_cleaned = output_dir + '/' + genome_id + '_cleaned.gff'
    outfile = open(filename_cleaned,'w')          
    outfile.write(cleaned_gffs)
    outfile.close()

    return