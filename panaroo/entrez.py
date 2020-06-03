#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to retrieve fasta and gff files
"""

def fasta_downloader(email, Accessions): #, complete):
    from Bio import Entrez
    import pandas as pd
    Entrez.email = email
    #Download the genomes in fasta format  
    search = " ".join(Accessions)
    
    handle = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
    genome_ids = handle['IdList']
    print (f'found {len(genome_ids)} ids')

    for genome_id in genome_ids:
        record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        filename = '{}.fasta'.format(genome_id)
        directory = 'genomes/' + filename
        print('Writing:{}'.format(filename))
        #sequences.append(record.read())
        with open(directory, 'w') as f:
            f.write(record.read())
    
    headers = pd.DataFrame(genome_ids, columns = ["genome_ids"])
    return headers
    
def gff_downloader(email, accessions):
    from Bio import Entrez
    import os
    import urllib.request
    import subprocess

    Entrez.email = email
    
    assembly_ids = []
    for x in accessions:
        handle = Entrez.read(Entrez.esearch(db="assembly", term=x))
        assembly_ids.append(handle['IdList'])
        
    print (f'found {len(assembly_ids)} ids')
            
    os.chdir('gffs_raw')
    
    files = []
    for assembly_id in assembly_ids:
            
        esummary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
        
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        files.append(label)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.gff.gz') ##https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
        file = f'{label}.gff.gz'
        print('Writing:{}'.format(file))
        urllib.request.urlretrieve(link, file)
        gunzip_command = 'gunzip ' + file
        subprocess.run(gunzip_command, shell = True)
    
    os.chdir("..")
    return files
        
