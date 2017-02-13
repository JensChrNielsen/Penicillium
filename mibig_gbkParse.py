#!/usr/bin/env python
"""
#Extracts genes from mibig genbank files if BGC is of fungal origin
Usage: pyhton gbkClusterParse.py /path/to/genbankfiles/*
path should contain the MIBiG database as genbank files
"""
def printFasta(seq):
    newseq = ''
    count = 0
    for i in seq:
        if count == 80:
            newseq += "\n"
            count = 0
        newseq += i
        count += 1
    return newseq

import sys
import re
from Bio import SeqIO
import glob

# Get location of genbank files 
path = sys.argv[1]
genbank_files = glob.glob(path)
synthase = 0
out = open("fungi_mibig.txt",'w')
for genbank_file in genbank_files:

    fh = open(genbank_file,'r')
    genbank_file = SeqIO.parse(fh,'genbank')

    for cluster in genbank_file:
        accession = cluster.id
        cluster_definition = cluster.description
        
        try:
            kingdom = cluster.annotations["taxonomy"][1]
            organism = cluster.annotations["taxonomy"][-1]
        except IndexError:
            pass
    
        if kingdom == "Fungi":
            for feat in cluster.features:    
            
                if feat.type == "source": 
                    organism = feat.qualifiers["organism"][0]
                
                if feat.type == "CDS":
                    synthase += 1
                    seq = feat.qualifiers["translation"][0]
                    
                    header = ">"+accession+";"+organism+";"+cluster_definition+";synthase_"+str(synthase)
                    header = header.replace(" ","_")
                    out.write(header+"\n")
                    out.write(printFasta(seq)+"\n")
                    
out.close()