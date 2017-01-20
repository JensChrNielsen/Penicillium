#!/usr/bin/env python
"""
Extracts PKSs NRPSs and other relevant metadata from genbank files generated as output by antismash
Usage: python genbankParse.py /path/to/genbankfiles/*
"""

import sys
import re
from Bio import SeqIO
import glob


def get_smcog():
    """returns TRUE/FALSE if smCOG exists 
    and returns smCOG information"""
    increment = 0
    while True:
        try:
            if feat.qualifiers["note"][increment][:5] =="smCOG": 
                return True, feat.qualifiers["note"][increment].split(":")[2].split("(Scor")[0].replace(" ","")
        except (IndexError, KeyError):
            return False,""
        increment += 1


def printFasta(seq):
    """print fasta nicely"""
    line_length = 80
    newseq = ''
    count = 0
    for i in seq:
        if count == line_length:
            newseq += "\n"
            count = 0
        newseq += i
        count += 1
    return newseq


def getDomains():
	"""Get domains in synthase"""
	domains = []
	try:
		feat.qualifiers["sec_met"]

		for sec_met in feat.qualifiers["sec_met"]:
			domain = re.search(r'NRPS/PKS Domain:\s(\w*)\s\((\d*)-(\d*)\)',sec_met)
			if domain:
				domain_type = domain.group(1)
				domains.append( domain_type )
			
		return domains
	
	except KeyError:
        
		return domains


# Get location of genbank files 
path = sys.argv[1]
genbank_files = glob.glob(path)

# Open output files
cluster_out = open("clusters.tsv",'w')
synthase_out = open("synthases.faa",'w')
architecture_out = open("architecture.txt","w")

#Initialization
#all_cluster_types = {"pks": [], "nrps": [], "terpene": [], "hybrid": [], "siderophore": [], "indole": [], "other": []}
all_cluster_types = {"pks": [], "nrps": [], "terpene": [], "hybrid": [], "indole": [], "other": []}
order = []
cluster_counts = {}
synthase_counter = 0
no_hybrids = 0
no_hybrids_PKSNRPS = 0

for genbank_file in genbank_files:

    genome = SeqIO.parse(genbank_file,'genbank')
    
    #Get species name
    species = genbank_file[genbank_file.rfind("/")+1:genbank_file.index("_")]
    print species,"\n","#"*len(species)
    order.append(species)

    [all_cluster_types[key].append(0) for key in all_cluster_types] 

    for scaffold in genome:

        cluster_flag = 0
        for feat in scaffold.features:

		    # With or without antSMASH class "other"
            if feat.type == "cluster":# and feat.qualifiers["product"][0] != "other":
                
				# Get cluster location
                min_cluster_loc = min(feat.location)
                max_cluster_loc = max(feat.location)
				 
				# Extract info
                cluster_number = feat.qualifiers["note"][0].replace(" ","_").replace(":","")
                cluster_type = feat.qualifiers["product"][0]
				
				# Write info
                cluster_out.write("#"+species+"\t"+cluster_number+"\tCluster type:"+cluster_type+"\n")
                
                # Get counts of different cluters types
                try:
                    cluster_counts[cluster_type] += 1
                except KeyError:
                    cluster_counts[cluster_type] = 1

                # Get counts of grouped cluster types
                try:
                    all_cluster_types[cluster_type][-1] += 1
                except KeyError:
                    if cluster_type[2:] == "pks":
                        all_cluster_types["pks"][-1] += 1
                    elif '-' in cluster_type:
                        all_cluster_types["hybrid"][-1] += 1
                        no_hybrids += 1
                        print no_hybrids
                        if "pks" in cluster_type or "nrps" in cluster_type:
                            no_hybrids_PKSNRPS += 1
                            print no_hybrids_PKSNRPS
                    else:
                        print cluster_type
                        all_cluster_types["other"][-1] += 1
                
                cluster_flag = 1
                
            if feat.type == "CDS" and cluster_flag == 1:
				
                # Get CDS location
                max_cds_loc = max(feat.location)
                min_cds_loc = min(feat.location)

				# If CDS is within cluster
                if max_cds_loc <= max_cluster_loc and min_cds_loc >= min_cluster_loc:
					
                    gene_id = feat.qualifiers["locus_tag"][0]

                    is_smcog , smcog = get_smcog()
                    
                    cluster_out.write(gene_id+"\t"+smcog+"\n")
                    
                    # If smCOG info is described for gene
                    if is_smcog:
                        
						# If synthase, extract info
                        if ( "t1pks" in cluster_type or "t3pks" in cluster_type or "nrps" in cluster_type ) and ( re.search(r'Beta-ketoacyl_synthase',smcog) or re.search(r'malonyl_CoA-acyl_carrier_protein_transacylase',smcog) or re.search(r'AMP-dependent_synthetase_and_ligase',smcog) ):

                            sequence = feat.qualifiers["translation"][0]
                            
                            # Write fasta of of synthase
                            synthase_counter += 1
                            synthase_out.write(">"+str(species)+"|"+str(gene_id)+"|"+str(cluster_type)+"|"+str(smcog)+"|"+str(cluster_number)+"|synthase_"+str(synthase_counter)+"\n")
                            synthase_out.write(printFasta(sequence)+"\n")
	
                                
                            # Extract information on domain architecture
                            domains = getDomains()
							
							# Write architecture to file
                            #if domains != []:
                            if re.search(r'Beta-ketoacyl_synthase',smcog):
                                architecture_out.write(">"+str(species)+"|"+str(gene_id)+"|"+str(cluster_type)+"|"+str(smcog)+"|"+str(cluster_number)+"|synthase_"+str(synthase_counter)+"\n")
                                architecture_out.write( ";".join(domains)+"\n" )
                        

                else:
                    cluster_flag = 0

    
SPECIES_ORDER_PRINT = "species"+";"+";".join( order )
print SPECIES_ORDER_PRINT
for sm_class in all_cluster_types:
    toPrint = sm_class+";"+";".join(str(x) for x in all_cluster_types[sm_class])
    print toPrint

print cluster_counts

cluster_out.close()
synthase_out.close()
architecture_out.close()
 
