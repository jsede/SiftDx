import os
import sys
import pandas as pd
import numpy as np

pathogen_db = pd.read_csv(os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "databases/pathogen_list.csv"), header = 0)
contam_db = pd.read_csv(os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "databases/known_contaminants.csv"), header = 0)

synthetic_list = [
    "shuttle vector",
    "synthetic construct",
    "cloning vector",
    "expression vector",
    "transformation vector",
    "reporter vector",
    "homo sapiens",
    "unclassified",
    "root"
]

phylum_exclusion = [
    "chordata", # vertebraes
    "arthropoda", # bugs
    "cnidaria", # aquatic animals
    "mollusca", # molluscs
    "ctenophora", # marine invertebraes
    "placozoa", # marine blobs
    "porifera", # sea sponges
    "rotifera", # microscopic wheel animals
    "echinodermata" #starfish
    "annelida" # seaworms
    ]

final_column_order = [
    "taxon",
    "numreads",
    "final_taxid",
    "accession",
    "pident",
    "alnlen",
    "evalue",	
    "bitscore",
    "rpm_sample",
    "rpm_ctrl",
    "zscore",
    "rank",
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
]

def calc_total_reads(summary_input):
    with open(summary_input, 'r') as f:
        for line in f:
            if line.startswith('fastp_before_total_reads:'):
                total_reads = int(line.strip().split(':')[1].replace(',', ''))
                break
    return total_reads

def calc_rpm_only(sample_input, summary_input):
    # dirpath
    dirpath = os.path.dirname(os.path.abspath(sample_input))

    # calculate the total reads
    sample_tr = calc_total_reads(summary_input)
    
    # calculate the scale factor for RPM
    taxon_sum = sample_tr/1000000
    
    # read in the zscore_input files.
    taxon_counts = pd.read_csv(sample_input, sep="\t", header=0)

    # calculate the RPM for sample and negative
    taxon_counts['rpm_sample'] = taxon_counts['numreads']/taxon_sum

    # export list with rpm only
    rpm_output = dirpath + "/rpm_only.tsv"
    taxon_counts.to_csv(rpm_output, sep="\t", index=None)

    # extract detected pathogens list
    known_pgs = dirpath + "/detected_pathogens.tsv"
    pathogen_db_list = pathogen_db['Species'].to_list() + pathogen_db['AltNames'].to_list()
    match_species = taxon_counts[taxon_counts['species'].isin(pathogen_db_list) | taxon_counts['taxon'].isin(pathogen_db_list)]
    match_species = match_species[(match_species['taxon'] != 'No Hit')] # remove the No hit because z-score of no hit is never 0
    match_species.to_csv(known_pgs, sep="\t", index=None)

if __name__ == "__main__":
    calc_rpm_only(sys.argv[1], sys.argv[2])