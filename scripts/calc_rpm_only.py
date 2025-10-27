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
    "avgseqlen",
    "numreads",
    "final_taxid",
    "accession",
    "pident",
    "alnlen",
    "evalue",	
    "bitscore",
    "estimated_pg",
    "rpm_sample",
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

def pg_calculation(summary_tr, summary_fr, spikein_input):
    ercc_reads = spikein_input['plus_reads'].sum() + spikein_input['minus_reads'].sum()
    if ercc_reads != 0:
        # the picograms of ERCC will always be 10
        ercc_pg = 10

        # so the final calculatoion is
        total_pg = ercc_pg/(ercc_reads/summary_tr)

        scale_factor = summary_fr/summary_tr
        scaled_pg = scale_factor * total_pg
    else:
        scaled_pg = 0

    return scaled_pg

def calc_total_reads(summary_input):
    total_reads = 0
    fullyqc_reads = 0
    ercc_reads = 0
    sequins_reads = 0

    with open(summary_input, 'r') as f:
        for line in f:
            if line.startswith('fastp_before_total_reads:'):
                total_reads = int(line.split(':', 1)[1].replace(',', '').strip())
            elif line.startswith('fullyqc_reads:'):
                fullyqc_reads = int(line.split(':', 1)[1].replace(',', '').strip())
            elif line.startswith('ercc:'):
                ercc_reads = int(line.split(':', 1)[1].replace(',', '').strip())
            elif line.startswith('sequins:'):
                sequins_reads = int(line.split(':', 1)[1].replace(',', '').strip())
    
    if ercc_reads > 0:
        pg_reads = ercc_reads
    elif sequins_reads > 0:
        pg_reads = sequins_reads
    else:
        pg_reads = 0
    return total_reads, fullyqc_reads, pg_reads

def calc_rpm_only(sample_input, summary_input, spikein_input):
    # dirpath
    dirpath = os.path.dirname(os.path.abspath(sample_input))

    # calculate the total reads
    sample_tr, sample_fr, sample_pg = calc_total_reads(summary_input)
    
    # calculate the scale factor for RPM
    taxon_sum = sample_tr/1000000
    
    # read in the zscore_input files.
    taxon_counts = pd.read_csv(sample_input, sep="\t", header=0)

    # calculate ERCC/Sequins RPM scale factor if applicable
    if sample_pg > 0:
        if os.path.isfile(os.path.abspath(spikein_input)) is True and os.stat(os.path.abspath(spikein_input)).st_size != 0:
            ercc_info = pd.read_csv(
                spikein_input, 
                sep="\t", 
                header = None, 
                names=["plus_reads", "minus_reads"],
                usecols=[6,7],
                dtype="int64"
                )
            scaled_pg = pg_calculation(sample_tr, sample_fr, ercc_info)
        else:
            scaled_pg = 0
        taxon_counts['scale_factor'] = taxon_counts['numreads']/sample_fr
        taxon_counts['estimated_pg'] = scaled_pg * taxon_counts['scale_factor']
        taxon_counts.to_csv(dirpath + "/pg_per_hit.csv", index=None)
        taxon_counts = taxon_counts.drop(columns='scale_factor')
    else:
        taxon_counts['estimated_pg'] = 0
    
    # calculate the RPM for sample and negative
    taxon_counts['rpm_sample'] = taxon_counts['numreads']/taxon_sum

    # export list with rpm only
    taxon_counts = taxon_counts[final_column_order]
    rpm_output = dirpath + "/rpm_only.tsv"
    taxon_counts.to_csv(rpm_output, sep="\t", index=None)

    # extract detected pathogens list
    known_pgs = dirpath + "/detected_pathogens.tsv"
    pathogen_db_list = pathogen_db['Species'].to_list() + pathogen_db['AltNames'].to_list()
    match_species = taxon_counts[taxon_counts['species'].isin(pathogen_db_list) | taxon_counts['taxon'].isin(pathogen_db_list)]
    match_species = match_species[(match_species['taxon'] != 'No Hit')] # remove the No hit because z-score of no hit is never 0
    match_species.to_csv(known_pgs, sep="\t", index=None)

if __name__ == "__main__":
    calc_rpm_only(sys.argv[1], sys.argv[2], sys.argv[3])