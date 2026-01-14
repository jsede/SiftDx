import os
import sys
import math
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
    "cellular root",
    "cellular organisms",
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

def set_zscore(row):
    for word in synthetic_list:
        if word in row['taxon'].lower():  # Convert both to lowercase for case-insensitive matching
            return 0
    for word in phylum_exclusion:
        if row['phylum'] is not np.nan:
            if word in row['phylum'].lower(): # remove chordata and arthropoda
                return 0
        else:
            return row['zscore']
    if row['kingdom'] is not np.nan:
        if "Viridiplantae" in row['kingdom']: # remove plant kingdom
            return 0
    return row['zscore']

def zscore_calculation(sample_input, negative_folder, summary_input, spikein_input):
    # dirpath
    dirpath = os.path.dirname(os.path.abspath(sample_input))


    # nav-ing to the summary file in the negative
    negative_input = os.path.join(os.path.abspath(negative_folder), "analysis/zscore_input.tsv")
    neg_summary_input = os.path.join(os.path.abspath(negative_folder), "summary/summary.txt")
    
    # calculate the total reads
    sample_tr, sample_fr, sample_pg = calc_total_reads(summary_input)
    negative_tr, negative_fr, negative_pg = calc_total_reads(summary_input)
    
    # calculate the scale factor for RPM
    taxon_sum = sample_tr/1000000
    neg_sum = negative_tr/1000000
    
    # read in the zscore_input files.
    taxon_counts = pd.read_csv(sample_input, sep="\t", header=0)
    neg_counts = pd.read_csv(negative_input, sep="\t", header=0)

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
        neg_counts['estimated_pg'] = 0

    # calculate the RPM for sample and negative
    taxon_counts['rpm_sample'] = taxon_counts['numreads']/taxon_sum
    neg_counts['rpm_ctrl'] = neg_counts['numreads']/neg_sum

    # calculate zscore for everything that it can (which is only the samples with both rpm_sample and rpm_ctrl)
    zscore_df = pd.merge(taxon_counts, neg_counts, on='taxon', how='inner')
    zscore_df = zscore_df[['taxon', 'rpm_sample', 'rpm_ctrl']]
    zscore_df['zscore'] = (zscore_df['rpm_sample'] - zscore_df['rpm_ctrl'].mean())/zscore_df['rpm_ctrl'].std()
    zscored = pd.merge(zscore_df, taxon_counts, on='taxon', how='inner').drop(columns="rpm_sample_y").rename(columns={'rpm_sample_x': 'rpm_sample'})
    if not zscored.empty:
        zscored['zscore'] = zscored.apply(set_zscore, axis=1)

        # need to check if the zscore calculations get rescaled if the numbers are too high.
        zscore_max = zscore_df['zscore'].max()
        zscore_min = zscore_df['zscore'].min()
        if zscore_max >= 99 and zscore_min != zscore_max:
            zscored['zscore'] = (
                (zscored['zscore'] - zscore_min) / (zscore_max - zscore_min) * (99 - 1) + 1
            )
            zscore_max = zscored['zscore'].max()
            zscore_min = zscored['zscore'].min()
        elif zscore_max >= 99 and zscore_min == zscore_max:
            zscored['zscore'] = 99
        
        # set zscore_max to 1 if <1
        if zscore_max <= 1:
            zscore_max = 1

    # grabbing only the taxon that are not in each other
    sample_only = taxon_counts[~taxon_counts['taxon'].isin(neg_counts['taxon'])].dropna(subset=['taxon']).reset_index(drop=True)
    sample_only['rpm_ctrl'] = np.nan
    neg_only = neg_counts[~neg_counts['taxon'].isin(taxon_counts['taxon'])].dropna(subset=['taxon']).reset_index(drop=True)
    neg_only['rpm_sample'] = np.nan
    
    # if taxon is not in negative, then score is immediately 100
    sample_only['zscore'] = 100

    # if taxon is in negative ONLY, then score is immediately -100
    neg_only['zscore'] = -100
    
    # clean up the basefile to be merged.
    sample_only = sample_only[final_column_order]
    neg_only = neg_only[final_column_order]
    zscored = zscored[final_column_order]

    # stack the files, and merge the final zscores in
    combined_counts = pd.concat([sample_only, neg_only, zscored], ignore_index=True)
    combined_counts['zscore'] = combined_counts.apply(set_zscore, axis=1)

    # separate the sample only & negative only 
    hundo_only = combined_counts[combined_counts['zscore'] == 100].copy()
    
    # Get RPM bounds from 100-zscore group
    hundo_rpm_max = hundo_only['rpm_sample'].max()
    hundo_rpm_min = hundo_only['rpm_sample'].min()
        
    # Rescale to be between zscore_max → 100
    if hundo_rpm_max != hundo_rpm_min:
        hundo_only['zscore'] = (
            (hundo_only['rpm_sample'] - hundo_rpm_min) /
            (hundo_rpm_max - hundo_rpm_min) *
            (100 - zscore_max) + zscore_max
        )
    else:
        hundo_only['zscore'] = 100
    
    # apply the new scores
    combined_counts.update(hundo_only)
    for col in combined_counts.select_dtypes(include='object').columns:
        combined_counts = combined_counts.fillna('-')
    zscore_output = dirpath + "/zscore.tsv"
    combined_counts.to_csv(zscore_output, sep="\t", index=None)

    # extract detected pathogens list
    known_pgs = dirpath + "/detected_pathogens.tsv"
    pathogen_db_list = pathogen_db['Species'].to_list() + pathogen_db['AltNames'].to_list()
    match_species = taxon_counts[taxon_counts['species'].isin(pathogen_db_list) | taxon_counts['taxon'].isin(pathogen_db_list)]
    match_species = match_species[(match_species['taxon'] != 'No Hit')] # remove the No hit because z-score of no hit is never 0
    match_species.to_csv(known_pgs, sep="\t", index=None)

    # set up for tpm krona chart
    krona_headers = {
        'taxon': '#queryID',
        'final_taxid': '#taxID',
        'zscore': '#score',
        'rpm': '#magnitude'
    }  
    krona_path = dirpath + "/krona_input.tsv"
    krona_csv = combined_counts[['taxon', 'final_taxid', 'zscore', 'rpm_sample']].copy()
    krona_csv.rename(columns=krona_headers, inplace=True)
    krona_csv.to_csv(krona_path, sep = '\t', index=False)

if __name__ == "__main__":
    zscore_calculation(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])