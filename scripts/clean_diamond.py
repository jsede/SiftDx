import os
import sys
import glob
import shutil
import assists
import logging
import taxidTools
import cleaners
from pathlib import Path
import pandas as pd
import acc2taxid as at

def diamond_cleanup(database, taxdump, input_file, mm2_data, entrez_cred):
    # Load taxdump for TaxidTools
    node = taxdump + "/nodes.dmp"
    lineage = taxdump + "/rankedlineage.dmp"
    logging.info("Loading Taxdump")
    taxdump = taxidTools.read_taxdump(node, lineage)

    # Assigning paths
    dirpath = os.path.dirname(os.path.abspath(input_file))
    input_paths = {}
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                key, path = parts
                input_paths[key] = os.path.join(dirpath, path)

    for json in Path(mm2_data).glob("*.json"):
        shutil.copy2(json, Path(dirpath) / json.name)

    # Now assign each required file using the expected labels
    fasta_file = input_paths["combined_lr_contigs_fa"]
    covstats_file = input_paths["cov_stats"]
    diamond_file = input_paths["nr_alignments_file"]
    
    check_files = [
        fasta_file,
        covstats_file,
        diamond_file
    ]

    for f in check_files:
        assists.check_files(f)

    # clean up the fasta file so that its now a dataframe column
    headers, lengths = cleaners.convert_fasta_to_list(fasta_file)
    fasta_df = pd.DataFrame({"Fasta_Headers": headers, "Seq_Length": lengths})
    covstats_df = cleaners.extract_name_numreads(covstats_file)
    if not fasta_df.empty or not covstats_df.empty:
        fasta_df = fasta_df.merge(
            covstats_df, left_on="Fasta_Headers", right_on="#rname", how='left'
        )
        fasta_df['numreads'] = fasta_df['numreads'].fillna("1").astype(int)
        logging.info("Created Fasta Headers")
    else:
        logging.warning("Fasta or Covstats DataFrame is empty, skipping merge.")
        diamond_df = pd.DataFrame(columns = ["NR", "accession"])
        full_diamond_df = pd.DataFrame(columns = ["NR", "sscinames", "staxids", "pident", "qlen", "evalue", "bitscore"])
        diamond_data = dirpath + "/diamond_data.tsv"
        diamond_df.to_csv(diamond_data, sep="\t", index=None)
        full_diamond_data = dirpath + "/full_diamond_data.tsv"
        full_diamond_df.to_csv(full_diamond_data, sep="\t", index=None)
        logging.info("Exporting Diamond cleaned data")
        sys.exit(0)
    fasta_df['numreads'] = fasta_df['numreads'].fillna("1").astype(int)
    logging.info("Created Fasta Headers")

    logging.info("Starting to clean up Diamond output")
    diamond_df = pd.DataFrame(columns = ["NR", "accession"])
    full_diamond_df = pd.DataFrame(columns = ["NR", "sscinames", "staxids", "pident", "qlen", "evalue", "bitscore"])
    if os.path.isfile(diamond_file) is True and os.path.getsize(os.path.realpath(diamond_file)) != 0:
        diamond_df, full_diamond_df = cleaners.diamond_cleanup(diamond_file, "NR", taxdump, dirpath, entrez_cred) # filtering the percentage identity to above 75.
    full_diamond_df = full_diamond_df.rename(columns={"NR": "Fasta_Headers"})
    diamond_acc_list = diamond_df["accession"].dropna().to_list()
    diamond_acc_list = list(dict.fromkeys(diamond_acc_list))
    logging.info("Cleaned up Diamond output and ready to assign taxonomy")

    logging.info("Searching prot databases for taxids")
    at.prot_accession_search(database, diamond_acc_list, dirpath, entrez_cred)

    
    diamond_data = dirpath + "/diamond_data.tsv"
    diamond_df.to_csv(diamond_data, sep="\t", index=None)
    full_diamond_data = dirpath + "/full_diamond_data.tsv"
    full_diamond_df.to_csv(full_diamond_data, sep="\t", index=None)
    logging.info("Exporting Diamond cleaned data")

if __name__ == "__main__":
    diamond_cleanup(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], [sys.argv[5], sys.argv[6]])