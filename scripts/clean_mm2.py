import os
import sys
import assists
import logging
import taxidTools
import cleaners
import pandas as pd
import acc2taxid as at


def mm2_cleanup(database, taxdump, input_file, entrez_cred):
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

    # Now assign each required file using the expected labels
    fasta_file = input_paths["combined_lr_contigs_fa"]
    covstats_file = input_paths["cov_stats"]
    minimap_contigs = input_paths["minimap2_contig_out"]
    minimap_lr = input_paths["minimap2_reads_out"]

    check_files = [
        fasta_file,
        covstats_file,
        #minimap_contigs,
        minimap_lr,
    ]

    for f in check_files:
        assists.check_files(f)

        # clean up the fasta file so that its now a dataframe column
    headers, lengths = cleaners.convert_fasta_to_list(fasta_file)
    fasta_df = pd.DataFrame({"Fasta_Headers": headers, "Seq_Length": lengths})
    covstats_df = cleaners.extract_name_numreads(covstats_file)
    fasta_df = fasta_df.merge(
        covstats_df, left_on="Fasta_Headers", right_on="#rname", how='left'
    )
    fasta_df['numreads'] = fasta_df['numreads'].fillna("1").astype(int)
    logging.info("Created Fasta Headers")

    logging.info("Starting to clean up Minimap2 output")
    minimap_contigs_df = pd.DataFrame(columns = ["MM2_NT", "accession"])
    full_mm2_contigs_df = pd.DataFrame(columns = ["MM2_NT", "accession", "MM2_taxid", "pident", "alnlen", "evalue", "bitscore"])
    if os.path.isfile(minimap_contigs) is True and os.path.getsize(os.path.realpath(minimap_contigs)) != 0:
        minimap_contigs_df, full_mm2_contigs_df = cleaners.minimap_cleanup(minimap_contigs, "MM2_NT", taxdump, database, dirpath, entrez_cred)
    full_mm2_contigs_df = full_mm2_contigs_df.rename(columns={"MM2_NT": "Fasta_Headers"})
    minimap_lr_df = pd.DataFrame(columns = ["MM2_NT", "accession"])
    full_mm2_reads_df = pd.DataFrame(columns = ["MM2_NT", "accession", "MM2_taxid", "pident", "alnlen", "evalue", "bitscore"])
    if os.path.isfile(minimap_lr) is True and os.path.getsize(os.path.realpath(minimap_lr)) != 0:
        minimap_lr_df, full_mm2_reads_df = cleaners.minimap_cleanup(minimap_lr, "MM2_NT", taxdump, database, dirpath, entrez_cred)
    full_mm2_reads_df = full_mm2_reads_df.rename(columns={"MM2_NT": "Fasta_Headers"})

    # Filter out empty or all-NA DataFrames before concatenation
    minimap_parts = [
        ("minimap_contigs_df", minimap_contigs_df),
        ("minimap_lr_df", minimap_lr_df),
    ]
    minimap_parts = [(k, df) for k, df in minimap_parts if not df.empty and not df.isna().all().all()]
    if minimap_parts:
        minimap_df = pd.concat([df for _, df in minimap_parts], keys=[k for k, _ in minimap_parts])
    else:
        minimap_df = pd.DataFrame(columns=["MM2_NT", "accession"])

    full_mm2_parts = [
        ("full_mm2_contigs_df", full_mm2_contigs_df),
        ("full_mm2_reads_df", full_mm2_reads_df),
    ]
    full_mm2_parts = [(k, df) for k, df in full_mm2_parts if not df.empty and not df.isna().all().all()]
    if full_mm2_parts:
        full_minimap_df = pd.concat([df for _, df in full_mm2_parts], keys=[k for k, _ in full_mm2_parts])
    else:
        full_minimap_df = pd.DataFrame(columns=["MM2_NT", "accession", "MM2_taxid", "pident", "alnlen", "evalue", "bitscore"])
    mm2_accession_list = minimap_df["accession"].dropna().to_list()
    mm2_accession_list = list(dict.fromkeys(mm2_accession_list))
    logging.info("Merged and cleaned up Minimap output and ready to assign taxonomy")
    
    at.nucl_accession_search(database, mm2_accession_list, dirpath, entrez_cred)

    mm2_data = dirpath + "/mm2_data.tsv"
    minimap_df.to_csv(mm2_data, sep="\t", index=None)
    full_mm2_data = dirpath + "/full_mm2_data.tsv"
    full_minimap_df.to_csv(full_mm2_data, sep="\t", index=None)
    logging.info("Exporting Minimap2 cleaned data")

if __name__ == "__main__":
    mm2_cleanup(sys.argv[1], sys.argv[2], sys.argv[3], [sys.argv[4], sys.argv[5]])