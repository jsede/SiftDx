import os
import sys
import logging
import taxidTools
import assists
import cleaners
import pandas as pd
import acc2taxid as at
import taxid2species as tl
import decision_tree as dt
import species2lineage as sl


def merge_nt_and_nr(database, taxdump, input_file, entrez_cred):
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
    blast_sr_file = input_paths["nt_alignments_sr_file"]
    diamond_file = input_paths["nr_alignments_file"]
    minimap_contigs = input_paths["minimap2_contig_out"]
    minimap_lr = input_paths["minimap2_reads_out"]
    kraken_file = input_paths["kraken_pluspf_file"]

    check_files = [
        fasta_file,
        covstats_file,
        blast_sr_file,
        diamond_file,
        #minimap_contigs,
        minimap_lr,
        kraken_file,
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

    
    # clean up the blast_sr file by dropping duplicates
    logging.info("Starting to clean up BLAST output")
    blast_sr_df = pd.DataFrame(columns = ["NT", "accession", "BLAST_Species"])
    full_blast_df = pd.DataFrame(columns = ["NT", "accession", "BLAST_Species", "staxids", "pident", "alnlen", "evalue", "bitscore"])
    if os.path.isfile(blast_sr_file) is True and os.path.getsize(os.path.realpath(blast_sr_file)) != 0:
        blast_sr_df, full_blast_df = cleaners.blast_cleanup(blast_sr_file, "NT", taxdump, dirpath, entrez_cred)
    full_blast_df = full_blast_df.rename(columns={"NT": "Fasta_Headers"})
    logging.info("Cleaned up BLAST output and ready to merge")

    # clean up minimap, join them together and merge it with
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
    minimap_df = pd.concat([df for _, df in minimap_parts], keys=[k for k, _ in minimap_parts])

    full_mm2_parts = [
        ("full_mm2_contigs_df", full_mm2_contigs_df),
        ("full_mm2_reads_df", full_mm2_reads_df),
    ]
    full_mm2_parts = [(k, df) for k, df in full_mm2_parts if not df.empty and not df.isna().all().all()]
    full_minimap_df = pd.concat([df for _, df in full_mm2_parts], keys=[k for k, _ in full_mm2_parts])
    mm2_accession_list = minimap_df["accession"].dropna().to_list()
    mm2_accession_list = list(dict.fromkeys(mm2_accession_list))
    logging.info("Merged and cleaned up Minimap output and ready to assign taxonomy")
    
    at.nucl_accession_search(database, mm2_accession_list, dirpath, entrez_cred)
    
    # clean up diamond
    logging.info("Starting to clean up Diamond output")
    diamond_df = pd.DataFrame(columns = ["NR", "accession"])
    full_diamond_df = pd.DataFrame(columns = ["NR", "sscinames", "staxids", "pident", "qlen", "evalue", "bitscore"])
    if os.path.isfile(diamond_file) is True and os.path.getsize(os.path.realpath(diamond_file)) != 0:
        diamond_df, full_diamond_df = cleaners.diamond_cleanup(diamond_file, "NR", taxdump, database, dirpath, entrez_cred) # filtering the percentage identity to above 75.
    full_diamond_df = full_diamond_df.rename(columns={"NR": "Fasta_Headers"})
    diamond_acc_list = diamond_df["accession"].dropna().to_list()
    diamond_acc_list = list(dict.fromkeys(diamond_acc_list))
    logging.info("Cleaned up Diamond output and ready to assign taxonomy")

    logging.info("Searching prot databases for taxids")
    at.prot_accession_search(database, diamond_acc_list, dirpath, entrez_cred)

    kraken_df = pd.read_csv(
        kraken_file,
        sep="\t",
        header=None,
        names=["ID", "Kraken_Species"],
        usecols=[1, 2],
    )
    kraken_df["Kraken_taxid"] = kraken_df["Kraken_Species"].apply(cleaners.get_taxid).astype(int)
    kraken_df["Kraken_Species"] = kraken_df["Kraken_Species"].apply(
        cleaners.remove_brackets
    )
    desired_cols = ["ID", "Kraken_Species", "Kraken_taxid"]
    kraken_df = kraken_df.reindex(columns=desired_cols)
    logging.info("Cleaned up Kraken output and ready merge")
    full_kraken_df = kraken_df.rename(columns={"ID": "Fasta_Headers", "Kraken_Species": "sscinames", "Kraken_taxid": "staxids"})
    full_kraken_df['staxids'] = pd.to_numeric(full_kraken_df['staxids'], errors='coerce')
    full_kraken_df = full_kraken_df[
        (full_kraken_df['staxids'] != 0) &
        (full_kraken_df['sscinames'].str.lower() != 'unclassified')
    ]
    add_cols = ["accession", "pident", "alnlen", "evalue", "bitscore"]
    for col in add_cols:
        full_kraken_df[col] = "-"

    # merge kraken dataframe with fasta headers
    merged_df = fasta_df.merge(
        kraken_df, left_on="Fasta_Headers", right_on="ID", how="outer"
    )



    # clean up the kraken info, get rid of unclassified, root and "other sequences" and make it dashes.
    merged_df["Kraken_Species"] = (
        merged_df["Kraken_Species"]
        .replace(["unclassified ", "root ", "other sequences "], "-")
        .str.strip()
    )

    # get lineage information, for each species.
    tl.taxid2lineage(dirpath, taxdump, entrez_cred)
    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict = loaded_cache[3]
    species_dict = loaded_cache[4]
    
    # clean up taxid
    df_map = [
        (minimap_df, "MM2_Species", "MM2_taxid"),
        (full_minimap_df, "MM2_Species", "MM2_taxid"),
        (blast_sr_df, "BLAST_Species", "staxids"),
        (full_blast_df, "sscinames", "staxids"),
        (diamond_df, "NR_Species", "NR_taxid"),
        (full_diamond_df, "sscinames", "staxids")
    ]

    for df, species_col, taxid_col in df_map:
        df[species_col] = df['accession'].map(species_dict)
        df[taxid_col] = df['accession'].map(taxid_dict)
    full_minimap_df = full_minimap_df.rename(columns={"MM2_Species": "sscinames", "MM2_taxid": "staxids"})
    

    # merge blast short reads dataframe with fasta headers + kraken
    blast_sr_df = blast_sr_df[["NT", "BLAST_Species"]]
    merged_df = merged_df.merge(
        blast_sr_df, left_on="Fasta_Headers", right_on="NT", how="outer"
    )
    merged_df = merged_df.merge(
        minimap_df, left_on="Fasta_Headers", right_on="MM2_NT", how="outer"
    )

    merged_df.drop(columns=["accession"], inplace=True) # get rid of accession column.

    # merge diamond dataframe with fasta headers + blast + kraken + minimap
    merged_df = merged_df.merge(
        diamond_df, left_on="Fasta_Headers", right_on="NR", how="outer"
    )
    merged_df.drop_duplicates(subset="Fasta_Headers", keep="first", inplace=True) # remove duplicates in fasta header to make sure we dont have extra.
    merged_df = merged_df.fillna("-").drop(
        columns=["NT", "ID", "accession", "NR", "MM2_NT"]
    ) # drop the unnecessary columns from the merge and fill blank spaces with dashes.
    merged_df["NT_Species"] = merged_df.apply(
        lambda row: row["BLAST_Species"]
        if row["MM2_Species"] == "-"
        else row["MM2_Species"],
        axis=1,
    )
    merged_df.drop(columns=["BLAST_Species", "MM2_Species"], inplace=True) # drop more unncessary columns from merge.
    
    # clean up the columns again make sure everything is a string.
    col_list = ["NT_Species", "NR_Species"]
    for col in col_list:
        merged_df[col] = merged_df[col].apply(cleaners.remove_num_str)
    merged_df = merged_df[merged_df['Fasta_Headers'] != '-']

    logging.info(f"Kraken, NT & NR merge complete")
    sample_df, full_sample_df = dt.decision(merged_df, taxdump, dirpath, entrez_cred)
    dataframes = [full_kraken_df, full_blast_df, full_minimap_df, full_diamond_df]
    non_empty_dataframes = [df for df in dataframes if not df.empty and not df.isna().all().all()]
    all_hits_df = pd.concat(non_empty_dataframes)

    desired_cols = [
        'taxid', 'superkingdom', 'kingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species', 'subspecies'
    ]

    all_hits_df.sort_values(by=["Fasta_Headers", "evalue", "bitscore"], ascending=[True, True, False], inplace=True)
    all_hits_df = all_hits_df.drop_duplicates(subset="Fasta_Headers", keep="first")
    all_hits_list = all_hits_df["staxids"].dropna().to_list()
    final_taxid_list = full_sample_df[
        (full_sample_df["final_taxid"].notna()) & (full_sample_df["final_taxid"] != "-")
    ]["final_taxid"].to_list()
    all_hits_list = final_taxid_list + all_hits_list
    all_hits_list = list(dict.fromkeys(all_hits_list))
    cleaners.get_lineage(all_hits_list, taxdump, dirpath, entrez_cred)
    lineage_cache = assists.check_lineage_json(loaded_cache[10])
    dfs = [pd.DataFrame(entry, index=[0]) for entry in lineage_cache]
    lineage_df = pd.concat(dfs, ignore_index=True, sort=False)
    lineage_df['rank'] = lineage_df['no rank'].replace('no rank', 'root')
    if "superkingdom" not in lineage_df.columns:
        lineage_df["superkingdom"] = None  # start fresh
        if "domain" in lineage_df.columns:
            lineage_df["superkingdom"] = lineage_df["acellular root"]
        if "acellular root" in lineage_df.columns:
            lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["domain"])
    for col in desired_cols:
        if col not in lineage_df.columns:
            lineage_df[col] = ""
    lineage_df = lineage_df[desired_cols].dropna(how='all').drop_duplicates(subset='taxid')
    lineage_df['taxid'] = (
    lineage_df['taxid']
        .astype(str)
        .str.replace('.0', '', regex=False)
    )
    full_sample_df['final_taxid'] = (
        full_sample_df['final_taxid']
        .astype(str)
        .str.replace('.0', '', regex=False)
    )
    full_sample_df = full_sample_df.merge(
        all_hits_df, left_on="Fasta_Headers", right_on="Fasta_Headers", how='left'
    )
    full_sample_df = full_sample_df.merge(
        lineage_df, left_on="final_taxid", right_on='taxid', how='left'
    )   
    full_sample_df = full_sample_df.drop(columns=["taxid","sscinames","staxids", "BLAST_Species"]).fillna("-")
    
    full_read_contig_info = dirpath + "/full_read_contig_info.tsv"
    full_sample_df.to_csv(full_read_contig_info, sep="\t", index=None)
    missing_lineage = full_sample_df[
        (full_sample_df['final_taxid'] != '-') & 
        (full_sample_df['final_taxid'] != '1') &
        (full_sample_df['superkingdom'] == '-')
    ]
    logging.info(f"Rows with missing lineage: {len(missing_lineage)}")
    collapse_df = cleaners.collapse_same_species(full_sample_df, taxdump)
    collapse_df['rank'] = collapse_df.apply(
        lambda row: sl.get_taxon_rank(row, 'final_taxid', taxdump),
        axis=1,
        result_type="expand"
    )
    zscore_input = dirpath + "/zscore_input.tsv"
    collapse_df.to_csv(zscore_input, sep="\t", index=None)
    
    return sample_df


if __name__ == "__main__":
    merge_nt_and_nr(sys.argv[1], sys.argv[2], sys.argv[3], [sys.argv[4], sys.argv[5]])