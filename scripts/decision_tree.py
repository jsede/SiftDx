import taxidTools
import logging
import pandas as pd
import cleaners
import entrez_search as es
import assists
import json

def decision(species_df, taxdump, dirpath, entrez_cred):
    """
    After all the contigs and unassembled reads have been annotated using the different software
    We need to create a final column that decides what the species is
    """
    # load cache files!
    loaded_cache = assists.load_cache_files(dirpath)

    # lets firstly check if the row contains "-"
    species_df["match"] = species_df.apply(lambda row: row[1:].nunique() == 1, axis=1)
    logging.info(f"Finding identical species matches...")

    # apply the determine final species function for every row.
    species_df[["final", "final_taxid"]] = species_df.apply(
        lambda row: determine_final_species(row, taxdump, entrez_cred),
        axis=1,
        result_type="expand"
    )
    logging.info(f"Determining best match species hit for unmatching species...")
    
    # extract the rows without final
    pd.set_option('mode.chained_assignment', None)
    logging.info("Assigning complete mismatches to their most common ancestor...")
    filtered_df = species_df[species_df['final'].isna() | (species_df['final'] == '')]
    if not filtered_df.empty:
        filtered_df.reset_index(inplace=True, drop=True)
        filtered_df[["final", "final_taxid"]] = filtered_df.apply(
            lambda row: taxid_assignment(row, taxdump),
            axis=1,
            result_type="expand"
        )
    
    stacked_df = pd.concat([species_df, filtered_df])
    stacked_df = stacked_df[stacked_df['final'].notna() & (stacked_df['final'] != '')]
    stacked_df.sort_values(by='Fasta_Headers', inplace=True, ascending=False)
    stacked_df['final_taxid'] = stacked_df["final_taxid"].replace("0", '-')
    stacked_df[["final", "final_taxid"]] = stacked_df.apply(
        lambda row: fix_taxid(row, taxdump),
        axis=1,
        result_type="expand"
    )
    full_stacked_df = stacked_df[['Fasta_Headers', 'Seq_Length', 'numreads', 'Kraken_Species', 'NT_Species', 'NR_Species', 'final', 'final_taxid']]
    full_stacked_df.to_csv(dirpath + "/final_decisions.tsv",sep="\t", index=None)
    stacked_df = stacked_df[['Fasta_Headers', 'Seq_Length', 'numreads', 'final', 'final_taxid']]

    
    # store the species + taxid into a json file.
    selected_columns = stacked_df[['final', 'final_taxid']]
    json_data = selected_columns.to_dict(orient='records')
    with open(loaded_cache[8], "w") as file:
        for item in json_data:
            json.dump({item["final"]: str(item['final_taxid']).replace('.0', '')}, file)
            file.write('\n')
    cleaners.rm_dupes_json(loaded_cache[8])
    return stacked_df, full_stacked_df


def determine_final_species(row, taxdump, entrez_cred):
    kraken_genus = cleaners.extract_genus(row, 'Kraken_taxid', 'Kraken_Species', taxdump, entrez_cred)
    nr_genus = cleaners.extract_genus(row, 'MM2_taxid', 'NT_Species', taxdump, entrez_cred)
    nt_genus = cleaners.extract_genus(row, 'NR_taxid', 'NR_Species', taxdump, entrez_cred)
    kraken_species = cleaners.extract_species(row["Kraken_Species"])
    nr_species = cleaners.extract_species(row["NR_Species"])
    nt_species = cleaners.extract_species(row["NT_Species"])

    primate_genus_list = [
        "Pan",
        "Gorilla",
        "Pongo",
        "Hylobates",
        "Symphalangus",
        "Nomascus",
        "Macaca",
        "Saimiri",
    ]
    
    synthetic_list = [
        "shuttle vector",
        "synthetic construct",
        "cloning vector",
        "expression vector",
        "transformation vector",
        "reporter vector",
    ]
    
    # If all rows match, then just print one of the columns (they all match anyway so it doesnt matter which we pick)
    if row["match"] is True:
        return row["Kraken_Species"], row["Kraken_taxid"]
    
    # Filtering the human reads
    elif (kraken_genus == "Homo" or row["Kraken_Species"] == "-") and \
        (nt_genus in primate_genus_list or
        nr_genus in primate_genus_list):
        return "Homo sapiens", 9606
        # If only one column contains valid data and the others are "-"
    
    elif (
        row["Kraken_Species"] != "-"
        and row["NR_Species"] == "-"
        and row["NT_Species"] == "-"
    ):
        return row["Kraken_Species"], row["Kraken_taxid"]
    elif (
        row["Kraken_Species"] == "-"
        and row["NR_Species"] != "-"
        and row["NT_Species"] == "-"
    ):
        return row["NR_Species"], row["NR_taxid"]
    elif (
        row["Kraken_Species"] == "-"
        and row["NR_Species"] == "-"
        and row["NT_Species"] != "-"
    ):
        return row["NT_Species"], row["MM2_taxid"]
    
    # If two columns match and the other column is not "-", then majority wins.
    elif row["Kraken_Species"] == row["NR_Species"] and row["Kraken_Species"] != "-":
        return row["Kraken_Species"], row["Kraken_taxid"]
    elif row["Kraken_Species"] == row["NT_Species"] and row["Kraken_Species"] != "-":
        return row["Kraken_Species"], row["Kraken_taxid"]
    elif row["NR_Species"] == row["NT_Species"] and row["NR_Species"] != "-":
        return row["NR_Species"], row["NR_taxid"]

    # If Kraken contains "Homo Sapiens" and NT contains "Naegleria fowleri" or other primates then its "Homo Sapiens"
    elif row["Kraken_Species"] == "Homo sapiens" and \
        (row["NT_Species"] == "Naegleria fowleri" or
        row["NR_Species"] =="Naegleria fowleri"):
        return row["Kraken_Species"], row["Kraken_taxid"]
    
    # Need to remove some more "Naegleria fowleri", so make it have to match two columns to call it "Naegleria fowleri"
    elif row["NT_Species"] == "Naegleria fowleri" and \
        row["Kraken_Species"] == "-" and \
        row["NR_Species"] == "-":
        return "-", "-"
    
    elif row["NR_Species"] == "Naegleria fowleri" and \
        row["Kraken_Species"] == "-" and \
        row["NT_Species"] == "-":
        return "-", "-"
    
    # If 2 columns contain "-", then print the column with stuff in it
    elif (
        (row["Kraken_Species"] == "-" or kraken_species == "-")
        and (row["NR_Species"] == row["NT_Species"] or nr_species == nt_species)
        and row["NR_Species"] != "-"
    ):
        return row["NR_Species"], row["NR_taxid"]  # or row["NT_Species"], row["NT_taxid"]
    elif (
        (row["NR_Species"] == "-" or nr_species == "-")
        and (row["Kraken_Species"] == row["NT_Species"] or kraken_species == nt_species)
        and row["Kraken_Species"] != "-"
    ):
        return row["Kraken_Species"], row["Kraken_taxid"]  # or row["NT_Species"], row["NT_taxid"]
    elif (
        (row["NT_Species"] == "-" or nt_species == "-")
        and (row["Kraken_Species"] == row["NR_Species"] or kraken_species == nr_species)
        and row["Kraken_Species"] != "-"
    ):
        return row["Kraken_Species"], row["Kraken_taxid"]  # or row["NR_Species"], row["NR_taxid"]
    

        
    # If Kraken is "-", and either NR or NT contains stuff within the "synthetic list", then return NT or NR thats not synthetic, if both synthetic then final is "-"
    elif row["Kraken_Species"] == '-':
        best_match = None
        for species in [row["NT_Species"], row["NR_Species"]]:
            for substring in synthetic_list:
                if substring.lower() not in species.lower():
                    best_match = species
                    best_match_taxid = row["MM2_taxid"] if species == row["NT_Species"] else row["NR_taxid"]
        if best_match is not None:
            return best_match, best_match_taxid
        else:
            return "-", "-"
        
    else:
        # Extract and compare genus if all previous conditions fail
        if kraken_genus and not nr_genus and not nt_genus:
            return kraken_genus, None
        elif not kraken_genus and nr_genus and not nt_genus:
            return nr_genus, None
        elif not kraken_genus and not nr_genus and nt_genus:
            return nt_genus, None
        elif kraken_genus == nr_genus:
            return kraken_genus, None
        elif kraken_genus == nt_genus:
            return kraken_genus, None
        elif nr_genus == nt_genus:
            return nr_genus, None
        elif (nt_genus == "Homo" or nr_genus == "Homo"):
            return "Homo sapiens", 9606
        else:
            return None, None

def taxid_assignment(row, taxdump):
    """
    Determining the consensus lineage if all three do not match. Basically searches the taxonomic tree 
    for the most common ancestor.
    """
    kraken_taxid = str(row["Kraken_taxid"]).replace('.0', '') if row["Kraken_taxid"] != "-" else "-"
    nr_taxid = str(row["NR_taxid"]).replace('.0', '') if row["NR_taxid"] != "-" else "-"
    nt_taxid = nt_taxid = str(row["MM2_taxid"]).replace('.0', '') if row["MM2_taxid"] != "-" else "-"

    taxid_list = [taxid for taxid in [kraken_taxid, nr_taxid, nt_taxid] if taxid != "-"]
    consensus_lineage = None
    consensus_taxid = None
    human_lineage = ['Homo', 'Hominidae', 'Primates', 'Mammalia', 'Chordata', 'Metazoa', 'Eukaryota']
    try:
        consensus_lineage = taxdump.lca(taxid_list).name
        consensus_taxid = taxdump.getTaxid(consensus_lineage)
    except KeyError as exc:
        # Log the error and continue without the offending taxid
        logging.error("Generated a KeyError: %s. Removing the offending taxid and retrying.", exc)
        off_taxid = exc.args[0]
        taxid_list = [taxid for taxid in taxid_list if taxid != off_taxid]

        # Retry with the updated taxid_list
        try:
            if len(taxid_list) > 1:
                consensus_lineage = taxdump.lca(taxid_list).name
                consensus_taxid = taxdump.getTaxid(consensus_lineage)
                logging.info("Yay!: %s removed successfully and ignored.", exc)
            else:
                consensus_lineage = "-"
                consensus_taxid = "-"
        except KeyError as exc:
            logging.error("Generated another KeyError: %s. Unable to find consensus taxid.", exc)
    
    # swap in Kraken species if consensus lineage is within the human lineage
    if row['Kraken_Species'] == 'Homo sapiens' and consensus_lineage in human_lineage:
        consensus_lineage = row['Kraken_Species']
        consensus_taxid = row['Kraken_taxid']

    return consensus_lineage, consensus_taxid

def fix_taxid(row, taxdump):
    if row['final_taxid'] == "-" and row['final'] != '-':
        og_taxid = row['final_taxid']
        species = None
        try:
            species = taxdump.getTaxid(row['final'])
            row['final_taxid'] = species
            logging.info(f"{row['final']} was assigned {og_taxid}, replacing with {species} via Taxdump")
        except KeyError as exc:
            logging.error("Generated an exception: %s", exc)
        if species is None:
            species = es.search_species(row['final'])
            row['final_taxid'] = species
            logging.info(f"{row['final']} was assigned {og_taxid}, replacing with {species} via Taxonomy Entrez")
    return row['final'], row['final_taxid']

