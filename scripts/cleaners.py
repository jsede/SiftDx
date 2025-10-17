import re
import pandas as pd
from Bio import SeqIO
import json
import assists
import logging
import species2lineage as sl
import acc2taxid as at
import entrez_search as es

exc_kingdom = [
    "Bacteria",
    "Eukaryota",
    "Viruses",
    'Archaea',
]
def read_csv(file, name, col):
    clean_df = pd.read_csv(
        file, sep="\t", header=None, names=[name, "accession"], usecols=[0, col]
    )
    return clean_df

def blast_cleanup(file, name, taxdump, dirpath, entrez_cred):
    loaded_cache = assists.load_cache_files(dirpath)
    blast_df = pd.read_csv(
            file, sep="\t", 
            header=None, 
            names=[
                name, 
                "accession",
                "pident", 
                "length", 
                "evalue",
                "bitscore",
                "qlen",
                "staxids"
                ], 
                usecols=[0,1,2,3,10,11,12,13],
                dtype = {
                    'accession': str,
                    'length': int, 
                    'pident': float, 
                    'qlen': int,
                    'evalue': float,
                    "bitscore": float,
                    'staxids': str}
        )
    taxid_df = pd.read_csv(
            file, sep="\t", header=None, names=["taxid"], usecols=[13], dtype=str
        )
    
    # Split semicolon-separated values and flatten into a single Series
    taxid_series = (
        taxid_df['taxid']
        .dropna()
        .str.split(';')
        .explode()
        .str.strip() 
    )
    taxids = taxid_series.drop_duplicates().tolist()
    get_lineage(taxids, taxdump, dirpath, entrez_cred)
    lineage_cache = assists.check_lineage_json(loaded_cache[10])

    dfs = [pd.DataFrame(entry, index=[0]) for entry in lineage_cache]
    lineage_df = pd.concat(dfs, ignore_index=True, sort=False)
    lineage_df['no rank'] = lineage_df['no rank'].replace('no rank', 'root')

    # Create 'superkingdom' column based on 'domain' and 'acellular root', occurs in the newer versions of taxdump
    if "superkingdom" not in lineage_df.columns:
        lineage_df["superkingdom"] = None
    if "acellular root" in lineage_df.columns:
        lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["acellular root"])
    if "domain" in lineage_df.columns:
        lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["domain"])

    lineage_df = lineage_df[['taxid', 'superkingdom', 'species']].dropna(how='all')
    lineage_df = lineage_df.astype(str).apply(lambda x: x.str.strip()) # convert all to string and strip whitespace
    taxid_to_kingdom = dict(zip(lineage_df["taxid"], lineage_df["superkingdom"]))
    taxid_to_species = dict(zip(lineage_df["taxid"], lineage_df["species"]))

    blast_df["superkingdom"] = blast_df["staxids"].map(taxid_to_kingdom)
    blast_df["BLAST_Species"] = blast_df["staxids"].map(taxid_to_species)

    blast_df['alnlen'] = blast_df['length']/blast_df['qlen']
    
    bacteria_condition = (blast_df['superkingdom'] == 'Bacteria') & (blast_df['alnlen'] > 0.5) & (blast_df['pident'] > 95)
    eukaryota_condition = (blast_df['superkingdom'] == 'Eukaryota') | (blast_df['superkingdom'] == 'Archaea') & (blast_df['alnlen'] > 0.5) & (blast_df['pident'] > 98)
    viruses_condition = (blast_df['superkingdom'] == 'Viruses') & (blast_df['alnlen'] > 0.5) & (blast_df['pident'] > 75)
    nothing_condition = ~blast_df['superkingdom'].isin(exc_kingdom)

    blast_df = blast_df[(bacteria_condition | eukaryota_condition | viruses_condition | nothing_condition)]
    blast_df.drop_duplicates(subset="NT", keep="first", inplace=True)
    blast_df['accession'] = blast_df['accession'].str.extract(r'([A-Z_]+\d+\.\d+)')
    full_blast_df = blast_df[["NT", "accession", "BLAST_Species", "staxids", "pident", "alnlen", "evalue", "bitscore"]]
    blast_df = blast_df[["NT", "accession", "BLAST_Species"]]
    return blast_df, full_blast_df


def minimap_cleanup(file, name, taxdump, database, dirpath, entrez_cred):
    loaded_cache = assists.load_cache_files(dirpath)
    minimap_df = pd.read_csv(
            file, sep="\t", header=None, 
            names=[name, "accession", "pident", "alen", "evalue", "bitscore", "qlen"], 
            usecols=[0,1,2,3,10,11,12],
            dtype = {"accession":str, "pident": float, "alen": int, "evalue": float, "bitscore": float, "qlen": int}
        )
    accession_df = pd.read_csv(
            file, sep="\t", header=None, names=["accession"], usecols=[1], dtype=str
        )
    accession_df.drop_duplicates(subset="accession", keep="first", inplace=True)
    mm2_accession_list = accession_df["accession"].dropna().to_list()
    at.nucl_accession_search(database, mm2_accession_list, dirpath, entrez_cred)

    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict = loaded_cache[3]
    minimap_df["MM2_taxid"] = minimap_df["accession"].map(taxid_dict)
    taxids = list(taxid_dict.values())
    taxids = list(set(taxids))
    
    get_lineage(taxids, taxdump, dirpath, entrez_cred)
    lineage_cache = assists.check_lineage_json(loaded_cache[10])
    dfs = [pd.DataFrame(entry, index=[0]) for entry in lineage_cache]
    lineage_df = pd.concat(dfs, ignore_index=True, sort=False)
    lineage_df['no rank'] = lineage_df['no rank'].replace('no rank', 'root')

    if "superkingdom" not in lineage_df.columns:
        lineage_df["superkingdom"] = None
    if "acellular root" in lineage_df.columns:
        lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["acellular root"])
    if "domain" in lineage_df.columns:
        lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["domain"])

    lineage_df = lineage_df[['taxid', 'superkingdom']].dropna(how='all')
    mapping_dict = dict(zip(lineage_df['taxid'], lineage_df['superkingdom']))
    minimap_df['superkingdom'] = minimap_df['MM2_taxid'].map(mapping_dict)
    minimap_df["alnlen"] = minimap_df['alen']/minimap_df['qlen']
    bacteria_condition = (minimap_df['superkingdom'] == 'Bacteria') & (minimap_df['alnlen'] > 0.5) & (minimap_df['pident'] > 95)
    eukaryota_condition = (minimap_df['superkingdom'] == 'Eukaryota') | (minimap_df['superkingdom'] == 'Archaea') & (minimap_df['alnlen'] > 0.5) & (minimap_df['pident'] > 98)
    viruses_condition = (minimap_df['superkingdom'] == 'Viruses') & (minimap_df['alnlen'] > 0.5) & (minimap_df['pident'] > 75)
    nothing_condition = ~minimap_df['superkingdom'].isin(exc_kingdom)

    minimap_df = minimap_df[(bacteria_condition | eukaryota_condition | viruses_condition | nothing_condition)]
    minimap_df.drop_duplicates(subset="MM2_NT", keep="first", inplace=True)
    full_minimap_df = minimap_df[["MM2_NT", "accession", "MM2_taxid", "pident", "alnlen", "evalue", "bitscore"]]
    minimap_df = minimap_df[["MM2_NT", "accession"]]
    return minimap_df, full_minimap_df
    
def diamond_cleanup(file, name, taxdump, database, dirpath, entrez_cred):
    loaded_cache = assists.load_cache_files(dirpath)
    diamond_df = pd.read_csv(
            file, sep="\t", header=None, 
            names=[name, "accession", "pident", "length", "evalue", "bitscore", "qlen", "staxids"], 
            usecols=[0,1,2,3,10,11,12,13], 
            dtype = {"pident":float, "length":int, "evalue": float, "bitscore": float, "qlen":int}
        )
    # accession_df = pd.read_csv(
    #         file, sep="\t", header=None, names=["accession"], usecols=[1], dtype=str
    #     )
    
    # accession_df.drop_duplicates(subset="accession", keep="first", inplace=True)
    # diamond_list = accession_df["accession"].dropna().to_list()
    # at.prot_accession_search(database, diamond_list, dirpath, entrez_cred)

    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict = loaded_cache[3]
    diamond_df["NR_taxid"] = diamond_df["accession"].map(taxid_dict)
    taxids = list(taxid_dict.values())
    taxids = list(set(taxids))

    get_lineage(taxids, taxdump, dirpath, entrez_cred)
    lineage_cache = assists.check_lineage_json(loaded_cache[10])
    dfs = [pd.DataFrame(entry, index=[0]) for entry in lineage_cache]
    lineage_df = pd.concat(dfs, ignore_index=True, sort=False)
    lineage_df["no rank"] = lineage_df["no rank"].replace("no rank", "root")
    
    if "superkingdom" not in lineage_df.columns:
        lineage_df["superkingdom"] = None
    if "acellular root" in lineage_df.columns:
        lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["acellular root"])
    if "domain" in lineage_df.columns:
        lineage_df["superkingdom"] = lineage_df["superkingdom"].fillna(lineage_df["domain"])
        
    lineage_df = lineage_df[["taxid", "superkingdom", "species"]].dropna(how="all").astype(str)

    taxid_to_kingdom = dict(zip(lineage_df["taxid"], lineage_df["superkingdom"]))
    taxid_to_species = dict(zip(lineage_df["taxid"], lineage_df["species"]))
    
    diamond_df["superkingdom"] = diamond_df["staxids"].map(taxid_to_kingdom)
    diamond_df["sscinames"] = diamond_df["staxids"].map(taxid_to_species)
    
    diamond_df['alnlen'] = (diamond_df['length']*3)/diamond_df['qlen']

    bacteria_condition = (diamond_df['superkingdom'] == 'Bacteria') & (diamond_df['alnlen'] > 0.5) & (diamond_df['pident'] > 95)
    eukaryota_condition = (diamond_df['superkingdom'] == 'Eukaryota') | (diamond_df['superkingdom'] == 'Archaea') & (diamond_df['alnlen'] > 0.5) & (diamond_df['pident'] > 98)
    viruses_condition = (diamond_df['superkingdom'] == 'Viruses') & (diamond_df['alnlen'] > 0.5) & (diamond_df['pident'] > 75)
    nothing_condition = ~diamond_df['superkingdom'].isin(exc_kingdom)

    diamond_df = diamond_df[(bacteria_condition | eukaryota_condition | viruses_condition | nothing_condition)]
    diamond_df.drop_duplicates(subset="NR", keep="first", inplace=True)
    full_diamond_df = diamond_df[["NR", "accession", "sscinames", "staxids", "pident", "alnlen", "evalue", "bitscore"]]
    diamond_df = diamond_df[["NR", "accession"]]
    return diamond_df, full_diamond_df


def get_lineage(taxon_list, taxdump, dirpath, entrez_cred):
    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict, prev_failed = loaded_cache[11], loaded_cache[13]
    taxid_lookup = {entry['taxid']: entry for entry in taxid_dict}
    taxid_lookup = {str(k): v for k, v in taxid_lookup.items()}
    
    for taxon in taxon_list:
        if taxon in taxid_lookup:
            logging.info(f"Taxid {taxon} with lineage already found in JSON file.")
            continue
        elif taxon in prev_failed:
            if taxon in taxid_lookup:
                logging.info(f"Taxid {taxon} was previously failed but is now found in JSON file.")
                continue
            else:
                logging.info(f"Taxid {taxon} has failed all attempts previously.")
                continue
        
        result = sl.get_taxon_lineage(int(taxon), taxdump, entrez_cred)
        if result is not None:
            with open(loaded_cache[10], "a") as file:
                json.dump(result, file)
                file.write("\n")
            taxid_lookup[taxon] = result
        else:
            logging.info(f"{taxon} not found, writing to failed_lineage.json")
            with open(loaded_cache[12], "a") as file:
                json.dump({taxon: "Not found"}, file)
                file.write("\n")

def convert_fasta_to_list(fasta_file):
    headers = []
    seq_lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.id
        length = len(record)
        headers.append(header)
        seq_lengths.append(length)
    return headers, seq_lengths

def extract_name_numreads(covstats_file):
    covstats_df = pd.read_csv(covstats_file, sep="\t", header=0, usecols=[0, 3], dtype = {"#rname": str ,"numreads": int})
    covstats_df['#rname'] = covstats_df['#rname'].str.extract(r'(\S+)')
    return covstats_df

def remove_brackets(text):
    return re.sub(r"\([^)]*\)", "", text)


def get_taxid(text):
    taxid = "-"
    taxid_pattern = r"\(taxid\s+(\d+)\)"
    matches = re.search(taxid_pattern, text)
    if matches:
        taxid = matches.group(1)
    return taxid


def remove_num_str(text):
    if isinstance(text, list):
        if isinstance(text[0], int):
            return '-'
        text = text[0] if len(text) > 0 else ''
    pattern = r"str\..*"
    result = re.sub(pattern, "", text)
    return result.strip()

def extract_genus(row, colname, name, taxdump, entrez_cred):
    result_dict = {}
    result = row[name]  # default fallback

    try:
        rank = sl.get_taxon_rank(row, colname, taxdump)

        if rank != 'no rank':
            if rank in ['species', 'subspecies']:
                parent = taxdump.getAncestry(row[colname])
                lineage = [node for node in parent if node.rank == 'genus']
                if lineage:
                    result = lineage[0].name
                else:
                    logging.warning(f"No genus found in lineage for {row[colname]}")
                    result_dict = es.search_lineage(row[colname], entrez_cred)
                    result = result_dict.get('genus', row[name])

            
    except Exception as e:
        logging.error(f"Error extracting genus for taxon {row[colname]}: {e}")

    return result

def extract_species(species):
    return " ".join(species.split()[:2]) if species != "-" else None

def rm_dupes_json(file):
    unique_species = {}
    with open(file, "r") as infile:
        for line in infile:
            json_data = json.loads(line)
            species_name = list(json_data.keys())[0]
            taxid = json_data[species_name]
            if species_name not in unique_species:
                unique_species[species_name] = taxid
            elif taxid != '-' and unique_species[species_name] == '-':
                unique_species[species_name] = taxid
        
    with open(file, "w") as outfile:
        for species_name, taxid in unique_species.items():
            json_data = {species_name: taxid}
            json.dump(json_data, outfile)
            outfile.write("\n")

def collapse_same_species(df, taxdump):
    # Columns to average
    numeric_cols = ['pident', 'alnlen', 'evalue', 'bitscore', 'Seq_Length', 'numreads']

    # Coerce problematic entries to NaN for numeric columns
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    agg_dict = {
        'numreads': 'sum',
        'Seq_Length': 'mean',
        'pident': 'mean',
        'alnlen': 'mean',
        'evalue': 'mean',
        'bitscore': 'mean',
        'final_taxid': lambda x: ','.join(sorted(set(map(str, x)))),
        'accession': lambda x: ','.join(sorted(set(map(str, x))))
    }

        # First, find rows with duplicate 'final' entries
    duplicate_final = df[df.duplicated(subset='final', keep=False)]
    if duplicate_final.empty:
        return df

    # Columns to include in aggregation
    aggregation_columns = list(agg_dict.keys())
    aggregation_columns.append('final')

    # Keep only the relevant columns for aggregation
    aggregation_df = duplicate_final[aggregation_columns]

    # Perform aggregation
    aggregated = aggregation_df.groupby('final').agg(agg_dict).reset_index()

    # Add back any taxonomy or metadata columns using first occurrence
    non_numeric_cols = [
        'superkingdom', 'kingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'subspecies'
    ]
    for col in non_numeric_cols:
        if col in df.columns:
            aggregated[col] = duplicate_final.groupby('final')[col].first().values

    # Drop original duplicated rows
    df = df[~df.duplicated(subset='final', keep=False)]

    # Concatenate cleaned df and aggregated one
    collapsed_df = pd.concat([df, aggregated], ignore_index=True)
    columns_to_drop = [
        'Fasta_Headers', 'Kraken_Species', 'NT_Species', 'NR_Species'
    ]

    collapsed_df.drop(columns=[col for col in columns_to_drop if col in collapsed_df.columns], inplace=True)
    collapsed_df = collapsed_df.fillna('-')
    collapsed_df.rename(columns={'final': 'taxon', 'Seq_length': 'avgseqlen'}, inplace=True)
    collapsed_df['taxon'] = collapsed_df['taxon'].replace('-', 'Unclassified')
    collapsed_df['rank'] = collapsed_df.apply(
        lambda row: sl.get_taxon_rank(row, 'final_taxid', taxdump),
        axis=1,
        result_type="expand"
    )

    desired_order = [
        'taxon', 'avgseqlen', 'numreads', 'final_taxid', 'accession',
        'pident', 'alnlen', 'evalue', 'bitscore',
        'rank', 'superkingdom', 'kingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species', 'subspecies'
    ]

    # Keep only the columns that exist in the dataframe
    existing_columns = [col for col in desired_order if col in collapsed_df.columns]
    collapsed_df = collapsed_df[existing_columns]
    return collapsed_df