import logging
import entrez_search as es
import pandas as pd
import taxidTools

def get_taxon_rank(row, taxid_name, taxdump):
    rank = None
    taxonomy_levels = [
        'subspecies', 'species', 'genus', 'family', 'order',
        'class', 'phylum', 'kingdom', 'superkingdom'
    ]
    try:
        taxid = row[taxid_name]
        if taxid not in ['-', '1', '0']:
            rank = taxdump.getRank(int(float(taxid)))
        else:
            rank = "Not assigned"
    except (KeyError, ValueError, TypeError) as exc:
        logging.error(f"Could not find rank for {taxid_name}: %s", exc)
    
    # Infer rank if not assigned or is 'no rank'
    if rank is None:
        for level in taxonomy_levels:
            val = row.get(level, '')
            if pd.notna(val) and str(val).strip() != '-':
                rank = level
                break
    
    return rank

def get_taxon_lineage(taxid, taxdump, entrez_cred):
    result_dict = {}
    try:
        lineage = taxdump.getAncestry(taxid)
        rank = [node.rank for node in lineage]
        taxonomy = [node.name for node in lineage]

        result_dict['taxid'] = taxid

        rank_len = len(rank)
        taxonomy_len = len(taxonomy)
        if rank_len == taxonomy_len:
            result_dict.update({rank[i]: taxonomy[i] for i in range(rank_len)})
        else:
            logging.error("The Rank and Taxonomy are different")
            
    except Exception as exc:
        logging.error("Could not find lineage via Taxdump: %s", exc)
        logging.info("Attempting to find lineage via Entrez")
        result_dict = es.search_lineage(taxid, entrez_cred)
    return result_dict