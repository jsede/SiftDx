import logging
import entrez_search as es
import taxidTools

def get_taxon_rank(row, taxdump):
    rank = None
    taxonomy_levels = [
        'subspecies', 'species', 'genus', 'family', 'order',
        'class', 'phylum', 'kingdom', 'superkingdom'
    ]
    try:
        taxid = row['final_taxid']
        if taxid != '-':
            rank = taxdump.getRank(int(float(taxid)))
    except (KeyError, ValueError, TypeError) as exc:
        logging.error("Generated an exception: %s", exc)
    
    # Infer rank if not assigned or is 'no rank'
    if rank == 'no rank':
        for level in taxonomy_levels:
            val = row.get(level, '')
            if isinstance(val, str) and val.strip() != '':
                rank = level
                break
    if rank is None:
        rank = "Not assigned"
    
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
        logging.error("Generated an exception: %s", exc)
        result_dict = es.search_lineage(taxid, entrez_cred)
    return result_dict