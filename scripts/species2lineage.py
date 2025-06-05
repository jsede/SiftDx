import logging
import entrez_search as es
import taxidTools

def get_taxon_rank(row, taxdump):
    rank = None
    try:
        rank = taxdump.getRank(row['final_taxid'])
    except KeyError as exc:
        logging.error("Generated an exception: %s", exc)
    if rank is None:
        rank = "Not assigned"
    return rank

def get_taxon_lineage(taxid, taxdump, entrez_cred):
    result_dict = {}
    success = False
    try:
        lineage = taxdump.getAncestry(taxid)
        rank = [node.rank for node in lineage]
        taxonomy = [node.name for node in lineage]

        result_dict['taxid'] = taxid

        rank_len = len(rank)
        taxonomy_len = len(taxonomy)
        if rank_len == taxonomy_len:
            result_dict.update({rank[i]: taxonomy[i] for i in range(rank_len)})
            success = True
        else:
            logging.error("The Rank and Taxonomy are different")
            
    except Exception as exc:
        logging.error("Generated an exception: %s", exc)
        result_dict = es.search_lineage(taxid, entrez_cred)
    return result_dict