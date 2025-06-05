import sqlite3
import logging
import json
import assists
import entrez_search as es

def nucl_accession_search(database, accession_list, dirpath, entrez_cred):
    db_list = [
        "nucl_gb",
        "nucl_wgs"
    ]
    nucl_gb_path = database + "/" + db_list[0] + ".sqlite"
    conn = sqlite3.connect(nucl_gb_path)
    cursor = conn.cursor()
    unfound_taxid = []

    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict, prev_failed = loaded_cache[3], loaded_cache[5]

    for query in accession_list:
        if query in taxid_dict:
            logging.info(f"Accession {query} already found in JSON file.")
            continue  # Skip fetching species name if it's already present
        elif query in prev_failed:
                logging.info(f"Accession {query} has failed all attempts previously")
                continue  # Skip because it has failed in previous attempts.
        else:
            accession, taxid = sql_search(query, cursor, db_list[0])
            if taxid is not None:
                logging.info(f"Found {accession} in {db_list[0]}! taxid is {taxid}")
                taxid_dict[accession] = taxid
                with open(loaded_cache[0], "a") as file:
                    json.dump({accession: int(taxid)}, file)
                    file.write("\n")
            else:
                unfound_taxid.append(
                    accession
                )
    
    if len(unfound_taxid) > 0:
        nucl_wgs_path = database + "/" + db_list[1] + ".sqlite"
        conn = sqlite3.connect(nucl_wgs_path)
        cursor = conn.cursor()
        for query in unfound_taxid:
            if query in taxid_dict:
                logging.info(f"Accession {query} already found in JSON file.")
                continue  # Skip fetching species name if it's already present
            elif query in prev_failed:
                    logging.info(f"Accession {query} has failed all attempts previously")
                    continue  # Skip because it has failed in previous attempts.
            else:
                accession, taxid = sql_search(query, cursor, db_list[1])
                if taxid is not None:
                    logging.info(f"Found {accession} in {db_list[1]}! taxid is {taxid}")
                    taxid_dict[accession] = taxid
                    with open(loaded_cache[0], "a") as file:
                        json.dump({accession: int(taxid)}, file)
                        file.write("\n")
                else:
                    logging.info(f"Did not find {accession} in {db_list[1]}! Searching Entrez")
                    accession, taxid = es.search_nt_taxid(query, entrez_cred)
                    if taxid is not None:
                        logging.info(f"Found {accession} in Entrez! taxid is {taxid}")
                        with open(loaded_cache[0], "a") as file:
                            json.dump({accession: int(taxid)}, file)
                            file.write("\n")
                    else:
                        logging.info(f"{accession} not found, writing to failed_accessions.json")
                        with open(loaded_cache[2], "a") as file:
                            json.dump({accession: "Not found"}, file)
                            file.write("\n")

def prot_accession_search(database, prot_acc_list, dirpath, entrez_cred):
    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict, prev_failed = loaded_cache[3], loaded_cache[5]
    dead_prot = database + "/dead_prot.sqlite"
    db_list = [
        "prot",
        "dead_prot"
    ]
    prot_path = database + "/" + db_list[0] + ".sqlite"
    conn = sqlite3.connect(prot_path)
    cursor = conn.cursor()
    unfound_taxid = []

    loaded_cache = assists.load_cache_files(dirpath)
    taxid_dict, prev_failed = loaded_cache[3], loaded_cache[5]

    for query in prot_acc_list:
        if query in taxid_dict:
            logging.info(f"Accession {query} already found in JSON file.")
            continue  # Skip fetching species name if it's already present
        elif query in prev_failed:
                logging.info(f"Accession {query} has failed all attempts previously")
                continue  # Skip because it has failed in previous attempts.
        else:
            accession, taxid = sql_search(query, cursor, db_list[0])
            if taxid is not None:
                logging.info(f"Found {accession} in {db_list[0]}! taxid is {taxid}")
                taxid_dict[accession] = taxid
                with open(loaded_cache[0], "a") as file:
                    json.dump({accession: int(taxid)}, file)
                    file.write("\n")
            else:
                unfound_taxid.append(
                    accession
                )
    
    if len(unfound_taxid) > 0:
        dead_prot = database + "/" + db_list[1] + ".sqlite"
        conn = sqlite3.connect(dead_prot)
        cursor = conn.cursor()
        for query in unfound_taxid:
            if query in taxid_dict:
                logging.info(f"Accession {query} already found in JSON file.")
                continue  # Skip fetching species name if it's already present
            elif query in prev_failed:
                    logging.info(f"Accession {query} has failed all attempts previously")
                    continue  # Skip because it has failed in previous attempts.
            else:
                accession, taxid = sql_search(query, cursor, db_list[1])
                if taxid is not None:
                    logging.info(f"Found {accession} in {db_list[1]}! taxid is {taxid}")
                    taxid_dict[accession] = taxid
                    with open(loaded_cache[0], "a") as file:
                        json.dump({accession: int(taxid)}, file)
                        file.write("\n")
                else:
                    logging.info(f"Did not find {accession} in {db_list[1]}! Searching Entrez")
                    accession, taxid = es.search_prot_taxid(query, entrez_cred)
                    if taxid is not None:
                        logging.info(f"Found {accession} in Entrez! taxid is {taxid}")
                        with open(loaded_cache[0], "a") as file:
                            json.dump({accession: int(taxid)}, file)
                            file.write("\n")
                    else:
                        logging.info(f"{accession} not found, writing to failed_accessions.json")
                        with open(loaded_cache[2], "a") as file:
                            json.dump({accession: "Not found"}, file)
                            file.write("\n")

def sql_search(accession, cursor, db_type):
    try:
        sql_query = f'SELECT * FROM ({db_type}) WHERE "accession.version" = ?'
        cursor.execute(sql_query, (accession,))
        result = cursor.fetchall()
        if result:
            if db_type != "prot":
                taxid = result[0][2]
            elif db_type == "prot":
                taxid = result[0][1]
        else:
            taxid = None
        return accession, taxid
    except Exception as exc:
        logging.error("%s generated an exception: %s", exc)
