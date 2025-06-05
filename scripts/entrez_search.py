import time
import logging
import signal
from Bio import Entrez
from urllib.request import HTTPError

timeout_minutes = 5

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Timed out")

def search_nt_taxid(accession, entrez_cred):
    Entrez.email = entrez_cred[0]
    Entrez.api_key = entrez_cred[1]

    max_retries = 3
    taxid = None
    search_handle, search_record = None, None
    for retry in range(max_retries):
        try:
            search_handle = Entrez.esummary(db="nuccore", id=accession)
            search_record = Entrez.read(search_handle)
            logging.info(f"Searching for {accession} in Nuccore in Entrez")
        
        except Exception as e:
            logging.info(f"Search for {accession} in Nuccore FAILED, {e}")
            if retry < max_retries - 1:
                time.sleep(2)
            else:
                logging.error(f"Maximum retries reached for {accession}. Giving up.")
                raise  # Reraise the exception if max retries reached
        if search_record:
            taxid = search_record[0]["TaxId"]
            taxid = int(str(taxid).split("(")[1].split(",")[0])
        else:
            taxid = None
        if search_handle:
            search_handle.close()

        return accession, taxid

def search_prot_taxid(accession, entrez_cred):
    Entrez.email = entrez_cred[0]
    Entrez.api_key = entrez_cred[1]
    
    max_retries = 3
    timeout_seconds = timeout_minutes * 60
    search_handle, search_record = None, None
    for retry in range(max_retries):
        try:
            logging.info(f"Searching for {accession} in Protein Entrez")
            
            # timeout handling.
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(timeout_seconds)

            search_handle = Entrez.esummary(db="protein", id=accession)
            try:
                search_record = Entrez.read(search_handle)
            except TimeoutException:
                logging.info(f"Search for {accession} timed out.")
                return accession, None
            except Exception as e:
                logging.info(f"Could not find {accession}: {e}")
                return accession, None
            finally:
                signal.alarm(0)

            if search_record:
                taxid = search_record[0]["TaxId"]
                taxid = int(str(taxid).split("(")[1].split(",")[0])
            else:
                taxid = None
            if search_handle:
                search_handle.close()

            return accession, taxid

        except Exception as e:
            logging.info(f"Search for {accession} in Protein FAILED, {e}")
            if retry < max_retries - 1:
                time.sleep(2)
            else:
                logging.error(f"Maximum retries reached for {accession}. Giving up.")
                return accession, None


def search_taxid(taxid, entrez_cred):
    Entrez.email = entrez_cred[0]
    Entrez.api_key = entrez_cred[1]
    
    max_retries = 3
    species = None
    search_handle, search_record = None, None
    for retry in range(max_retries):
        try:
            logging.info(f"Searching for {taxid} in Taxonomy Entrez")
            try:
                search_handle = Entrez.efetch(db="taxonomy", id=taxid)
                search_record = Entrez.read(search_handle)
            except Exception as e:
                logging.info(f"Could not find {taxid}: {e}")
                return taxid, None
            if search_record:
                species = search_record[0]["ScientificName"]
            if search_handle:
                search_handle.close()

            return species

        except Exception as e:
            logging.info(f"Search for {taxid} in Taxonomy FAILED, {e}")
            if retry < max_retries - 1:
                time.sleep(2)
            else:
                logging.error(f"Maximum retries reached for {taxid}. Giving up.")
                raise  # Reraise the exception if max retries reached

def search_species(name, entrez_cred):
    Entrez.email = entrez_cred[0]
    Entrez.api_key = entrez_cred[1]
    
    max_retries = 3
    taxid = '-'
    search_handle, search_record = None, None
    for retry in range(max_retries):
        try:
            logging.info(f"Searching for {name} in Taxonomy Entrez")
            try:
                search_handle = Entrez.esearch(db="taxonomy", term=name)
                search_record = Entrez.read(search_handle)
            except Exception as e:
                logging.info(f"Could not find {name}: {e}")
                return name, None
            if search_record:
                if search_record["IdList"]:
                    taxid = search_record["IdList"][0]
                else:
                    taxid = '-'
            else:
                taxid = '-'
            if search_handle:
                search_handle.close()
            return taxid

        except Exception as e:
            logging.info(f"Search for {name} in Taxonomy FAILED, {e}")
            if retry < max_retries - 1:
                time.sleep(2)
            else:
                logging.error(f"Maximum retries reached for {name}. Giving up.")
                raise  # Reraise the exception if max retries reached

def search_lineage(taxid, entrez_cred):
    desired_cols = [
        'taxid', 'superkingdom', 'kingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species', 'subspecies'
    ]
    Entrez.email = entrez_cred[0]
    Entrez.api_key = entrez_cred[1]
    
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    records = Entrez.read(handle)
    if not records:
        return None

    record = records[0]
    lineage_ex = record.get("LineageEx", [])
    
    lineage_dict = {'taxid': int(taxid)}
    for item in lineage_ex:
        rank = item['Rank']
        name = item['ScientificName']
        if rank and rank != 'no rank':
            lineage_dict[rank] = name
        elif rank == 'no rank' and 'no rank' not in lineage_dict:
            lineage_dict['no rank'] = name  # e.g. "root"
    
    # Add the current species itself
    if record['Rank'] == 'species':
        lineage_dict['species'] = record['ScientificName']
    elif record['Rank'] and record['ScientificName']:
        lineage_dict[record['Rank']] = record['ScientificName']
    if "domain" in lineage_dict:
        lineage_dict["superkingdom"] = lineage_dict.pop("domain")
    normalised = {key: lineage_dict.get(key, '-') for key in desired_cols}

    return normalised
