import taxidTools
import logging
import concurrent.futures
import json
import assists
import entrez_search as es

def taxid2lineage(dirpath, taxdump, entrez_cred):

    loaded_cache = assists.load_cache_files(dirpath)

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        species_search = []
        taxid_data = loaded_cache[3]
        species_dict = loaded_cache[4]

        # attempt running of new database
        for accession, taxid in taxid_data.items():
            # Check if the accession is already in the existing data
            if accession in species_dict:
                logging.info(f"Accession {accession} already found in JSON file.")
                continue  # Skip fetching species name if it's already present
            elif accession in loaded_cache[7]:
                logging.info(
                    f"Accession {accession} has failed all attempts previously"
                )
                continue  # Skip because it has failed in previous attempts.
            task = executor.submit(get_species, accession, taxid, taxdump, entrez_cred)
            species_search.append(task)

        for task in concurrent.futures.as_completed(species_search):
            result = task.result()
            try:
                if result:
                    accession, taxid, species = result
                    logging.info(f"Found {accession}! Species is {species}")
                    species_dict[accession] = species
                    if species is not None:
                        with open(loaded_cache[1], "a") as file:
                            json.dump({accession: species}, file)
                            file.write("\n")
                    elif species is None:
                        with open(loaded_cache[6], "a") as file:
                            json.dump({accession: "Not found"}, file)
                            file.write("\n")
            except Exception as exc:
                logging.error("%s generated an exception: %s", result, exc)


def get_species(accession, taxid, taxdump, entrez_cred):
    try:
        species = None
        species = taxdump.getName(taxid)
        lineage = taxdump.getAncestry(taxid)
        lineage.filter()
        [node.name for node in lineage]
    except KeyError as exc:
        logging.error("The following taxid generated an exception: %s", exc)
    if species is None:
        species = es.search_taxid(taxid, entrez_cred)
    return accession, taxid, species

