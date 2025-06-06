import os
import json
import subprocess
import sys
import base64
import logging
import shutil
import os.path
import species2lineage as sl


def run_cmd(command):
    """
    Run commands with error outputs.
    """
    logging.info("Running command: %s", command)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if result.returncode != 0:
        logging.critical("Failed to run command: %s", result.args)
        logging.critical("stdout: %s", result.stdout)
        logging.critical("stderr: %s", result.stderr)
        sys.exit(1)
    return result


def check_files(file):
    """
    Check input files if they exist and have contents
    """

    if os.path.isfile(file) is True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = (
            file
            + " either file is does not exist or is empty, please check files. Please Note this."
        )
        logging.critical(msg)


def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) is True:
        truemsg = folder + " folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making folder"
        logging.info(msg)


def check_dependencies(cmd_exec):
    cmd_path = shutil.which(cmd_exec)

def base64encode(file_path):
    with open(file_path, "rb") as file:
        file_data = file.read()
        
    base64file = base64.b64encode(file_data).decode("utf-8")
    return base64file

def copy_neg2model(name, tmp_path, neg_path, model_path):
    preprocessing_files = [
        "/combined_reads_contigs_file.fa",
        "/ercc_coverage.txt",
        "/reads_mapped_to_contigs.cov_stats",
        "/fastp.json",
        "/fullyQc_fwd.fq",
    ]

    alignments_files = [
        '/accession_species.json',
        '/accession_taxid.json',
        '/failed_accessions.json',
        '/taxid_species.json',
        "/minimap2_lr_out.paf",
        "/minimap2_contig_out.paf",
        "/nr_alignments_file.tsv",
        "/nt_alignments_sr_blast.tsv",
        "/combined_rc.kraken.txt",
    ]

    shutil.copytree(model_path, tmp_path, dirs_exist_ok=True)
    try:
        os.mkdir(tmp_path + f"/{name}")
        os.mkdir(tmp_path + f"/{name}")
    except FileExistsError:
        pass
    for file in preprocessing_files:
        if os.path.exists(neg_path + file):        
            shutil.copy(neg_path + file, tmp_path + f"/{name}")
        else:
            logging.info(f"Missing {file}")
    for file in alignments_files:
        if os.path.exists(neg_path + file):
            shutil.copy(neg_path + file, tmp_path + f"/{name}")
        else:
            logging.info(f"Missing {file}")

def total_reads(file):
    with open(file) as json_file:
        fastp_json = json.load(json_file)
    total_reads = fastp_json['summary']['before_filtering']['total_reads']
    return total_reads

def check_previous_json(cache_file):
    existing_taxonomy_data = {}
    try:
        with open(cache_file, "r") as file:
            for line in file:
                data = json.loads(line)
                existing_taxonomy_data.update(data)
    except FileNotFoundError:
        existing_taxonomy_data = {}
    return existing_taxonomy_data

def check_lineage_json(cache_file):
    existing_lineage_data = {}
    try:
        with open(cache_file, "r") as file:
            json_data = file.read().split('\n')
        data = [json.loads(entry) for entry in json_data if entry.strip()]
        existing_lineage_data = data
    except FileNotFoundError:
        existing_lineage_data = {}
    return existing_lineage_data

def load_cache_files(dirpath):
    cache_file = dirpath + "/accession_taxid.json"
    taxid_file = dirpath + "/accession_species.json"
    failed_file = dirpath + "/failed_accessions.json"
    failed_taxid = dirpath + "/failed_taxids.json"
    taxid_species_file = dirpath + "/taxid_species.json"
    lineage_file = dirpath + "/taxid_lineage.json"
    failed_lineage = dirpath + "/failed_lineage.json"

    prev_json = check_previous_json(cache_file)
    prev_taxid_json = check_previous_json(taxid_file)
    prev_failed = check_previous_json(failed_file)
    prev_failed_taxid = check_previous_json(failed_taxid)
    prev_taxid_species_json = check_previous_json(taxid_species_file)
    prev_lineage_json = check_lineage_json(lineage_file)
    prev_failed_lineage = check_previous_json(failed_lineage)

    return (
        cache_file, #0
        taxid_file, #1
        failed_file, #2
        prev_json, #3
        prev_taxid_json, #4
        prev_failed, #5
        failed_taxid, #6
        prev_failed_taxid, #7
        taxid_species_file, #8
        prev_taxid_species_json, #9
        lineage_file, #10
        prev_lineage_json, #11
        failed_lineage, #12
        prev_failed_lineage #13
    )


