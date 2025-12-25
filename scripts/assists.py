import os
import json
import subprocess
import sys
import base64
import logging
import shutil
import os.path
import species2lineage as sl
import pandas as pd


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

def taxon_split(zscore_df, taxon, header):
    taxon_df = zscore_df[zscore_df[header] == taxon]
    return taxon_df

def get_parasites(zscore_df, parasite_list):
    parasite_df = zscore_df[zscore_df['species'].isin(parasite_list['#Organism/Name'])]
    return parasite_df

def flag_pathogens(row, pathogen_db, contam_db):
    taxon = row['taxon']
    match_species = pathogen_db[(pathogen_db['Species'].str.lower() == taxon.lower()) | (pathogen_db['AltNames'].str.lower() == taxon.lower())]
    contam_species = contam_db[(contam_db['Species'] == taxon.lower())]
    if not contam_species.empty:
            taxon = f"{taxon}"
            if contam_species['Reason'].isna().all() is True:
                contam_status = contam_status.iloc[0]['Statement']
            else:
                contam_status = ". " + taxon + " " + contam_species.iloc[0]['Statement'] + " " + contam_species.iloc[0]['Reason']
    else:
        contam_status = ""

    if not match_species.empty:
        taxon = f"{taxon} \u2757"
        status = match_species.iloc[0]['Status']
        additions = []
        if not pd.isnull(match_species.iloc[0]['Disease_type']):
            additions.append(" causing a " + str(match_species.iloc[0]['Disease_type']) + " infection")
        if not pd.isnull(match_species.iloc[0]['Disease_Name']):
            additions.append(" called " + str(match_species.iloc[0]['Disease_Name']))
        if not pd.isnull(match_species.iloc[0]['Commensal']):
            commensal_locations = match_species.iloc[0]['Commensal'].split('; ')
            if len(commensal_locations) > 1:
                commensal_locations[-2] += " and " + commensal_locations[-1]
                commensal_locations = commensal_locations[:-1]
            additions.append(". " + match_species.iloc[0]['Species'] + " is a commensal microbe found in the " + ', '.join(commensal_locations))
        if not pd.isnull(match_species.iloc[0]['NNDSS_Notifiable']):
            additions.append(", it is also a notifiable disease in NSW")
        if additions:
            status += "".join(additions)
    elif taxon == "-":
        status = ""
    else:
        status = " is not a known pathogen"

    return taxon, status

def html_loop(taxon_df, negative):
    html_output = []
 
    # sort by zscore or rpm
    sort_col = "zscore" if negative != "None" else "rpm_sample"
    taxon_df = taxon_df.sort_values(by=sort_col, ascending=False)
    
    # identify genus-level rows and species-level rows
    genus_rows = taxon_df[(taxon_df['rank'] == 'genus') & (taxon_df['genus'] != '-')]
    species_rows = taxon_df[taxon_df['species'].notna() & (taxon_df['species'] != '-')]

    # identify species that have genus-level representatives
    genus_set = set(genus_rows['genus'])
    species_with_genus = species_rows[species_rows['genus'].str.strip().isin(genus_set)]
    nested_labels = set(species_with_genus['taxon_label'])
    
    # drop the species with genus representatives from the main taxon_df
    taxon_df = taxon_df[~taxon_df['taxon_label'].isin(nested_labels)]

    # Process the rows without genus-level representatives first
    for _, row in taxon_df.iterrows():
        # combine the lineage info into one line.
        lineage_cols = [
            'superkingdom', 'kingdom', 'phylum', 'class',
            'order', 'family', 'genus', 'species', 'subspecies'
        ]
        lineage_parts = []
        for col in lineage_cols:
            val = row.get(col, '')
            if pd.isna(val) or val == '-':
                continue  # skip this level but keep going
            lineage_parts.append(val)
        lineage = ' > '.join(lineage_parts)

        # modify the percentage identity
        try:
            pident_val = float(row['pident'])
            pident_str = f"{pident_val:.2f}%"
        except (ValueError, TypeError):
            pident_str = '-'

        if row['genus'] in genus_set:
            # get all species under this genus from species_with_genus
            species_rows = species_with_genus[species_with_genus['genus'] == row['genus']]
            # build nested species accordions
            accordion_html = ""

            # start genus accordion
            sp_metric = f"Z-score: {row['zscore']}" if negative != "None" else f"RPM: {row['rpm_sample']}"
            genus_html = f"""
            <button class="accordion">
                <span class="button-content">
                    <span class="taxon">{row['genus']}</span>
                    <span class="zscore">{sp_metric}</span>
                </span>
            </button>
            <div class="panel">
                <p><b>Lineage:</b> {lineage}</p>
            """
            
            for _, sp_row in species_rows.iterrows():
                try:
                    sp_pident_str = f"{float(sp_row['pident']):.2f}%"
                    lineage_parts = []
                    for col in lineage_cols:
                        val = sp_row.get(col, '')
                        if pd.isna(val) or val == '-':
                            continue  # skip this level but keep going
                        lineage_parts.append(val)
                    sp_lineage = ' > '.join(lineage_parts)
                except:
                    sp_pident_str = "-"
                
                sp_metric = f"Z-score: {sp_row['zscore']}" if negative != "None" else f"RPM: {sp_row['rpm_sample']}"

                genus_html += f"""
                <div style="margin-left:20px;">
                    <button class="accordion" style="background-color:#f1f1f1;">
                        <span class="button-content">
                            <span class="taxon">{sp_row['taxon_label']}</span>
                            <span class="zscore">{sp_metric}</span>
                        </span>
                    </button>
                    <div class="panel">
                        <p>
                            {sp_row['taxon']} {sp_row['taxon_status']}<br>
                            Closest NCBI Taxonomy ID: {sp_row['final_taxid']}<br>
                            Closest NCBI Accession: {sp_row['accession']}<br>
                            Percentage Identity: {sp_pident_str}<br>
                            Average Alignment Length: {sp_row['alnlen']}<br>
                            E-Value: {sp_row['evalue']}<br>
                            Bitscore: {sp_row['bitscore']}<br>
                """
                if negative != "None":
                    genus_html += f"RPM in Sample: {sp_row['rpm_sample']} &nbsp;&nbsp;&nbsp; RPM in NegCtrl: {sp_row['rpm_ctrl']}<br>"
                else:
                    genus_html += f"RPM: {sp_row['rpm_sample']}<br>"

                genus_html += f"""
                            Lineage: {sp_lineage}<br>
                        </p>
                    </div>
                </div>
                """
            genus_html += "</div>"  # close genus panel
            accordion_html = genus_html

        else:
            if negative != "None":
                accordion_html = f"""
                <button class="accordion">
                    <span class="button-content">
                        <span class="taxon">{row['taxon_label']}</span>
                        <span class="zscore">Z-score: {row['zscore']}</span>
                    </span>
                </button>
                <div class="panel">
                    <p>
                        {row['taxon']} {row['taxon_status']}<br>
                        Closest NCBI Taxonomy ID: {row['final_taxid']}<br>
                        Closest NCBI Accession: {row['accession']}<br>
                        Percentage Identity: {pident_str}<br>
                        Average Alignment Length: {row['alnlen']}<br>
                        E-Value: {row['evalue']}<br>
                        Bitscore: {row['bitscore']}<br>
                        RPM in Sample: {row['rpm_sample']}&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;RPM in NegCtrl: {row['rpm_ctrl']}<br>
                        Lineage: {lineage}<br>
                    </p>
                </div>
                """
            else:
                accordion_html = f"""
                <button class="accordion">
                    <span class="button-content">
                        <span class="taxon">{row['taxon_label']}</span>
                        <span class="zscore">RPM: {row['rpm_sample']}</span>
                    </span>
                </button>
                <div class="panel">
                    <p>
                        {row['taxon']} {row['taxon_status']}<br>
                        Closest NCBI Taxonomy ID: {row['final_taxid']}<br>
                        Closest NCBI Accession: {row['accession']}<br>
                        Percentage Identity: {pident_str}<br>
                        Average Alignment Length: {row['alnlen']}<br>
                        E-Value: {row['evalue']}<br>
                        Bitscore: {row['bitscore']}<br>
                        Lineage: {lineage}<br>
                    </p>
                </div>
                """
        html_output.append(accordion_html)
    return html_output

def get_first_taxid_rank(row, taxdump):
    taxid_field = row.get("final_taxid", "")
    if pd.isna(taxid_field) or str(taxid_field).strip() in ["-", "", "None"]:
        return "-"
    
    # Split on ";" and take the first non-empty taxid
    first_taxid = str(taxid_field).split(";")[0].strip()
    if not first_taxid:
        return "-"
    
    # Replace the row's 'final_taxid' temporarily
    temp_row = row.copy()
    temp_row["final_taxid"] = first_taxid
    
    return sl.get_taxon_rank(temp_row, "final_taxid", taxdump)