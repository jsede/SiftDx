import os
import sys

def gen_pipeline_summary(summary_file):
    dirpath = os.path.dirname(os.path.abspath(summary_file))
    workdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    template_html = workdir + "/assets/pipeline_template.html"
    summary_data = {}
    ercc_found = False
    sortmerna_found = False
    with open(summary_file, 'r') as summary_info:
        for line in summary_info:
            if ':' in line:
                key, value = line.strip().split(':', 1)
                summary_data[key.strip()] = value.strip()
            if 'ercc' in line.lower():
                ercc_found = True
            if 'sortmerna' in line.lower():
                sortmerna_found = True
    
    # special case for ERCC and Sortmerna
    ercc_html = "-"
    if ercc_found:
        ercc_html = f"""
            <b>Number of ERCC Reads Removed: {summary_data.get("ercc_reads", "N/A")},</b> <br>
            <b>Number of non-ERCC Host Reads: {summary_data.get("non-ercc_reads", "N/A")}</b> <br>
        """
    
    sortmerna_html = "DNA mode activated, this step was skipped. <br>"
    if sortmerna_found:
        sortmerna_html = f"""
            <b>Number of Reads After rRNA Depletion: {summary_data.get("rrna_reads_post_depletion", "N/A")}</b> <br>
            <b>Number of Non-host rRNA Reads: {summary_data.get("non_host_rrna_reads", "N/A")}</b> <br>
        """

    # Read HTML template
    with open(template_html, 'r', encoding='utf-8') as f:
        html_content = f.read()

    # Map summary keys to HTML placeholders
    replace_dict = {
        "py_prefastp_ph": summary_data.get("fastp_before_total_reads", "N/A"),
        "py_fastp_after_ph": summary_data.get("fastp_after_total_reads", "N/A"),
        "py_duplication_rate_ph": summary_data.get("duplication_rate", "N/A"),
        "py_low_quality_reads_ph": summary_data.get("low_quality_reads", "N/A"),
        "py_prinseq_good_ph": summary_data.get("prinseq_good", "N/A"),
        "py_lc_reads_ph": summary_data.get("lc_reads", "N/A"),
        "py_kraken2_human_depleted_ph": summary_data.get("kraken2_human_depleted", "N/A"),
        "py_bowtie2_human_depleted_ph": summary_data.get("bowtie2_human_depleted", "N/A"),
        "py_ercc_ph": ercc_html,
        "py_sortmerna_ph": sortmerna_html,
        "py_number_of_contigs_ph": summary_data.get("number_of_contigs", "N/A"),
        "py_shortest_contig_ph": summary_data.get("shortest_contig", "N/A"),
        "py_longest_contig_ph": summary_data.get("longest_contig", "N/A"),
        "py_avg_contig_length_ph": summary_data.get("avg_contig_length", "N/A"),
        "py_n50_ph": summary_data.get("n50", "N/A"),
        "py_total_unassembled_reads_ph": summary_data.get("total_unassembled_reads", "N/A"),
        "py_total_unassembled_shorter_reads_ph": summary_data.get("total_unassembled_shorter_reads", "N/A"),
        "py_total_unassembled_longer_reads_ph": summary_data.get("total_unassembled_longer_reads", "N/A"),
        "py_total_assembled_reads_ph": summary_data.get("total_assembled_reads", "N/A"),
        "py_number_of_identified_taxa_ph": summary_data.get("number_of_identified_taxa", "N/A"),
        "py_number_of_contigs_unassigned_ph": summary_data.get("number_of_contigs_unassigned", "N/A"),
        "py_number_of_contigs_assigned_ph": summary_data.get("number_of_contigs_assigned", "N/A"),
        "py_number_of_reads_unassigned_ph": summary_data.get("number_of_reads_unassigned", "N/A"),
        "py_number_of_reads_assigned_ph": summary_data.get("number_of_reads_assigned", "N/A"),
        "py_number_of_kraken_assigned_taxa_ph": summary_data.get("number_of_kraken_assigned_taxa", "N/A"),
        "py_number_of_nt_assigned_taxa_ph": summary_data.get("number_of_nt_assigned_taxa", "N/A"),
        "py_number_of_nr_assigned_taxa_ph": summary_data.get("number_of_nr_assigned_taxa", "N/A"),
        "py_taxa_id_total_ph": summary_data.get("taxa_id_total", "N/A"),
        "py_taxa_id_sample_ph": summary_data.get("taxa_id_sample", "N/A"),
        "py_taxa_id_negative_ph": summary_data.get("taxa_id_negative", "N/A"),
        "py_taxa_id_sample_and_negative_ph": summary_data.get("taxa_id_sample_and_negative", "N/A"),
    }

    # Replace placeholders in HTML
    for placeholder, value in replace_dict.items():
        html_content = html_content.replace(placeholder, value)

    # Write the updated HTML
    output_file = dirpath + "/pipeline_summary.html"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

if __name__ == "__main__":
    gen_pipeline_summary(sys.argv[1])