process merge_nt_nr {
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "summary.txt"
    publishDir "${output}/${pair_id}/analysis", mode: 'copy', pattern: "final_decisions.tsv"
    publishDir "${output}/${pair_id}/analysis", mode: 'copy', pattern: "full_read_contig_info.tsv"
    publishDir "${output}/${pair_id}/analysis", mode: 'copy', pattern: "zscore_input.tsv"

    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats),
            path(fqc_txt)
        tuple path(mm2_contigs), // minimap2_contig_out.paf
            path(mm2_contigs_m8), // minimap2_contig_out_frompaf.m8
            path(mm2_reads), // minimap2_reads_out.paf
            path(mm2_reads_m8), // minimap2_reads_out_frompaf.m8
            path(k2_pluspf), // combined_rc.kraken.txt
            path(nr_alignments), // nr_alignments_file.tsv
            path(nt_alignments) //n t_alignments_sr_blast.tsv
        val entrez_email
        val entrez_api_key
        val output
        

    output:
        tuple path("final_decisions.tsv"),
        path("full_read_contig_info.tsv"),
        path("zscore_input.tsv"),
        path("merge_nt_nr_summary.txt")

    script:
    """
    cat <<EOF > inputs.txt
    combined_lr_contigs_fa ${combined_lr_contigs_fa}
    cov_stats ${cov_stats}
    minimap2_contig_out ${mm2_contigs_m8}
    minimap2_reads_out ${mm2_reads_m8}
    kraken_pluspf_file ${k2_pluspf}
    nr_alignments_file ${nr_alignments}
    nt_alignments_sr_file ${nt_alignments}
    EOF

    ${params.python} ${baseDir}/scripts/merge_nt_nr.py ${params.database} ${params.taxdump} inputs.txt ${entrez_email} ${entrez_api_key}

    awk 'NR > 1 && \$1 != "Unclassified" { count++ } END { print "number_of_identified_taxa:", count }' zscore_input.tsv > temp_id_taxa.txt
    awk -F'\\t' 'NR == 1 {for (i=1; i<=NF; i++) if (\$i=="Fasta_Headers") h=i; else if (\$i=="final") f=i}
        NR > 1 && \$h ~ /^k14/ {
            if (\$f == "-") unassigned++;
            else assigned++;
        }
        END {
            print "number_of_contigs_unassigned:", unassigned ? unassigned : 0;
            print "number_of_contigs_assigned:", assigned ? assigned : 0;
        }' full_read_contig_info.tsv > temp_contig_stats.txt
    awk -F'\\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if (\$i == "Fasta_Headers") h = i;
            else if (\$i == "final") f = i;
        }
    }
    NR > 1 && \$h !~ /^k14/ {
        if (\$f == "-") unassigned++;
        else assigned++;
    }
    END {
        print "number_of_reads_unassigned:", unassigned ? unassigned : 0;
        print "number_of_reads_assigned:", assigned ? assigned : 0;
    }' full_read_contig_info.tsv > temp_read_stats.txt
    awk -F'\\t' '
    NR==1 {
        for (i=1; i<=NF; i++) {
            if (\$i == "Kraken_Species") kraken_col = i
            if (\$i == "NT_Species") nt_col = i
            if (\$i == "NR_Species") nr_col = i
        }
        next
    }
    {
        if (kraken_col && \$kraken_col != "-") kraken[\$kraken_col] = 1
        if (nt_col && \$nt_col != "-") nt[\$nt_col] = 1
        if (nr_col && \$nr_col != "-") nr[\$nr_col] = 1
    }
    END {
        printf "number_of_kraken_assigned_taxa: %d\\n", length(kraken)
        printf "number_of_nt_assigned_taxa: %d\\n", length(nt)
        printf "number_of_nr_assigned_taxa: %d\\n", length(nr)
    }
    ' full_read_contig_info.tsv > temp_tc_stats.txt

    cat temp_id_taxa.txt temp_contig_stats.txt temp_read_stats.txt temp_tc_stats.txt > merge_nt_nr_summary.txt
    """
}

workflow merge_taxonomy {
    take:
        pair_id
        preprocessing_data
        tax_class_data
        entrez_email
        entrez_api_key
        output
    
    main:
        merge_data = merge_nt_nr(
            pair_id,
            preprocessing_data,
            tax_class_data,
            entrez_email,
            entrez_api_key,
            output
        )

    emit:
        merge_data
}