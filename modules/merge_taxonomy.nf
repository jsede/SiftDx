process merge_nt_nr {
    input:
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats)
        tuple path(mm2_contigs), // minimap2_contig_out.paf
            path(mm2_contigs_m8), // minimap2_contig_out_frompaf.m8
            path(mm2_reads), // minimap2_reads_out.paf
            path(mm2_reads_m8), // minimap2_reads_out_frompaf.m8
            path(k2_pluspf), // combined_rc.kraken.txt
            path(nr_alignments), // nr_alignments_file.tsv
            path(nt_alignments) //n t_alignments_sr_blast.tsv
        path database
        path taxdump
        val entrez_email
        val entrez_api_key
        

    output:
        tuple path("final_decisions.csv"),
        path("full_read_contig_info.tsv"),
        path("zscore_input.tsv")

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

    python ${baseDir}/scripts/merge_nt_nr.py ${database} ${taxdump} inputs.txt ${entrez_email} ${entrez_api_key}
    rm inputs.txt
    """
}

workflow merge_taxonomy {
    take:
        pair_id
        preprocessing_data
        tax_class_data
        database
        taxdump
        entrez_email
        entrez_api_key
    
    main:
        merge_data = merge_nt_nr(
            preprocessing_data,
            tax_class_data,
            database,
            taxdump,
            entrez_email,
            entrez_api_key
        )

    emit:
        merge_data
}