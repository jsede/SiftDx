process gen_summary {
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "summary.txt"
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "pipeline_summary.html"

    input:
        val pair_id
        tuple path(final_decisions),
            path(full_read_contig_info),
            path(zscore_input),
            path(fqc_txt)
        tuple path(zscore),
            path(detected_pathogens)
        val pipeline_template
        val output

    output:
        tuple path("pipeline_summary.html"),
            path(fqc_txt)



    script:
    """
    python ${baseDir}/scripts/gen_summary_html.py ${fqc_txt} ${pipeline_template}
    """
}

workflow pipeline_summary {
    take:
        pair_id
        finalisation_data
        zscore_data
        pipeline_template
        output

    main:
        gen_summary(
            pair_id,
            finalisation_data,
            zscore_data,
            pipeline_template,
            output
        )
    
}