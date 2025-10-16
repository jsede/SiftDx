process gen_summary {
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "summary.txt"
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "pipeline_summary.html"

    input:
        val pair_id
        val summary_data
        val pipeline_template
        val output

    output:
        tuple path("pipeline_summary.html"),
            path("summary.txt")



    script:
    """
    python3 ${baseDir}/scripts/gen_summary_html.py ${summary_data}
    cp ${summary_data} .
    """
}

workflow pipeline_summary {
    take:
        pair_id
        summary_data
        pipeline_template
        output

    main:
        gen_summary(
            pair_id,
            summary_data,
            pipeline_template,
            output
        )
    
}