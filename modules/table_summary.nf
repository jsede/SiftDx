process gen_table {
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "table_summary.html"
    
    input:
        val pair_id
        tuple path(zscore),
            path(detected_pathogens)
        path table_template
        val output

    output:
        path "table_summary.html"

    script:
    """
    python3 ${baseDir}/scripts/gen_table_summary.py ${zscore} ${table_template}
    """
}

workflow table_summary {
    take:
        pair_id
        taxonomy_data
        table_template
        output

    main:
        gen_table(
            pair_id,
            taxonomy_data,
            table_template,
            output
        )
    
}