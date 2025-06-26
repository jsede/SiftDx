process gen_table {
    publishDir "${output}/${pair_id}/summary", mode: 'copy'
    
    input:
        val pair_id
        path zscore
        path table_template
        val output

    output:
        path "table_summary.html"

    script:
    """
    python ${baseDir}/scripts/gen_table_summary.py ${zscore} ${table_template}
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