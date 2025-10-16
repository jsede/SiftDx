process gen_zscore_table {
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "table_summary.html"
    
    input:
        val pair_id
        tuple path(zscore),
            path(detected_pathogens),
            path(summary)
        path(negative)
        path table_template
        val output

    output:
        path "table_summary.html"

    script:
    """
    ${params.python} ${baseDir}/scripts/gen_table_summary.py ${zscore} ${negative}
    """
}

process gen_rpm_table {
    publishDir "${output}/${pair_id}/summary", mode: 'copy', pattern: "table_summary.html"
    
    input:
        val pair_id
        tuple path(zscore),
            path(detected_pathogens),
            path(summary)
        path table_template
        val output

    output:
        path "table_summary.html"

    script:
    """
    ${params.python} ${baseDir}/scripts/gen_table_summary.py ${zscore} "None"
    """
}

workflow table_summary {
    take:
        pair_id
        taxonomy_data
        negative
        table_template
        output

    main:
        if (negative) {
            gen_zscore_table(
                pair_id,
                taxonomy_data,
                negative,
                table_template,
                output
            )
        } else {
            gen_rpm_table(
                pair_id,
                taxonomy_data,
                table_template,
                output
            )
        }
}