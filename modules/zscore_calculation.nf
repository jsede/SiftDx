include { table_summary } from './table_summary.nf'

process gen_zscore {
    publishDir "${output}/${pair_id}/results", mode: 'copy', pattern: "final_decisions.tsv"
    publishDir "${output}/${pair_id}/results", mode: 'copy', pattern: "zscore.tsv"
    input:
        val pair_id
        path (negative, stageAs: 'negative')
        tuple path(final_decisions),
            path(full_read_contig_info),
            path(zscore_input),
            path(fqc_txt)
        val output

    output:
        tuple path("zscore.tsv"),
        path("detected_pathogens.tsv")

    script:
    """
    python ${baseDir}/scripts/zscore.py ${zscore_input} ${negative} ${fqc_txt}
    """
    }

workflow zscore_calculation {
    take:
        pair_id
        negative
        finalisation_data
        table_template
        output

    main:
        zscore_data = gen_zscore(
            pair_id,
            negative,
            finalisation_data,
            output
        )

        table_summary(
            pair_id,
            zscore_data,
            table_template,
            output
        )

}