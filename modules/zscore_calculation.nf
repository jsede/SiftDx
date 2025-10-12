include { table_summary } from './table_summary.nf'
include { pipeline_summary } from './summary_html.nf'

process gen_zscore {
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
    python3 ${baseDir}/scripts/zscore.py ${zscore_input} ${negative} ${fqc_txt}
    awk -F'\t' '
        NR == 1 { next }  # skip header
        { total++ }
        \$9 != "-" { sample++ }
        \$10 != "-" { negative++ }
        \$9 != "-" && \$10 != "-" { both++ }
        END {
            printf "taxa_id_total: %d\\n", total
            printf "taxa_id_sample: %d\\n", sample
            printf "taxa_id_negative: %d\\n", negative
            printf "taxa_id_sample_and_negative: %d\\n", both
        }
    ' zscore.tsv >> ${fqc_txt}

    """
    }

workflow zscore_calculation {
    take:
        pair_id
        negative
        finalisation_data
        table_template
        pipeline_template
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

        pipeline_summary(
            pair_id,
            finalisation_data,
            zscore_data,
            pipeline_template,
            output
        )
}