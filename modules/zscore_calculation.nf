process gen_zscore {
    publishDir "${output}/${pair_id}/results", mode: 'copy', pattern: "zscore.tsv"
    
    input:
        val pair_id
        path (negative, stageAs: 'negative')
        tuple path(final_decisions),
            path(full_read_contig_info),
            path(zscore_input),
            path(fqc_txt)
        path(preprocessing_summary)
        val output

    output:
        tuple path("zscore.tsv"),
        path("detected_pathogens.tsv"),
        path("zscore_summary.txt")

    script:
    """
    ${params.python} ${baseDir}/scripts/zscore.py ${zscore_input} ${negative} ${preprocessing_summary}
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
    ' zscore.tsv > zscore_summary.txt

    """
    }

process no_negative {
    input:
        val pair_id
        tuple path(final_decisions),
            path(full_read_contig_info),
            path(zscore_input),
            path(fqc_txt)
        path(preprocessing_summary)
        val output
    
    output:
        tuple path("rpm_only.tsv"),
        path("detected_pathogens.tsv"),
        path("zscore_summary.txt")
    
    script:
    """
    ${params.python} ${baseDir}/scripts/calc_rpm_only.py ${zscore_input} ${preprocessing_summary}
    awk -F'\t' '
        NR == 1 { next }
        { total++ }
        END {
            printf "taxa_id_total: %d\\n", total
            printf "taxa_id_sample: %d\\n", total
            printf "taxa_id_negative: %d\\n", 0
            printf "taxa_id_sample_and_negative: %d\\n", 0
        }
    ' "${zscore_input}" > zscore_summary.txt
    """
}

workflow zscore_calculation {
    take:
        pair_id
        negative
        finalisation_data
        preprocessing_summary
        output

    main:
        if (negative) {
            zscore_data = gen_zscore(
                pair_id,
                negative,
                finalisation_data,
                preprocessing_summary,
                output
            )
        } else {
            zscore_data = no_negative (
                pair_id,
                finalisation_data,
                preprocessing_summary,
                output
            )
        }

    emit:
        zscore_data
}