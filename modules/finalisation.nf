include { merge_taxonomy } from './merge_taxonomy.nf'
include { zscore_calculation } from './zscore_calculation.nf'
include { table_summary } from './table_summary.nf'
include { pipeline_summary } from './summary_html.nf'

process merge_all_summaries {
    input:
        val pair_id
        path preprocessing_summary
        path taxonomy_summary
        path zscore_summary
        path output

    output:
        path "summary.txt"

    script:
    """
    cat ${preprocessing_summary} ${taxonomy_summary} ${zscore_summary} > summary.txt
    """
}

workflow finalisation {
    take:
        pair_id
        preprocessing_data
        preprocessing_summary
        tax_class_data
        negative
        table_template
        pipeline_template
        output
        entrez_email
        entrez_api_key

    main:
        taxonomy_data = merge_taxonomy (
            pair_id,
            preprocessing_data,
            tax_class_data,
            entrez_email,
            entrez_api_key,
            output
        )

        zscore_data = zscore_calculation (
            pair_id,
            negative,
            taxonomy_data,
            preprocessing_summary,
            output
        )

        summary_data = merge_all_summaries(
            pair_id,
            preprocessing_summary,
            taxonomy_data.map { _final, _tsv, _zscore, summary -> summary },
            zscore_data.map { _result, _pathogens, summary -> summary },
            output
        )

        table_summary(
            pair_id,
            zscore_data,
            negative,
            table_template,
            output
        )

        pipeline_summary(
            pair_id,
            summary_data,
            pipeline_template,
            output
        )

    emit:
        taxonomy_data
}