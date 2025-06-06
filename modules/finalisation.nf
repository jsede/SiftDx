include { merge_taxonomy } from './merge_taxonomy.nf'
include { table_summary } from './table_summary.nf'

workflow finalisation {
    take:
        pair_id
        preprocessing_data
        tax_class_data
        database
        taxdump
        output
        entrez_email
        entrez_api_key
        table_template
        

    main:
        taxonomy_data = merge_taxonomy (
            pair_id,
            preprocessing_data,
            tax_class_data,
            database,
            taxdump,
            entrez_email,
            entrez_api_key
        )

        table_summary (
            pair_id,
            taxonomy_data,
            table_template
        )

    // emit:
    //     taxonomy_data
}