include { merge_taxonomy } from './merge_taxonomy.nf'

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

    main:
        taxonomy_data = merge_taxonomy (
            pair_id,
            preprocessing_data,
            tax_class_data,
            database,
            taxdump,
            entrez_email,
            entrez_api_key,
            output
        )

    emit:
        taxonomy_data
}