nextflow.enable.dsl=2

include { preprocessing } from './modules/preprocessing.nf'
include { taxonomic_classification } from './modules/taxonomic_classification.nf'
include { finalisation } from './modules/finalisation.nf'
include { zscore_calculation } from './modules/zscore_calculation.nf'


workflow {    
    def output = params.output ?: new File(params.r1).parent
    def pair_id = new File(params.r1).name.split('_S')[0]
    def reads = [file(params.r1), file(params.r2)]
    def negative = params.ndata ?: null
    def bowtie2_index = params.ercc ? params.bowtie2_index_ercc : (params.sequins ? params.bowtie2_index_sequins : params.bowtie2_index) // if ERCC/Sequins is flag provided, set bowtie2_index
    def ercc_config = params.ercc ? params.ercc_config : null // if ERCC is flag provided, set ercc_config to params.ercc_config, else null
    def sequins_config = params.sequins ? params.sequins_config : null // if the sequins flag is provided, set sequis config to params.sequins_config
    def cov_stats = params.cov_stats
    def entrez_email = params.entrez_email
    def entrez_api_key = params.entrez_api_key
    def table_summary = params.table_summary // the table summary html template
    def pipeline_template = params.pipeline_template // the pipeline summary html template

    if (!params.na || !['DNA', 'RNA'].contains(params.na.toUpperCase())) {
        error "Invalid value for --na. Must be either 'DNA' or 'RNA'."
    }


    log.info "Pair ID: ${pair_id}"
    log.info "Negative: ${negative}"
    log.info "Output: ${output}"
    log.info "Kraken2 database: ${params.kraken2_db}"
    log.info "Bowtie2 index: ${bowtie2_index}"
    log.info "ERCC config: ${params.ercc_config}"
    log.info "Sequins config: ${params.sequins_config}"

    // Call the preprocessing workflow
    log.info "Starting preprocessing."
    preprocessing_data = preprocessing(
        pair_id,
        reads,
        bowtie2_index,
        ercc_config,
        sequins_config,
        cov_stats,
        output
    )

    log.info "Preprocessing complete. Starting taxonomic classification."
    // Call the taxonomic classification workflow
    tax_class_data = taxonomic_classification(
        pair_id,
        preprocessing_data.aln_prep_data,
        output
    )
    
    log.info "Taxonomic classification complete. Starting finalisation."
    // Call the finalisation workflow
    finalisation_data = finalisation(
        pair_id,
        preprocessing_data.aln_prep_data,
        preprocessing_data.preprocessing_summary,
        preprocessing_data.spikein_data,
        tax_class_data,
        negative,
        table_summary,
        pipeline_template,
        output,
        entrez_email,
        entrez_api_key
    )

}