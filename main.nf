nextflow.enable.dsl=2

include { preprocessing } from './modules/preprocessing.nf'
include { taxonomic_classification } from './modules/taxonomic_classification.nf'

workflow {    
    def output = params.output ?: new File(params.r1).parent
    def pair_id = new File(params.r1).name.split('_')[0]
    def reads = [file(params.r1), file(params.r2)]
    def kdb = params.kraken2_db // Kraken2 HUMAN database
    def bowtie2_index = params.ercc ? params.bowtie2_index_ercc : params.bowtie2_index // if ERCC is flag provided, set bowtie2_index to params.bowtie2_index_ercc, else params.bowtie2_index
    def ercc_config = params.ercc ? params.ercc_config : null // if ERCC is flag provided, set ercc_config to params.ercc_config, else null
    def sortmerna_db = params.na == 'RNA' ? params.sortmerna_db : null //if NA is RNA, set sortmerna_db to params.sortmerna_db, else null
    def cov_stats = params.cov_stats
    def mm2_index = params.mm2_index
    def pluspf_db = params.pluspf_db // Kraken2 PLUSPF database
    def diamond_db = params.diamond_db // Diamond database
    def blast_db = params.blast_db // BLAST database
    
    // Call the preprocessing workflow
    preprocessing_data = preprocessing(
        pair_id,
        reads,
        kdb,
        bowtie2_index,
        ercc_config,
        sortmerna_db,
        cov_stats,
        output
    )
    
    
    // Call the taxonomic classification workflow
    tax_class_data = taxonomic_classification(
        pair_id,
        preprocessing_data,
        mm2_index,
        pluspf_db,
        diamond_db,
        blast_db,
        output
    )
}