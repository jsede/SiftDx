nextflow.enable.dsl=2

include { fastp; prinseq; kraken2 } from './modules/preprocessing.nf'

workflow {    
    def output = params.output ?: new File(params.r1).parent
    def pair_id = new File(params.r1).name.split('_')[0]
    def reads = [file(params.r1), file(params.r2)]
    def kdb = params.kraken2_db
    
    // Call fastp
    fastp_data = fastp(
        pair_id, 
        reads, 
        output
    )
    
    // Call prinseq
    prinseq_data = prinseq(
        fastp_data,
        output
     )

    // Call Kraken2
    kraken2_data = kraken2(
        prinseq_data,
        kdb,
        output
    )
}