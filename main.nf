nextflow.enable.dsl=2

include { preprocessing } from './modules/preprocessing.nf'

workflow {    
    def output = params.output ?: new File(params.r1).parent
    def pair_id = new File(params.r1).name.split('_')[0]
    def reads = [file(params.r1), file(params.r2)]
    def kdb = params.kraken2_db
    def bowtie2_index = params.use_ercc ? params.bowtie2_index_ercc : params.bowtie2_index
    
    // Call the preprocessing workflow
    preprocessing_data = preprocessing(
        pair_id,
        reads,
        kdb,
        bowtie2_index,
        output
    )

}