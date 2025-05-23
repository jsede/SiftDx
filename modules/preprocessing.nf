include { qc } from './qc.nf'
include { human_depletion } from './human_depletion.nf'
include { assembly } from './assembly.nf'

workflow preprocessing {
    take:
        pair_id
        reads  // the list of reads (R1 and R2)
        kdb    // Kraken2 database
        bowtie2_index  // Bowtie2 index
        ercc_config  // ERCC configuration file (optional)
        sortmerna_db  // SortMeRNA database
        cov_stats // Coverage statistics file
        output  // the output directory

    main:
        qc_data = qc(
            pair_id,
            reads,
            output
        )
        
        human_depletion_data = human_depletion(
            pair_id,
            qc_data,
            kdb,
            bowtie2_index,
            ercc_config,
            sortmerna_db,
            cov_stats,
            output
        )

        assembly_data = assembly(
            pair_id,
            human_depletion_data,
            cov_stats,
            output
        )

    emit:
        assembly_data
}
