include { qc } from './qc.nf'
include { host_depletion } from './host_depletion.nf'
include { assembly } from './assembly.nf'

process merged_preprocessing_summaries {
    input:
        path qc_summary
        path hd_summary
        path assembly_summary

    output:
        path "preprocessing_summary.txt"

    script:
    """
    cat ${qc_summary} ${hd_summary} ${assembly_summary} > preprocessing_summary.txt
    """
}
workflow preprocessing {
    take:
        pair_id
        reads  // the list of reads (R1 and R2)
        kdb    // Kraken2 database
        bowtie2_index  // Bowtie2 index
        ercc_config  // ERCC configuration file (optional)
        sequins_config // Sequins configuration file (optional)
        sortmerna_db  // SortMeRNA database
        cov_stats // Coverage statistics file
        output  // the output directory

    main:
        qc_data = qc(
            pair_id,
            reads,
            output
        )
        
        host_depletion_data = host_depletion(
            pair_id,
            qc_data.prinseq_data, // Use prinseq_data output from qc
            qc_data.qc_summary, // Use qc_summary output from qc
            kdb,
            bowtie2_index,
            ercc_config,
            sequins_config,
            sortmerna_db,
            cov_stats,
            output
        )

        assembly_data = assembly(
            pair_id,
            host_depletion_data.fullyqc_data, // Use fullyqc_data output from host_depletion
            cov_stats,
            output
        )

        preprocessing_summary = merged_preprocessing_summaries(
            qc_data.qc_summary,
            host_depletion_data.hd_summary,
            assembly_data.assembly_summary
        )
        
    emit:
        aln_prep_data = assembly_data.aln_prep_data
        preprocessing_summary = preprocessing_summary
}
