process fastp {
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy'
    
    input:
        val pair_id
        path reads  // Accepts the list of reads (R1 and R2)
        val output  // Accept the output directory

    output:
        tuple path("fastp_1.fastq.gz"), 
        path("fastp_2.fastq.gz"), 
        path("fastp.json")    

    script:
    """
    fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 fastp_1.fastq \
        --out2 fastp_2.fastq \
        --json fastp.json \
        -b ${params.read_length} -B ${params.read_length} \
        --qualified_quality_phred ${params.qualified_quality_phred} \
        --unqualified_percent_limit ${params.unqualified_percent_limit} \
        --length_required ${params.length_required} \
        --low_complexity_filter \
        --detect_adapter_for_pe \
        --thread ${params.threads}
    
    gzip fastp_1.fastq
    gzip fastp_2.fastq
    """
}

process prinseq {
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "prinseq_lc_good*"

    input:
        val pair_id
        tuple path(read1), path(read2), path(json)
        val output

    output:
        tuple path("prinseq_lc_good_out_R1.fastq.gz"),
        path("prinseq_lc_good_out_R2.fastq.gz"),
        path("prinseq_lc_single_out_R1.fastq.gz"),
        path("prinseq_lc_single_out_R2.fastq.gz"),
        path("prinseq_lc_bad_out_R1.fastq.gz"),
        path("prinseq_lc_bad_out_R2.fastq.gz")

    script:
    """
    prinseq++ \
        -fastq ${read1} \
        -fastq2 ${read2} \
        -lc_dust 0.07 \
        -out_name prinseq_lc

    gzip prinseq*.fastq
    """
}

workflow qc{
    take:
        pair_id
        reads  // the list of reads (R1 and R2)
        output  // the output directory
    
    main:
        // Call fastp
        fastp_data = fastp(
            pair_id, 
            reads, 
            output
        )
        
        // Call prinseq
        prinseq_data = prinseq(
            pair_id,
            fastp_data,
            output
        )
        
    emit:
        prinseq_data
}