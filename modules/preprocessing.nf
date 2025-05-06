process fastp {
    //tag "pair_id=${pair_id} reads=${reads} output=${output}"

    input:
        val(pair_id)
        path(reads)  // Accepts the list of reads (R1 and R2)
        val(output)  // Accept the output directory

    output:
        tuple val(pair_id), path("fastp_1.fastq.gz"), path("fastp_2.fastq.gz"), path("fastp.json")

    script:
    """
    mkdir -p ${output}/${pair_id}  

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
    //tag "pair_id=${pair_id}"

    input:
        tuple val(pair_id), path(read1), path(read2), path(json)
        val(output)


    output:
        tuple val(pair_id),
        path("prinseq_lc_good_out_R1.fastq.gz"),
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

process kraken2 {
    input:
        tuple val(pair_id), path(read1), path(read2), path(single1), path(single2), path(bad1), path(bad2)
        path(kdb)
        val(output)


    output:
        tuple path("kraken2.nonhuman.report"),
        path("kraken2.nonhuman.output"),
        path("kraken2_nonhuman_1.fastq.gz"),
        path("kraken2_nonhuman_2.fastq.gz")

    script:
    """
    kraken2 --db ${kdb} \
            --report kraken2.nonhuman.report \
            --output kraken2.nonhuman.output \
            --use-mpa-style \
            --paired ${read1} ${read2} \
            --unclassified-out kraken2_nonhuman#.fastq

    gzip kraken2_nonhuman_*.fastq

    """
}


