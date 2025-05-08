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
        path("kraken2_nonhuman_1.fastq"),
        path("kraken2_nonhuman_2.fastq")

    script:
    """
    kraken2 --db ${kdb} \
            --report kraken2.nonhuman.report \
            --output kraken2.nonhuman.output \
            --use-mpa-style \
            --paired ${read1} ${read2} \
            --unclassified-out kraken2_nonhuman#.fastq

    """
}

process bowtie2 {
    input:
        tuple val(pair_id), path(kraken2output), path(read1), path(read2)
        path(bowtie2_index)
        val(output)

    output:
        tuple val(pair_id), path("bowtie2_host.sam")

    script:
    """
    idx_base=\$(find -L ${bowtie2_index} -name '*.3.bt2' | head -n 1 | xargs readlink -f | sed 's/\\.3\\.bt2\$//')
    echo "Bowtie2 index base: \${idx_base}"
    bowtie2 -x \${idx_base} \
            --very-sensitive-local \
            -1 ${read1} \
            -2 ${read2} \
            -q -S bowtie2_host.sam

    """
}

process samtools {
    input:
        tuple val(pair_id), path(bowtie2_output)
    
    output:
        tuple val(pair_id),
        path("host_depleted_1.fastq.gz"),
        path("host_depleted_2.fastq.gz")
        path("bowtie2_host_sorted.bam")

    script:
    """
    samtools fastq \
        -1 host_depleted_1.fastq \
        -2 host_depleted_2.fastq \
        -0 /dev/null -s /dev/null -n -f 13 \
        ${bowtie2_output}
    
    gzip host_depleted_1.fastq
    gzip host_depleted_2.fastq

    samtools sort ${bowtie2_output} -o bowtie2_host_sorted.bam
    """
}

workflow preprocessing {
    take:
        pair_id
        reads  // the list of reads (R1 and R2)
        kdb    // Kraken2 database
        bowtie2_index  // Bowtie2 index
        output  // the output directory

    main:
        log.info "Starting preprocessing:"
        log.info "Pair ID: ${pair_id}"
        log.info "Reads: ${reads}"
        log.info "Output: ${output}"
        log.info "Kraken2 database: ${kdb}"
        log.info "Bowtie2 index: ${bowtie2_index}"
        
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

        // Call Bowtie2
        bowtie2_data = bowtie2(
            kraken2_data,
            bowtie2_index,
            output
        )

        // Call samtools
        samtools_data = samtools(
            bowtie2_data,
        )

    emit:
        kraken2_data
}
