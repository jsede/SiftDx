process kraken2 {
    input:
        val pair_id
        tuple path(read1), path(read2), path(single1), path(single2), path(bad1), path(bad2)
        path kdb
        val output

    output:
        tuple path("kraken2.nonhuman.report"),
        path("kraken2.nonhuman.output"),
        path("kraken2_nonhuman_1.fastq"),
        path("kraken2_nonhuman_2.fastq")
    
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "kraken2.nonhuman.*"

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
        val pair_id
        tuple path(kraken2report), path(kraken2output),path(read1), path(read2)
        path bowtie2_index
        val output

    output:
        tuple path("bowtie2_host.sam"),
        path("kraken2_nonhuman_1.fastq.gz"),
        path("kraken2_nonhuman_2.fastq.gz")
    
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "kraken2_nonhuman_*"

    script:
    """
    # Get the absolute paths of the input files
    read1_resolved=\$(readlink -f ${read1})
    read2_resolved=\$(readlink -f ${read2})
    
    idx_base=\$(find -L ${bowtie2_index} -name '*.3.bt2' | head -n 1 | xargs readlink -f | sed 's/\\.3\\.bt2\$//')
    echo "Bowtie2 index base: \${idx_base}"
    bowtie2 -x \${idx_base} \
            --very-sensitive-local \
            -1 ${read1} \
            -2 ${read2} \
            -q -S bowtie2_host.sam

    # Move the files to the current directory
    rm kraken2_nonhuman_*.fastq
    mv \${read1_resolved} .
    mv \${read2_resolved} .

    # Gzip the original files
    gzip \$(basename \${read1_resolved})
    gzip \$(basename \${read2_resolved})
    
    """
}

process samtools {
    input:
        val pair_id
        tuple path(bowtie2_output), path(read1), path(read2)
        val output
    
    output:
        tuple path("host_depleted_1.fastq.gz"),
        path("host_depleted_2.fastq.gz"),
        path("bowtie2_host_sorted.bam")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "host_depleted*"
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "bowtie2_host_sorted.bam"

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

process ercc {
    input:
        val pair_id
        tuple path(read1), path(read2), path(bowtie2_sorted)
        path ercc_config
        val output
    
    output:
        tuple path("bowtie2_coverage.txt"),
        path("ercc_coverage.txt"),
        path("ercc_plot.png")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "*.txt"
    publishDir "${output}/${pair_id}/final_plots", mode: 'copy', pattern: "ercc_plot.png"
    
    script:
    """
    pileup.sh in=${bowtie2_sorted} out=bowtie2_coverage.txt -Xmx8g secondary=false
    egrep -e "^ERCC" bowtie2_coverage.txt > ercc_coverage.txt
    samtools index ${bowtie2_sorted}
    regions=()
    while IFS=\$'\t' read -r region _; do
        regions+=("$region")
    done < ercc_coverage.txt
    region_list=\$(IFS=','; echo "\${regions[*]}")
    samtools view -F 260 ${bowtie2_sorted} \${region_list} -o bowtie2_host_Aligned.ERCC_only.out.sam
    python scripts/ercc_plot.py \${PWD} bowtie2_host_Aligned.ERCC_only.out.sam ercc_coverage.txt
    """
}

process sortmerna {
    input:
        val pair_id
        tuple path(host_depleted_1), path(host_depleted_2), path(bowtie2_sorted) 
        path sortmerna_db
        val output

    output:
        tuple path("host_depleted_1.fastq.gz"),
        path("host_depleted_2.fastq.gz")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "sortmerna_output"

    script:
    """
    mkdir -p ${output}/preprocessing/${pair_id}/sortmerna
    sortmerna \
        --ref ${sortmerna_db} \
        --aligned " + aligned +\
        --other noRna \
        --fastx \
        --reads host_depleted_1 --reads host_depleted_2 \
        --out2 TRUE \
        --paired_in TRUE \
        --workdir ${output}/preprocessing/${pair_id}/sortmerna \
    """
}

process fullyqc {
    input:
        val pair_id
        tuple path(host_depleted_1), path(host_depleted_2), path(bowtie2_sorted) 
        val output

    output:
        tuple path("fullyQc_1.fastq.gz"),
        path("fullyQc_2.fastq.gz")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "fullyqc_output"

    script:
    """
    cp ${host_depleted_1} fullyQc_1.fastq.gz
    cp ${host_depleted_2} fullyQc_2.fastq.gz
    """
}

workflow human_depletion {
    take:
        pair_id
        prinseq_data  // the list of reads (R1 and R2)
        kdb    // Kraken2 database
        bowtie2_index  // Bowtie2 index
        ercc_config  // ERCC configuration file (optional)
        sortmerna_db  // SortMeRNA database
        cov_stats // Coverage statistics file
        output  // the output directory

    main:
        if (!params.na || !['DNA', 'RNA'].contains(params.na.toUpperCase())) {
            error "Invalid value for --na. Must be either 'DNA' or 'RNA'."
        }

        log.info "Starting preprocessing:"
        log.info "Pair ID: ${pair_id}"
        log.info "Output: ${output}"
        log.info "Kraken2 database: ${kdb}"
        log.info "Bowtie2 index: ${bowtie2_index}"
        log.info "ERCC config: ${ercc_config}"
        

        // Call Kraken2
        kraken2_data = kraken2(
            pair_id,
            prinseq_data,
            kdb,
            output
        )

        // Call Bowtie2
        bowtie2_data = bowtie2(
            pair_id,
            kraken2_data,
            bowtie2_index,
            output
        )

        // Call samtools
        samtools_data = samtools(
            pair_id,
            bowtie2_data,
            output
        )
        
        if (ercc_config) {
            // Call ERCC
            ercc_data = ercc(
                pair_id,
                samtools_data,
                ercc_config,
                output // Added missing comma here
            )
        } else {
            log.info "Skipping ERCC process as --ercc is not provided."
        }
        
        if (params.na?.toUpperCase() == 'RNA') {
            // Call SortMeRNA
            fullyqc_data = sortmerna(
                pair_id,
                samtools_data,
                sortmerna_db,
                output
            )
        } else {
            log.info "Skipping SortMeRNA process as --rna is not provided."
            fullyqc_data = fullyqc(
                pair_id,
                samtools_data,
                output
            )
        }
    emit:
        fullyqc_data
}