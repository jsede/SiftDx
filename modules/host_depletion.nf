process kraken2 {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "kraken2.nonhost.*"
    
    input:
        val pair_id
        tuple path(read1), path(read2), path(single1), path(single2), path(bad1), path(bad2), path(fqc_txt)
        val output

    output:
        tuple path("kraken2.nonhost.report"),
        path("kraken2.nonhost.output"),
        path("kraken2_nonhost_1.fastq"),
        path("kraken2_nonhost_2.fastq"),
        path("kraken2_summary.txt")

    script:
    """
    kraken2 --db ${params.kraken2_db} \
            --report kraken2.nonhost.report \
            --output kraken2.nonhost.output \
            --use-mpa-style \
            --paired ${read1} ${read2} \
            --unclassified-out kraken2_nonhost#.fastq

    awk '{s++} END {printf "kraken2_host_depleted: %.0f\\n", (s/4)*2}' kraken2_nonhost_1.fastq > kraken2_summary.txt

    """
}

process bowtie2 {
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "kraken2_nonhost_*"

    input:
        val pair_id
        tuple path(kraken2report), path(kraken2output),path(read1), path(read2), path(fqc_txt)
        path bowtie2_index
        val output

    output:
        tuple path("bowtie2_host.sam"),
        path("kraken2_nonhost_1.fastq.gz"),
        path("kraken2_nonhost_2.fastq.gz")

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
    rm kraken2_nonhost_*.fastq
    cp -i \${read1_resolved} .
    cp -i \${read2_resolved} .

    # Gzip the original files
    gzip \$(basename \${read1_resolved})
    gzip \$(basename \${read2_resolved})
    
    """
}

process samtools {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "{host_depleted*,bowtie2_host_sorted.bam}"

    input:
        val pair_id
        tuple path(bowtie2_output), path(read1), path(read2)
        val output
    
    output:
        tuple path("host_depleted_1.fastq.gz"),
        path("host_depleted_2.fastq.gz"),
        path("bowtie2_host_sorted.bam"),
        path("bowtie2_summary.txt")

    
    script:
    """
    samtools fastq \
        -1 host_depleted_1.fastq \
        -2 host_depleted_2.fastq \
        -0 /dev/null -s /dev/null -n -f 13 \
        ${bowtie2_output}
    
    awk '{s++} END {printf "bowtie2_host_depleted: %.0f\\n", (s/4)*2}' host_depleted_1.fastq > bowtie2_summary.txt
    gzip host_depleted_1.fastq
    gzip host_depleted_2.fastq

    samtools sort ${bowtie2_output} -o bowtie2_host_sorted.bam
    """
}

process ercc {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "*.txt, ercc_plot.png"
    
    input:
        val pair_id
        tuple path(read1), path(read2), path(bowtie2_sorted), path(fqc_txt)
        val output
    
    output:
        tuple path("bowtie2_coverage.txt"),
        path("ercc_coverage.txt"),
        path("ercc_plot.png"),
        path("spikein_summary.txt")

    
    script:
    """
    pileup.sh in=${bowtie2_sorted} out=bowtie2_coverage.txt -Xmx8g secondary=false
    egrep -e "^ERCC" bowtie2_coverage.txt > ercc_coverage.txt
    samtools index ${bowtie2_sorted}
    regions=()
    while IFS=\$'\t' read -r region _; do
        regions+=("\$region")
    done < ercc_coverage.txt
    region_list=\$(IFS=','; echo "\${regions[*]}")
    samtools view -F 260 ${bowtie2_sorted} \${region_list} -o bowtie2_host_Aligned.ERCC_only.out.sam
    ${params.python} ${baseDir}/scripts/ercc_plot.py ercc_plot.png ${params.ercc_config} ercc_coverage.txt
    ercc_count=`cut -f 1 bowtie2_host_Aligned.ERCC_only.out.sam | sort -k 1 | uniq | wc -l`
    echo "ercc_reads: \$ercc_count" > spikein_summary.txt
    """
}

process sequins {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "*.txt, sequins_plot.png"
    
    input:
        val pair_id
        tuple path(read1), path(read2), path(bowtie2_sorted), path(fqc_txt)
        val output
    
    output:
        tuple path("bowtie2_coverage.txt"),
        path("sequins_coverage.txt"),
        path("sequins_plot.png"),
        path("spikein_summary.txt")

    
    script:
    """
    pileup.sh in=${bowtie2_sorted} out=bowtie2_coverage.txt -Xmx8g secondary=false
    egrep -e "^Sequins" bowtie2_coverage.txt > sequins_coverage.txt
    samtools index ${bowtie2_sorted}
    regions=()
    while IFS=\$'\t' read -r region _; do
        regions+=("\$region")
    done < sequins_coverage.txt
    region_list=\$(IFS=','; echo "\${regions[*]}")
    samtools view -F 260 ${bowtie2_sorted} \${region_list} -o bowtie2_host_Aligned.sequins_only.out.sam
    ${params.python} ${baseDir}/scripts/ercc_plot.py sequins_plot.png ${params.sequins_config} sequins_coverage.txt
    sequins_count=`cut -f 1 bowtie2_host_Aligned.sequins_only.out.sam | sort -k 1 | uniq | wc -l`
    echo "sequins_reads: \$sequins_count" > spikein_summary.txt
    """
}

process sortmerna {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "sortmerna_output"

    input:
        val pair_id
        tuple path(host_depleted_1), path(host_depleted_2), path(bowtie2_sorted), path(fqc_txt)
        val output

    output:
        tuple path("fullyQc_1.fastq.gz"),
        path("fullyQc_2.fastq.gz")

    script:
    """
    mkdir -p ${output}/${pair_id}/preprocessing/sortmerna
    sortmerna \
        --ref ${params.sortmerna_db} \
        --aligned ${output}/${pair_id}/aligned \
        --other fullyQc \
        --fastx \
        --reads ${host_depleted_1} --reads ${host_depleted_2} \
        --out2 TRUE \
        --paired_in TRUE \
        --workdir ${output}/${pair_id}/preprocessing/sortmerna


    gzip fullyQc_*.fastq
    """
}

process fullyqc {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "fullyqc_output"

    input:
        val pair_id
        tuple path(host_depleted_1), path(host_depleted_2), path(bowtie2_sorted), path(fqc_txt)
        val output

    output:
        tuple path("fullyQc_1.fastq.gz"),
        path("fullyQc_2.fastq.gz")

    script:
    """
    cp ${host_depleted_1} fullyQc_1.fastq.gz
    cp ${host_depleted_2} fullyQc_2.fastq.gz
    
    """
}

process merged_hd_summaries {
    input:
        tuple path(kraken2report), path(kraken2output),path(read1), path(read2), path(kraken2_summary)
        tuple path(host_depleted_1), path(host_depleted_2), path(bowtie2_sorted), path(bowtie2_summary)

    output:
        path "host_depletion_summary.txt"

    script:
    """
    cat ${kraken2_summary} ${bowtie2_summary} > host_depletion_summary.txt
    """
}

process merged_hd_ercc_summaries {
    input:
        tuple path(kraken2report), path(kraken2output),path(read1), path(read2), path(kraken2_summary)
        tuple path(host_depleted_1), path(host_depleted_2), path(bowtie2_sorted), path(bowtie2_summary)
        tuple path(bowtie2_coverage), path(ercc_coverage), path(ercc_plot), path(spikein_summary)

    output:
        path "host_depletion_summary.txt"

    script:
    """
    cat ${kraken2_summary} ${bowtie2_summary} ${spikein_summary}> host_depletion_summary.txt
    """
}

workflow host_depletion {
    take:
        pair_id
        prinseq_data  // the list of reads (R1 and R2)
        qc_summary // the qc summary file
        bowtie2_index // bowtie2 index for host depletion
        ercc_config  // ERCC configuration file (optional)
        sequins_config // Sequins configuration file (optional)
        output  // the output directory

    main:
        if (!params.na || !['DNA', 'RNA'].contains(params.na.toUpperCase())) {
            error "Invalid value for --na. Must be either 'DNA' or 'RNA'."
        }

        log.info "Starting preprocessing:"
        log.info "Pair ID: ${pair_id}"
        log.info "Output: ${output}"
        log.info "Kraken2 database: ${params.kraken2_db}"
        log.info "Bowtie2 index: ${bowtie2_index}"
        log.info "ERCC config: ${params.ercc_config}"
        log.info "Sequins config: ${params.sequins_config}"
        

        // Call Kraken2
        kraken2_data = kraken2(
            pair_id,
            prinseq_data,
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
            spikein_data = ercc(
                pair_id,
                samtools_data,
                output
            )
        } else if (sequins_config) {
            // Call Sequins
            spikein_data = sequins(
                pair_id,
                samtools_data,
                output
            )
        } else {
            log.info "Skipping ERCC/Sequins process as neither --ercc nor --sequins is provided."
            spikein_data = Channel.empty()
        }
        
        if (params.na?.toUpperCase() == 'RNA') {
            // Call SortMeRNA
            log.info "Temporarily skipping Sortmerna as software for macos is not available."
            fullyqc_data = sortmerna(
                pair_id,
                samtools_data,
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

        if (ercc_config || sequins_config) {
            merge_hd_summaries = merged_hd_ercc_summaries(
                kraken2_data,
                samtools_data,
                spikein_data
            )
        } else {
            merge_hd_summaries = merged_hd_summaries(
                kraken2_data,
                samtools_data
            )
        }
        
        


    emit:
        fullyqc_data = fullyqc_data
        qc_summary = qc_summary
        hd_summary = merge_hd_summaries
        spikein_data = spikein_data
        
}