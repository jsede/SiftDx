process megahit {
    input:
        val pair_id
        tuple path(fullyQc_1), path(fullyQc_2)
        val output

    output:
        path("megahit/final.contigs.fa")

    publishDir "${output}/${pair_id}", mode: 'copy'
    
    script:
    """
    megahit \
        -1 ${fullyQc_1} \
        -2 ${fullyQc_2} \
        -o megahit
    """
}

process unassembled_reads {
    input:
        val pair_id
        path megahit_contigs
        tuple path(fullyQc_1), path(fullyQc_2)
        val output

    output:
        tuple path("unassembled_reads_fwd.fq.gz"),
            path("unassembled_reads_rev.fq.gz"),
            path("reads_mapped_to_contigs.cov_stats")
        path("final.contigs.fq")
    
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "unassembled_reads*"
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "reads_mapped_to_contigs.*"

    script:
    """
    bbwrap.sh ref=${megahit_contigs} \
        in1=${fullyQc_1} \
        in2=${fullyQc_2} \
        out=reads_mapped_to_contigs_unsorted.sam \
        -Xmx=8g
    samtools sort -n -o reads_mapped_to_contigs.sam reads_mapped_to_contigs_unsorted.sam
    samtools fastq -f 12 -1 unassembled_reads_fwd.fq \
        -2 unassembled_reads_rev.fq \
        reads_mapped_to_contigs.sam
    samtools view -bS reads_mapped_to_contigs.sam > reads_mapped_to_contigs.bam
    samtools sort reads_mapped_to_contigs.bam > reads_mapped_to_contigs.sorted.bam
    samtools coverage reads_mapped_to_contigs.sorted.bam > reads_mapped_to_contigs.cov_stats
    seqtk seq -F '#' ${megahit_contigs} > final.contigs.fq

    gzip unassembled_reads_*.fq
    """
}

process megahit_fail {
    input:
        val pair_id
        path megahit_contigs
        tuple path(fullyQc_1), path(fullyQc_2)
        val cov_stats
        val output

    output:
        tuple path("unassembled_reads_fwd.fq.gz"),
            path("unassembled_reads_rev.fq.gz"),
            path("reads_mapped_to_contigs.cov_stats")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "unassembled_reads*"
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "reads_mapped_to_contigs.cov_stats"

    script:
    """
    cp ${fullyQc_1} unassembled_reads_fwd.fq.gz
    cp ${fullyQc_2} unassembled_reads_rev.fq.gz
    cp ${cov_stats} reads_mapped_to_contigs.cov_stats
    """
}

process alignment_prep {
    input:
        val pair_id
        path megahit_contigs
        
        tuple path(unassembled_reads_fwd), 
            path(unassembled_reads_rev),
            path(cov_stats)
        val output

    output:
        tuple path("combined_sr_file.fq"),
            path("combined_sr_file.fa"),
            path("combined_reads_contigs_file.fq"),
            path("combined_reads_contigs_file.fa")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy'
    
    script:
    """
    python ${baseDir}/scripts/separate_reads_by_size.py ${unassembled_reads_fwd} unassembled_reads_longer_fwd.fq unassembled_reads_shorter_fwd.fq
    python ${baseDir}/scripts/separate_reads_by_size.py ${unassembled_reads_rev} unassembled_reads_longer_rev.fq unassembled_reads_shorter_rev.fq 

    seqtk mergepe unassembled_reads_shorter_fwd.fq unassembled_reads_shorter_rev.fq > combined_sr_file.fq
    seqtk seq -a combined_sr_file.fq > combined_sr_file.fa
    seqtk mergepe ${unassembled_reads_fwd} ${unassembled_reads_rev} > merged_reads.fq
    cat ${megahit_contigs} merged_reads.fq > combined_reads_contigs_file.fq
    seqtk seq -a combined_reads_contigs_file.fq > combined_reads_contigs_file.fa
    """
}

workflow assembly {
    take:
        pair_id
        fully_qc
        cov_stats
        output

    main:
        // Run megahit process
        megahit(
            pair_id, 
            fully_qc, 
            output
        )

        valid_files = megahit.out.filter { file -> file.size() > 0 }
        invalid_files = megahit.out.filter { file -> file.size() == 0 }
        
        valid_files.ifEmpty { 
            log.info "No valid files found. Skipping unassembled_reads process."
        }

        // Handle empty invalid_files channel
        invalid_files.ifEmpty { 
            log.info "Valid files found. Skipping megahit_fail process."
        }
        
        unassembled_reads(
            pair_id,
            valid_files,
            fully_qc,
            output
        )

        megahit_fail(
            pair_id,
            invalid_files,
            fully_qc,
            cov_stats,
            output
        )
        
        final_reads = Channel.empty()
        final_reads = final_reads.mix(unassembled_reads.out[0], megahit_fail.out[0])
        final_reads.view()

        alignment_prep(
            pair_id,
            megahit.out,
            final_reads,
            output
        )


    emit:
        alignment_prep.out
        

}