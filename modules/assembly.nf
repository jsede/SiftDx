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
        tuple path("unassembled_reads_fwd.fq"),
        path("unassembled_reads_rev.fq"),
        path("reads_mapped_to_contigs.cov_stats")
    
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
        tuple path("unassembled_reads_fwd.fq"),
        path("unassembled_reads_rev.fq"),
        path("reads_mapped_to_contigs.cov_stats")

    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "unassembled_reads*"
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "reads_mapped_to_contigs.cov_stats"

    script:
    """
    cp ${fullyQc_1} unassembled_reads_fwd.fq
    cp ${fullyQc_2} unassembled_reads_rev.fq
    cp ${cov_stats} reads_mapped_to_contigs.cov_stats
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
        
        final_output = Channel.empty()
        final_output = final_output.mix(unassembled_reads.out, megahit_fail.out)
    
    emit:
        final_output

}