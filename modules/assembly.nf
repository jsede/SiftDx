process megahit {
    cache 'lenient'
    publishDir "${output}/${pair_id}", mode: 'copy', pattern: "megahit/final.contigs.fa"

    input:
        val pair_id
        tuple path(fullyQc_1), path(fullyQc_2), path(fullyqc_summary)
        val output

    output:
        tuple path("megahit/final.contigs.fa"),
            path("megahit_summary.txt")
    
    script:
    """
    megahit \
        -1 ${fullyQc_1} \
        -2 ${fullyQc_2} \
        -o megahit

    line=\$(tail -n 2 megahit/log | head -n 1)
    numContigs=\$(echo "\$line" | awk -F'[, ]+' '{for(i=1;i<=NF;i++) if(\$i=="contigs") print \$(i-1)}')
    min=\$(echo "\$line" | awk -F'[, ]+' '{for(i=1;i<=NF;i++) if(\$i=="min") print \$(i+1)}')
    max=\$(echo "\$line" | awk -F'[, ]+' '{for(i=1;i<=NF;i++) if(\$i=="max") print \$(i+1)}')
    avg=\$(echo "\$line" | awk -F'[, ]+' '{for(i=1;i<=NF;i++) if(\$i=="avg") print \$(i+1)}')
    n50=\$(echo "\$line" | awk -F'[, ]+' '{for(i=1;i<=NF;i++) if(\$i=="N50") print \$(i+1)}')

    echo -e "number_of_contigs: \$numContigs\nshortest_contig: \$min\nlongest_contig: \$max\navg_contig_length: \$avg\nn50: \$n50" > megahit_summary.txt
    """
}

process unassembled_reads {
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "unassembled_reads*"
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "reads_mapped_to_contigs.*"

    input:
        val pair_id
        tuple path (megahit_contigs), path(assembly_summary)
        tuple path(fullyQc_1), path(fullyQc_2), path(fullyqc_summary)
        val output

    output:
        tuple path("unassembled_reads_fwd.fq.gz"),
            path("unassembled_reads_rev.fq.gz"),
            path("reads_mapped_to_contigs.cov_stats"),
            path("unassembled_summary.txt"),
            path("final.contigs.fq")
    
    
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

    awk '{s++} END {printf "total_unassembled_reads: %.0f\\n", (s/4)*2}' unassembled_reads_fwd.fq > unassembled_summary.txt
    gzip unassembled_reads_*.fq
    """
}

process megahit_fail {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "unassembled_reads*"
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy', pattern: "reads_mapped_to_contigs.cov_stats"

    input:
        val pair_id
        tuple path (megahit_contigs), path(assembly_summary)
        tuple path(fullyQc_1), path(fullyQc_2), path(fullyqc_summary)
        val cov_stats
        val output

    output:
        tuple path("unassembled_reads_fwd.fq.gz"),
            path("unassembled_reads_rev.fq.gz"),
            path("reads_mapped_to_contigs.cov_stats"),
            path("unassembled_summary.txt"),
            path("final.contigs.fq")

    script:
    """
    cp ${fullyQc_1} unassembled_reads_fwd.fq.gz
    cp ${fullyQc_2} unassembled_reads_rev.fq.gz
    cp ${cov_stats} reads_mapped_to_contigs.cov_stats
    cp ${megahit_contigs} final.contigs.fq
    awk '{s++} END {printf "total_unassembled_reads: %.0f\\n", (s/4)*2}' unassembled_reads_fwd.fq.gz > unassembled_summary.txt
    """
}

process alignment_prep {
    cache 'lenient'
    publishDir "${output}/${pair_id}/preprocessing", mode: 'copy'
    
    input:
        val pair_id
        path megahit_contigs
        tuple path(unassembled_reads_fwd), 
            path(unassembled_reads_rev),
            path(cov_stats),
            path(fqc_txt),
            path(final_contigs_fq)
        val output

    output:
        tuple path(megahit_contigs), 
            path("combined_sr_file.fq"),
            path("combined_sr_file.fa"),
            path("combined_reads_contigs_file.fq"),
            path("combined_reads_contigs_file.fa"),
            path("unassembled_reads_longer_fwd.fq"),
            path("unassembled_reads_longer_rev.fq"),
            path(cov_stats),
            path("alignment_prep_summary.txt")

    script:
    """
    ${params.python} ${baseDir}/scripts/separate_reads_by_size.py ${unassembled_reads_fwd} unassembled_reads_longer_fwd.fq unassembled_reads_shorter_fwd.fq
    ${params.python} ${baseDir}/scripts/separate_reads_by_size.py ${unassembled_reads_rev} unassembled_reads_longer_rev.fq unassembled_reads_shorter_rev.fq 

    seqtk seq -F '#' ${megahit_contigs} > megahit_contigs.fq
    seqtk mergepe unassembled_reads_shorter_fwd.fq unassembled_reads_shorter_rev.fq > combined_sr_file.fq
    seqtk seq -a combined_sr_file.fq > combined_sr_file.fa
    seqtk mergepe ${unassembled_reads_fwd} ${unassembled_reads_rev} > merged_reads.fq
    cat megahit_contigs.fq merged_reads.fq > combined_reads_contigs_file.fq
    seqtk seq -a combined_reads_contigs_file.fq > combined_reads_contigs_file.fa

    awk '{s++} END {printf "total_unassembled_shorter_reads: %.0f\\n", (s/4)}' combined_sr_file.fq > tmp_short.txt
    awk '{s++} END {printf "total_unassembled_longer_reads: %.0f\\n", (s/4)*2}' unassembled_reads_longer_fwd.fq > tmp_long.txt
    awk 'NR > 1 {sum += \$((NF-5))} END {printf "total_assembled_reads: %d\\n", sum}' ${cov_stats} > tmp_assembled.txt

    cat tmp_short.txt tmp_long.txt tmp_assembled.txt > alignment_prep_summary.txt
    """
}

process merged_assembly_summaries {
    input:
        path(megahit_summary)
        path(unassembled_reads_summary)
        path(alignment_prep_summary)

    output:
        path "assembly_summary.txt"

    script:
    """
    cat ${megahit_summary} ${unassembled_reads_summary} ${alignment_prep_summary} > assembly_summary.txt
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
        
        valid_files = megahit.out.filter { contigs, _summary -> contigs.size() > 0 }
        invalid_files = megahit.out.filter { contigs, _summary -> contigs.size() == 0 }
        
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

        aln_prep_data = alignment_prep(
            pair_id,
            megahit.out.map { contigs, _summary -> contigs },
            final_reads,
            output
        )
        
        merged_assembly_summaries(
            megahit.out.map { _contigs, summary -> summary },
            final_reads.map { _fwd, _rev, _cov_stats, summary, _final_fq -> summary },
            aln_prep_data.map { _contigs, _combined_sr_fq, _combined_sr_fa, _combined_lr_contigs_fq, _combined_lr_contigs_fa, _unassembled_longer_fwd, _unassembled_longer_rev, _cov_stats, summary -> summary }
        )
    emit:
        aln_prep_data = aln_prep_data
        assembly_summary = merged_assembly_summaries.out


}