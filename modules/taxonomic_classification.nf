process minimap2_contigs {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "minimap2_contig_out*"
    
    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats)
        path mm2_index
        val output
    
    output:
        tuple path ("minimap2_contig_out.paf"),
        path ("minimap2_contig_out_frompaf.m8")
    
    script:
    """
    minimap2 -c --split-prefix mm2_contigs ${mm2_index} ${megahit_contigs} > minimap2_contig_out.paf
    python ${baseDir}/scripts/paf2blast6.py minimap2_contig_out.paf
    """
}

process minimap2_reads {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "minimap2_reads_out*"

    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats)
        path mm2_index
        val output
    
    output:
        tuple path ("minimap2_reads_out.paf"),
        path ("minimap2_reads_out_frompaf.m8")

    script:
    """
    minimap2 -c --split-prefix mm2_reads ${mm2_index} ${unassembled_reads_longer_fwd} ${unassembled_reads_longer_rev} > minimap2_reads_out.paf
    python ${baseDir}/scripts/paf2blast6.py minimap2_reads_out.paf
    """
}

process k2_pluspf {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "combined_rc.kraken.txt"

    input:
        val pair_id
        path pluspf_db
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats)
        val output
    
    output:
        path "combined_rc.kraken.txt"

    script:
    """
    cat ${combined_lr_contigs_fa} ${combined_sr_fa} > combined_rc.fa
    kraken2 --db ${pluspf_db} --use-names --threads 8 combined_rc.fa > combined_rc.kraken.txt
    """
}

process diamond {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "nr_alignments_file.tsv"

    input:
        val pair_id
        path diamond_db
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats)
        val output

    output:
        path "nr_alignments_file.tsv"
    
    script:
    """
    diamond blastx --db ${diamond_db} \
        --query ${combined_lr_contigs_fq} --mid-sensitive --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames\
        --masking 0 -c 1 -b 6 \
        --out nr_alignments_file.tsv
    """

}

process blast{
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "nt_alignments_sr_blast.tsv"

    input:
        val pair_id
        path blast_db
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev)
        val output
    
    output:
        path "nt_alignments_sr_blast.tsv"

    

    script:
    """
    idx_base=\$(readlink -f "${blast_db}")
    echo "Blast Database: \${idx_base}"

    blastn -task megablast \\
        -query "${combined_sr_fa}" \\
        -db "\${idx_base}" \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames" \\
        -out nt_alignments_sr_blast.tsv
    """
}

workflow taxonomic_classification {
    take:
        pair_id
        preprocessing_data
        mm2_index
        pluspf_db
        diamond_db
        blast_db
        output

    main:
        // Call Minimap2 on the contigs only.

        mm2_contigs = minimap2_contigs(
            pair_id,
            preprocessing_data,
            mm2_index,
            output
        )
        
        // Call Minimap2 on the longer reads.
        mm2_reads = minimap2_reads(
            pair_id,
            preprocessing_data,
            mm2_index,
            output
        )
        
        // Call Kraken2 with the PlusPF database
        k2_pluspf_data = k2_pluspf(
            pair_id,
            pluspf_db,
            preprocessing_data,
            output
        )

        // Call Diamond on the combined contigs and longer reads.
        diamond_data = diamond(
            pair_id,
            diamond_db,
            preprocessing_data,
            output
        )

        // Call Blast on the shorter reads.
        blast_data = blast(
            pair_id,
            blast_db,
            preprocessing_data,
            output
        )

    emit:
        mm2_contigs.concat(
            mm2_reads,
            k2_pluspf_data,
            diamond_data,
            blast_data
        ).view()
}