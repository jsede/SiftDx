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
            path(cov_stats),
            path(fqc_txt)
        val output
    
    output:
        tuple path ("minimap2_contig_out.paf"),
        path ("minimap2_contig_out_frompaf.m8")
    
    script:
    """
    minimap2 -c --split-prefix mm2_contigs ${params.mm2_index} ${megahit_contigs} > minimap2_contig_out.paf
    ${params.python} ${baseDir}/scripts/paf2blast6.py minimap2_contig_out.paf
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
            path(cov_stats),
            path(fqc_txt)
        val output
    
    output:
        tuple path ("minimap2_reads_out.paf"),
        path ("minimap2_reads_out_frompaf.m8")

    script:
    """
    minimap2 -c -x sr --split-prefix mm2_reads ${params.mm2_index} ${unassembled_reads_longer_fwd} ${unassembled_reads_longer_rev} > minimap2_reads_out.paf
    ${params.python} ${baseDir}/scripts/paf2blast6.py minimap2_reads_out.paf
    """
}

process k2_pluspf {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "combined_rc.kraken.txt"

    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats),
            path(fqc_txt)
        val output
    
    output:
        path "combined_rc.kraken.txt"

    script:
    """
    cat ${combined_lr_contigs_fa} ${combined_sr_fa} > combined_rc.fa
    kraken2 --db ${params.pluspf_db} --use-names --threads 8 combined_rc.fa > combined_rc.kraken.txt
    """
}

process diamond {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "nr_alignments_file.tsv"

    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats),
            path(fqc_txt)
        val output

    output:
        path "nr_alignments_file.tsv"
    
    script:
    """
    diamond blastx --db ${params.diamond_db} \
        --query ${combined_lr_contigs_fq} --mid-sensitive --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids\
        --masking 0 -c 1 -b 6 \
        --out nr_alignments_file.tsv
    """

}

process split_input {
    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats),
            path(fqc_txt)
        val output

    output:
        path "chunk_*"

    script:
    """
    ${params.python} ${baseDir}/scripts/split_fasta.py ${combined_sr_fa} 20000
    """
}

process split_blast {
    input:
        val pair_id
        path chunk
        val output

    output:
        path "${chunk.simpleName}.blast.tsv"

    script:
    """
    blastn -task megablast \
        -query ${chunk} \
        -db ${params.blast_db} \
        -max_target_seqs 10 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids" \
        -out ${chunk.simpleName}.blast.tsv
    """
}

process merge_blast {
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "nt_alignments_sr_blast.tsv"

    input:
        val pair_id
        path blast_chunks
        val output

    output:
        path "nt_alignments_sr_blast.tsv"

    script:
    """
    cat ${blast_chunks.join(' ')} > nt_alignments_sr_blast.tsv
    """
}
process blast{
    publishDir "${output}/${pair_id}/alignments", mode: 'copy', pattern: "nt_alignments_sr_blast.tsv"

    input:
        val pair_id
        tuple path(megahit_contigs),
            path(combined_sr_fq),
            path(combined_sr_fa),
            path(combined_lr_contigs_fq),
            path(combined_lr_contigs_fa),
            path(unassembled_reads_longer_fwd),
            path(unassembled_reads_longer_rev),
            path(cov_stats),
            path(fqc_txt)
        val output
    
    output:
        path "nt_alignments_sr_blast.tsv"

    script:
    """
    blastn -task megablast \\
        -query "${combined_sr_fa}" \\
        -db "${params.blast_db}" \\
        -max_target_seqs 10 \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids" \\
        -out nt_alignments_sr_blast.tsv
    """
}

workflow taxonomic_classification {
    take:
        pair_id
        preprocessing_data
        output

    main:
        // Call Minimap2 on the contigs only.

        mm2_contigs = minimap2_contigs(
            pair_id,
            preprocessing_data,
            output
        )
        
        // Call Minimap2 on the longer reads.
        mm2_reads = minimap2_reads(
            pair_id,
            preprocessing_data,
            output
        )
        
        // Call Kraken2 with the PlusPF database
        k2_pluspf_data = k2_pluspf(
            pair_id,
            preprocessing_data,
            output
        )

        // Call Diamond on the combined contigs and longer reads.
        diamond_data = diamond(
            pair_id,
            preprocessing_data,
            output
        )

        // Call Blast on the shorter reads.
        
        blast_file_ch = preprocessing_data.map { 
            _megahit_contigs, _combined_sr_fq, combined_sr_fa, _combined_lr_contigs_fq,
            _combined_lr_contigs_fa, _unassembled_fwd, _unassembled_rev, _cov_stats, _fqc_txt ->
                combined_sr_fa 
        }
        blast_file = blast_file_ch.toList().get(0)
        blast_size = file(blast_file).getSize()
        if (blast_size > 4_000_000) {
            blast_split = split_input(
                pair_id,
                preprocessing_data,
                output
            )

            blast_chunks = split_blast(
                pair_id,
                blast_split.out,
                output
            )

            blast_data = merge_blast(
                pair_id,
                blast_chunks.out,
                output
            )
        } else {
            blast_data = blast(
                pair_id,
                preprocessing_data,
                output
            )
        }
    emit:
        mm2_contigs
            .concat(mm2_reads)
            .concat(k2_pluspf_data)
            .concat(diamond_data)
            .concat(blast_data)
            .flatten().toList().view()
}