process SEQTK_SAMPLE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)
    val(num_reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seed = 42 // Fixed seed for reproducibility
    
    if (num_reads == 'all') {
        if (reads.size() == 2) {
            // Paired-end
            """
            echo "SEQTK_SAMPLE: Processing paired-end data (all reads)"
            echo "Input files: ${reads}"
            cp ${reads[0]} ${prefix}_subsampled_1.fastq.gz
            cp ${reads[1]} ${prefix}_subsampled_2.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                seqtk: \$(seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d' ' -f2)
            END_VERSIONS
            """
        } else {
            // Single-end
            """
            echo "SEQTK_SAMPLE: Processing single-end data (all reads)"
            echo "Input files: ${reads}"
            cp ${reads[0]} ${prefix}_subsampled_1.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                seqtk: \$(seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d' ' -f2)
            END_VERSIONS
            """
        }
    } else {
        if (reads.size() == 2) {
            // Paired-end
            """
            echo "SEQTK_SAMPLE: Processing paired-end data (${num_reads} reads)"
            echo "Input files: ${reads}"
            seqtk sample -s ${seed} ${reads[0]} ${num_reads} | gzip > ${prefix}_subsampled_1.fastq.gz
            seqtk sample -s ${seed} ${reads[1]} ${num_reads} | gzip > ${prefix}_subsampled_2.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                seqtk: \$(seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d' ' -f2)
            END_VERSIONS
            """
        } else {
            // Single-end
            """
            echo "SEQTK_SAMPLE: Processing single-end data (${num_reads} reads)"
            echo "Input files: ${reads}"
            seqtk sample -s ${seed} ${reads[0]} ${num_reads} | gzip > ${prefix}_subsampled_1.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                seqtk: \$(seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d' ' -f2)
            END_VERSIONS
            """
        }
    }
} 