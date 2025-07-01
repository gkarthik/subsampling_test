// Original Module: https://github.com/andersen-lab/flusra/blob/main/modules/nf-core/bwa/mem/main.nf
process BWA_MEM {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    bwa index -p "reference" ${reference}

    bwa mem \\
        -t ${task.cpus} \\
        "reference" \\
        ${reads} \\
        | samtools view \\
            -F 4 -b \\
            | samtools sort \\
                -o ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
