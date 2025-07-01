// Original Module: https://github.com/andersen-lab/flusra/blob/main/modules/nf-core/sratools/fasterqdump/main.nf
process SRATOOLS_FASTERQDUMP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path('*.fastq'), emit: reads
    path "*.fastq.gz"
    path "versions.yml", emit: versions

    script:
    """
    fasterq-dump \\
        --threads ${task.cpus} \\
        ${sra}

    pigz \\
        -p ${task.cpus} \\
        -k \\
        *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}