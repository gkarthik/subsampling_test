process COVERAGE_STATS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)
    path reference

    output:
    path "${meta.id}.stats.txt", emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${projectDir}/bin/calculate_coverage_stats.py \\
        --bam ${bam} \\
        --reference ${reference} \\
        --sample_id ${meta.id} \\
        --sra_accession ${meta.sra_accession} \\
        --subsample ${meta.subsample ?: 'all'} \\
        --total_reads ${meta.total_reads ?: 0} \\
        --output ${meta.id}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        samtools: \$(samtools --version 2>&1 | sed 's/samtools //g')
    END_VERSIONS
    """

    stub:
    """
    echo "${meta.id},${meta.sra_accession},${meta.subsample ?: 'all'},${meta.total_reads ?: 0},0,0,0,0,0,0" > ${meta.id}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        samtools: \$(samtools --version 2>&1 | sed 's/samtools //g')
    END_VERSIONS
    """
} 