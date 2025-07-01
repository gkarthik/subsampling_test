// Original Module: https://github.com/andersen-lab/flusra/blob/main/modules/nf-core/sratools/prefetch/main.nf
process SRATOOLS_PREFETCH {
    tag "${id}"
    label 'process_low'

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path(id), emit: sra
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args = '5 1 100'
    // <num retries> <base delay in seconds> <max delay in seconds>
    template('retry_with_backoff.sh')
}