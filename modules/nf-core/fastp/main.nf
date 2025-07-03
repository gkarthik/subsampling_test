// Original Module: https://github.com/andersen-lab/flusra/blob/main/modules/nf-core/fastp/main.nf
process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads), val(trim_front_read_1), val(trim_front_read_2), val(trim_tail_read_1), val(trim_tail_read_2)

    output:
    tuple val(meta), path('trimmed_*.fastq'), emit: trimmed_reads
    tuple val(meta), path("${meta.id}_fastp.json"), emit: json_reports
    path "versions.yml", emit: versions

    script:
    
    // Filter out any single files if paired files exist
    def filtered_reads = reads.sort()
    if (filtered_reads.any { it.toString().contains('_1.fastq') } && filtered_reads.any { it.toString().contains('_2.fastq') }) {
        filtered_reads = filtered_reads.findAll { it.toString().contains('_1.fastq') || it.toString().contains('_2.fastq') }
    }
    
    def args = [
        trim_front_read_1 ? "--trim_front1 ${trim_front_read_1}" : null,
        filtered_reads.size() > 1 && trim_front_read_2 ? "--trim_front2 ${trim_front_read_2}" : null,
        trim_tail_read_1 ? "--trim_tail1 ${trim_tail_read_1}" : null,
        filtered_reads.size() > 1 && trim_tail_read_2 ? "--trim_tail2 ${trim_tail_read_2}" : null,
        "-q ${params.fastp_qualified_quality_phred ?: 15}",
        "-u ${params.fastp_unqualified_percent_limit ?: 40}",
        "-n ${params.fastp_n_base_limit ?: 5}",
        params.fastp_average_qual ? "-e ${params.fastp_average_qual}" : null
    ].grep()
    
    // Determine whether the input reads are single-end or paired-end.
    
    if (filtered_reads.size() == 2) {
        """
        fastp \\
            -i ${filtered_reads[0]} \\
            -I ${filtered_reads[1]} \\
            --thread ${task.cpus} \\
            ${args.join(' ')} \\
            -o trimmed_${meta.id}_1.fastq \\
            -O trimmed_${meta.id}_2.fastq \\
            --json ${meta.id}_fastp.json \\
            --html ${meta.id}_fastp.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
    else if (filtered_reads.size() == 1) {
        """
        fastp \\
            -i ${filtered_reads[0]} \\
            --thread ${task.cpus} \\
            ${args.join(' ')} \\
            -o trimmed_${meta.id}_1.fastq \\
            --json ${meta.id}_fastp.json \\
            --html ${meta.id}_fastp.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
    else {
        error("Invalid number of input reads: ${filtered_reads}")
    }

    stub:
    """
    touch ${meta.id}_1.fastq
    touch ${meta.id}_2.fastq
    touch trimmed_${meta.id}_1.fastq
    touch trimmed_${meta.id}_2.fastq
    touch ${meta.id}_fastp.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}