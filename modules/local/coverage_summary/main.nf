process COVERAGE_SUMMARY {
    tag "coverage_summary"
    label 'process_low'

    input:
    path stats_files

    output:
    path "*_coverage_stats.csv", emit: summaries
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${projectDir}/bin/coverage_summary.py \\
        --stats_files ${stats_files} \\
        --output_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo "subsample,total_reads,mapped_reads,mapping_rate,breadth_coverage,mean_depth" > test_coverage_stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
} 