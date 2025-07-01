process COVERAGE_ANALYSIS {
    tag "coverage_analysis"
    label 'process_medium'

    input:
    path variant_files

    output:
    path "output/*", emit: charts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${projectDir}/bin/coverage_analysis.py .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output
    touch output/dummy_chart.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
} 