process COVERAGE_ANALYSIS {
    tag "coverage_analysis"
    label 'process_medium'

    input:
    path variant_files

    output:
    path "*/*.png", emit: charts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${projectDir}/bin/coverage_analysis.py .

    # Move sample result folders out of the intermediate "output" directory so that
    # charts reside directly under <sample_id>/.
    if [ -d "output" ]; then
        mv output/* ./ 2>/dev/null || true
        rmdir output 2>/dev/null || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dummy_sample
    touch dummy_sample/dummy_chart.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
} 