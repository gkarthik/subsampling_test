process COVERAGE_ANALYSIS {
    tag "coverage_analysis"
    label 'process_medium'

    input:
    path variant_files

    output:
    path "coverage_analysis_charts", emit: charts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p coverage_analysis_charts
    python ${projectDir}/bin/coverage_analysis.py . 
    
    # Move output to specific directory
    if [ -d "output" ]; then
        cp -r output/* coverage_analysis_charts/ 2>/dev/null || true
        rm -rf output 2>/dev/null || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p coverage_analysis_charts
    touch coverage_analysis_charts/dummy_chart.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
} 