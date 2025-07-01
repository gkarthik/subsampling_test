process CHECK_PAIRING {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (reads.size() == 2) {
        """
        # Check if files are properly paired by examining read headers
        echo "Checking if files are properly paired..."
        echo "File 1: ${reads[0]} (\$(stat -c%s ${reads[0]} 2>/dev/null || echo 0) bytes)"
        echo "File 2: ${reads[1]} (\$(stat -c%s ${reads[1]} 2>/dev/null || echo 0) bytes)"
        
        # Check if second file is empty or very small (indicating fake pairing)
        file2_size=\$(stat -c%s ${reads[1]} 2>/dev/null || echo 0)
        
        # Determine if files are compressed
        if [[ "${reads[0]}" == *.gz ]]; then
            read_cmd="zcat"
        else
            read_cmd="cat"
        fi
        
        if [ "\$file2_size" -le 50 ]; then
            echo "Second file is empty or tiny (\$file2_size bytes) - treating as single-end"
            \$read_cmd ${reads[0]} | gzip > ${prefix}_1.fastq.gz
        else
            # Both files have content - check if they're properly paired
            # Extract first read header from each file (safely)
            header1=\$(\$read_cmd ${reads[0]} | head -1 2>/dev/null || echo "")
            header2=\$(\$read_cmd ${reads[1]} | head -1 2>/dev/null || echo "")
            
            echo "Header 1: \$header1"
            echo "Header 2: \$header2"
            
            # Extract read IDs (remove /1, /2 suffixes)
            read_id1=\$(echo "\$header1" | cut -d' ' -f1 | sed 's|/1\$||' | sed 's|/2\$||')
            read_id2=\$(echo "\$header2" | cut -d' ' -f1 | sed 's|/1\$||' | sed 's|/2\$||')
            
            # More flexible pairing check - check if read IDs match (with or without /1, /2 suffixes)
            if [ "\$read_id1" = "\$read_id2" ] && [ -n "\$read_id1" ]; then
                # Additional check: see if headers have different suffixes or read numbers
                if [[ "\$header1" != "\$header2" ]]; then
                    echo "Files are properly paired - keeping as paired-end"
                    \$read_cmd ${reads[0]} | gzip > ${prefix}_1.fastq.gz
                    \$read_cmd ${reads[1]} | gzip > ${prefix}_2.fastq.gz
                else
                    echo "Files have identical headers - combining for single-end alignment"
                    \$read_cmd ${reads[0]} ${reads[1]} | gzip > ${prefix}_1.fastq.gz
                fi
            else
                echo "Files are not properly paired - combining for single-end alignment"
                # Combine files for single-end processing
                \$read_cmd ${reads[0]} ${reads[1]} | gzip > ${prefix}_1.fastq.gz
            fi
        fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bash: \$(bash --version | head -1 | cut -d' ' -f4)
        END_VERSIONS
        """
    } else {
        """
        echo "Single file detected - keeping as single-end"
        
        # Determine if file is compressed
        if [[ "${reads[0]}" == *.gz ]]; then
            read_cmd="zcat"
        else
            read_cmd="cat"
        fi
        
        \$read_cmd ${reads[0]} | gzip > ${prefix}_1.fastq.gz
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bash: \$(bash --version | head -1 | cut -d' ' -f4)
        END_VERSIONS
        """
    }
} 