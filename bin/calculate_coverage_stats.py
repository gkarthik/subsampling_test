#!/usr/bin/env python3

import subprocess
import argparse
import json
import os
import sys

def calculate_coverage_stats(bam_file, reference_file, total_reads):
    print(f"Calculating coverage for {bam_file}...")
    
    # Get basic stats with samtools
    flagstat_cmd = f"samtools flagstat {bam_file}"
    flagstat_result = subprocess.run(flagstat_cmd, shell=True, capture_output=True, text=True)
    
    mapped_reads = 0
    if flagstat_result.returncode == 0:
        for line in flagstat_result.stdout.split('\n'):
            if 'mapped (' in line:
                mapped_reads = int(line.split()[0])
    else:
        print(f"Error with samtools flagstat: {flagstat_result.stderr}")
        return None
    
    # Calculate coverage with samtools depth
    cmd = f"samtools depth {bam_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error with samtools depth: {result.stderr}")
        return None
    
    # Parse coverage
    total_bases = 0
    covered_bases = 0
    total_depth = 0
    
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) == 3:
                depth = int(parts[2])
                total_bases += 1
                total_depth += depth
                if depth > 0:
                    covered_bases += 1
    
    breadth_coverage = (covered_bases / total_bases * 100) if total_bases > 0 else 0
    mean_depth = (total_depth / total_bases) if total_bases > 0 else 0
    mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0
    
    return {
        'total_reads': total_reads,        # Total read count from SRA
        'mapped_reads': mapped_reads,
        'mapping_rate': round(mapping_rate, 2),
        'breadth_coverage': round(breadth_coverage, 2),
        'mean_depth': round(mean_depth, 2),
        'covered_bases': covered_bases,
        'total_bases': total_bases
    }

def main():
    parser = argparse.ArgumentParser(description='Calculate coverage statistics from BAM file')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--reference', required=True, help='Reference FASTA file')
    parser.add_argument('--sample_id', required=True, help='Sample identifier')
    parser.add_argument('--sra_accession', required=True, help='SRA accession')
    parser.add_argument('--subsample', required=True, help='Subsample level')
    parser.add_argument('--total_reads', type=int, required=True, help='Total read count from SRA')
    parser.add_argument('--output', required=True, help='Output stats file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam):
        print(f"Error: BAM file {args.bam} does not exist")
        sys.exit(1)
    
    if not os.path.exists(args.reference):
        print(f"Error: Reference file {args.reference} does not exist")
        sys.exit(1)
    
    print(f"Processing sample: {args.sample_id}")
    print(f"SRA: {args.sra_accession}, Subsample: {args.subsample}")
    print(f"Total SRA reads: {args.total_reads:,}")
    print(f"BAM file: {args.bam}")
    print(f"Reference: {args.reference}")
    
    # Calculate coverage statistics
    stats = calculate_coverage_stats(args.bam, args.reference, args.total_reads)
    
    if stats:
        print(f"\n=== COVERAGE STATISTICS for {args.sample_id} ===")
        print(f"Total SRA reads: {args.total_reads:,}")
        print(f"Mapped reads: {stats['mapped_reads']:,}")
        print(f"Mapping rate: {stats['mapping_rate']:.2f}%")
        print(f"Breadth coverage: {stats['breadth_coverage']:.2f}%")
        print(f"Mean depth: {stats['mean_depth']:.2f}x")
        print(f"Covered bases: {stats['covered_bases']:,}")
        print(f"Total bases: {stats['total_bases']:,}")
        
        # Write stats to output file (CSV format) - now includes original read count
        with open(args.output, 'w') as f:
            f.write(f"{args.sample_id},{args.sra_accession},{args.subsample},{args.total_reads},{stats['mapped_reads']},{stats['mapping_rate']},{stats['breadth_coverage']},{stats['mean_depth']},{stats['covered_bases']},{stats['total_bases']}\n")
        
        print(f"Results saved to: {args.output}")
        
    else:
        print(f"Failed to calculate coverage statistics for {args.sample_id}")
        sys.exit(1)

if __name__ == "__main__":
    main() 