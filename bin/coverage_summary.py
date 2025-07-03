#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Create CSV summaries grouped by SRA accession')
    parser.add_argument('--stats_files', nargs='+', required=True, help='List of stats files')
    parser.add_argument('--output_dir', required=True, help='Output directory for CSV files')
    
    args = parser.parse_args()
    
    stats_files = args.stats_files
    print(f"Processing {len(stats_files)} stats files")
    
    # Read all stats files
    all_data = []
    
    for stats_file in stats_files:
        if not os.path.exists(stats_file):
            print(f"Warning: Stats file {stats_file} does not exist, skipping")
            continue
        
        try:
            with open(stats_file, 'r') as f:
                line = f.read().strip()
                parts = line.split(',')
                if len(parts) == 10:
                    sample_id, sra_accession, subsample, total_reads, mapped_reads, mapping_rate, breadth_coverage, mean_depth, covered_bases, total_bases = parts
                    all_data.append({
                        'sample_id': sample_id,
                        'sra_accession': sra_accession,
                        'subsample': subsample,
                        'total_reads': int(total_reads),
                        'mapped_reads': int(mapped_reads),
                        'mapping_rate': float(mapping_rate),
                        'breadth_coverage': float(breadth_coverage),
                        'mean_depth': float(mean_depth),
                        'covered_bases': int(covered_bases),
                        'total_bases': int(total_bases)
                    })
                    print(f"  Loaded: {sample_id} (SRA: {sra_accession}, Subsample: {subsample})")
                else:
                    print(f"Error: Invalid format in {stats_file}")
        except Exception as e:
            print(f"Error reading {stats_file}: {e}")
    
    if not all_data:
        print("No data loaded")
        return
    
    # Group by SRA accession
    sra_groups = {}
    for data in all_data:
        sra_acc = data['sra_accession']
        if sra_acc not in sra_groups:
            sra_groups[sra_acc] = []
        sra_groups[sra_acc].append(data)
    
    print(f"\nGrouped into {len(sra_groups)} SRA accessions")
    
    # Create CSV for each SRA accession
    for sra_acc, sra_data in sra_groups.items():
        print(f"\nProcessing SRA: {sra_acc} ({len(sra_data)} subsamples)")
        
        # Create DataFrame
        df = pd.DataFrame(sra_data)
        
        # Select and order columns for per-SRA CSV
        columns = [
            'subsample',
            'total_reads', 
            'mapped_reads',
            'mapping_rate',
            'breadth_coverage',
            'mean_depth',
            'covered_bases',
            'total_bases'
        ]
        
        df = df[columns]
        
        # Sort by subsample (convert to int for proper sorting, handle 'all')
        def sort_key(x):
            if x == 'all':
                return float('inf')  # 'all' goes last
            try:
                return int(x)
            except:
                return 0
        
        df['sort_key'] = df['subsample'].apply(sort_key)
        df = df.sort_values('sort_key').drop('sort_key', axis=1)
        
        # Save to CSV
        output_file = os.path.join(args.output_dir, f"{sra_acc}_coverage_stats.csv")
        df.to_csv(output_file, index=False)
        
        print(f"  Saved to: {output_file}")
        print(f"  Subsamples: {', '.join(df['subsample'].astype(str))}")
        
        # Print some summary stats for this SRA
        print(f"  Mean mapping rate: {df['mapping_rate'].mean():.2f}%")
        print(f"  Mean breadth coverage: {df['breadth_coverage'].mean():.2f}%")
        print(f"  Mean depth range: {df['mean_depth'].min():.1f}x - {df['mean_depth'].max():.1f}x")
    
    print(f"\n=== SUMMARY ===")
    print(f"Created {len(sra_groups)} CSV files for {len(all_data)} total samples")

if __name__ == "__main__":
    main() 