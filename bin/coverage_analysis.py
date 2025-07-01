import pandas as pd
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import sys
import argparse

plt.switch_backend('Agg')

def find_variant_files(variants_dir):
    """Find all variant TSV files in the specified directory and subdirectories"""
    variant_files = []
    
    variant_files.extend(glob.glob(os.path.join(variants_dir, "*_genome_variants.tsv")))
    
    for subdir in glob.glob(os.path.join(variants_dir, "*/")):
        variant_files.extend(glob.glob(os.path.join(subdir, "*_genome_variants.tsv")))
    
    return variant_files

def extract_sample_info(filename):
    """Extract sample ID and subsample level from filename"""
    basename = os.path.basename(filename)
    parts = basename.replace('_genome_variants.tsv', '').split('_')
    
    sample_parts = []
    subsample_level = None
    
    for i, part in enumerate(parts):
        if part.isdigit() or part == 'all':
            subsample_level = part
            sample_parts = parts[:i]
            break
    
    sample_id = '_'.join(sample_parts) if sample_parts else 'unknown'
    
    return sample_id, subsample_level

def main():
    parser = argparse.ArgumentParser(description='Generate variant analysis charts from a specified directory')
    parser.add_argument('variants_dir', help='Directory containing variant TSV files')
    args = parser.parse_args()
    
    if not os.path.exists(args.variants_dir):
        print(f"Error: Directory '{args.variants_dir}' does not exist.")
        return
    
    variant_files = find_variant_files(args.variants_dir)
    print(f"Found {len(variant_files)} variant files in {args.variants_dir}")
    
    if not variant_files:
        print("No variant files found. Make sure the directory contains *_genome_variants.tsv files.")
        return
    
    # Group files by sample
    samples = {}
    for file in variant_files:
        sample_id, subsample_level = extract_sample_info(file)
        if sample_id not in samples:
            samples[sample_id] = {}
        samples[sample_id][subsample_level] = file
        print(f"Found: {sample_id} - {subsample_level} -> {file}")
    
    # Process each sample
    for sample_id, files_dict in samples.items():
        print(f"\n=== Processing {sample_id} ===")
        
        output_dir = os.path.join("output", sample_id)
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory: {output_dir}")
        
        if 'all' not in files_dict:
            print(f"No 'all' reference file found for {sample_id}, skipping")
            continue
        
        all_file = files_dict['all']
        try:
            df_all = pd.read_csv(all_file, sep='\t')
            df_all_filtered = df_all[df_all['ALT_FREQ'] > 0.1]
            print(f"Reference data: {len(df_all)} total variants, {len(df_all_filtered)} with ALT_FREQ > 0.1")
        except Exception as e:
            print(f"Error reading {all_file}: {e}")
            continue
        
        if len(df_all_filtered) == 0:
            print(f"No variants in reference for {sample_id}")
            continue
        
        # Process each subsampling level
        subsample_levels = [k for k in files_dict.keys() if k != 'all' and k.isdigit()]
        results = []
        
        for level in sorted(subsample_levels, key=int):
            file = files_dict[level]
            print(f"\nProcessing {level} reads...")
            
            try:
                # Read subsampled data
                df_sub = pd.read_csv(file, sep='\t')
                df_sub_filtered = df_sub[df_sub['ALT_FREQ'] > 0.1]
                
                print(f"  Subsample data: {len(df_sub)} total variants, {len(df_sub_filtered)} with ALT_FREQ > 0.1")
                
                if len(df_sub_filtered) > 0:
                    # Merge with reference data
                    df_merged = pd.merge(
                        df_sub_filtered, 
                        df_all_filtered, 
                        on=['REGION', 'POS', 'REF', 'ALT'], 
                        suffixes=(f'_{level}', '_all')
                    )
                    
                    print(f"  Merged: {len(df_merged)} common variants")
                    
                    if len(df_merged) > 1:
                        # Calculate R² score
                        r2 = r2_score(df_merged[f'ALT_FREQ_{level}'], df_merged['ALT_FREQ_all'])
                        print(f"  R² score: {r2:.4f}")
                        
                        # Create scatter plot
                        plt.figure(figsize=(8, 6))
                        plt.scatter(df_merged[f'ALT_FREQ_{level}'], df_merged['ALT_FREQ_all'], alpha=0.6)
                        plt.xlabel(f'ALT_FREQ_{level}')
                        plt.ylabel('ALT_FREQ_all')
                        plt.title(f'Variant Frequency Correlation\nSample: {sample_id}, Subsample: {level}, R² = {r2:.4f}')
                        plt.plot([0,1], [0,1], color='red', linestyle='--', alpha=0.8)
                        plt.xlim(0, 1)
                        plt.ylim(0, 1)
                        plt.grid(True, alpha=0.3)
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_dir, f'scatter_{level}.png'), dpi=150, bbox_inches='tight')
                        plt.close()
                        
                        # Create variant count comparison bar chart
                        plt.figure(figsize=(8, 6))
                        counts = [len(df_sub_filtered), len(df_all_filtered)]
                        labels = [f'{level} reads', 'All reads']
                        plt.bar(labels, counts, color=['lightblue', 'lightcoral'])
                        plt.ylabel('Number of Variants (ALT_FREQ > 0.1)')
                        plt.title(f'Variant Count Comparison\nSample: {sample_id}, Subsample: {level}')
                        for i, count in enumerate(counts):
                            plt.text(i, count + max(counts)*0.01, str(count), ha='center', va='bottom')
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_dir, f'bar_{level}.png'), dpi=150, bbox_inches='tight')
                        plt.close()
                        
                        # Create histogram of variant frequencies
                        plt.figure(figsize=(12, 5))
                        
                        plt.subplot(1, 2, 1)
                        plt.hist(df_sub_filtered['ALT_FREQ'], bins=30, alpha=0.7, color='lightblue', edgecolor='black')
                        plt.xlabel('Alternative Allele Frequency')
                        plt.ylabel('Count')
                        plt.title(f'Variant Frequencies\n{level} reads')
                        plt.grid(True, alpha=0.3)
                        
                        plt.subplot(1, 2, 2)
                        plt.hist(df_all_filtered['ALT_FREQ'], bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
                        plt.xlabel('Alternative Allele Frequency')
                        plt.ylabel('Count')
                        plt.title('Variant Frequencies\nAll reads')
                        plt.grid(True, alpha=0.3)
                        
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_dir, f'hist_{level}.png'), dpi=150, bbox_inches='tight')
                        plt.close()
                        
                        results.append({
                            'subsample_level': level,
                            'total_variants_subsample': len(df_sub),
                            'filtered_variants_subsample': len(df_sub_filtered),
                            'total_variants_all': len(df_all),
                            'filtered_variants_all': len(df_all_filtered),
                            'common_variants': len(df_merged),
                            'r2_score': r2
                        })
                    else:
                        print(f"  Not enough common variants for analysis")
                        
            except Exception as e:
                print(f"Error processing {file}: {e}")
        
        if results:
            print(f"\nSUMMARY FOR {sample_id}:")
            print(f"Reference variants (ALT_FREQ > 0.1): {len(df_all_filtered)}")
            for result in results:
                print(f"  {result['subsample_level']} reads: R² = {result['r2_score']:.4f}, "
                      f"Common variants = {result['common_variants']}")
            
            # Calculate mean R²
            r2_scores = [r['r2_score'] for r in results]
            if r2_scores:
                print(f"  Mean R² score: {np.mean(r2_scores):.4f}")
            
            print(f"Generated {len(results)} analysis sets for {sample_id} in {output_dir}")

if __name__ == "__main__":
    main() 