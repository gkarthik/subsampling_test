# Subsampling Test Pipeline

Nextflow pipeline for testing the effects of read subsampling on viral variant calling accuracy.

## Usage

### Docker
```bash
docker build -t subsampling-test:latest .
nextflow run test_subsampling.nf -profile docker
```

### Conda
```bash
mamba env create -f env.yml
conda activate subsampling-test
nextflow run test_subsampling.nf
```

## Pipeline

1. Downloads SRA data
2. Quality control with fastp
3. Subsamples reads at multiple levels (30K, 100K, 1M, 5M, 10M, all)
4. Aligns with BWA-MEM
5. Calls variants with ivar
6. Generates coverage analysis charts

## Output

Results in `subsampling_test_outputs/coverage_analysis/`:
- Scatter plots of variant frequency correlations
- Bar charts of variant counts
- Histograms of variant distributions
- R² correlation statistics 