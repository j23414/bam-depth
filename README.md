# bam-depth

A modular workflow for generating depth and coverage statistics for alignment files pre and post filtering

```bash
nextflow run j23414/bam-depth \
  --samplesheet [path/bam_samplesheet.csv] \
  --plot_depth_script [bin/plot_depth.R] \
  -profile stjude
```

<img src="https://raw.githubusercontent.com/j23414/bam-depth/refs/heads/main/diagram.png" width="350">

## Optional: post-process merging for coverage tables

```
# (option A) Just get table of coverage
cat filter-results/samtools_coverage/*tophit* > all-top-coverage.tsv
python bam-filter/bin/coverage_to_wide.py \
  --top-coverage all-top-coverage.tsv \
  --merged results-coverage.tsv
```