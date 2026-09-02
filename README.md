# bam-depth

A modular workflow for generating depth and coverage statistics for alignment files pre and post filtering

```bash
nextflow run j23414/bam-depth \
  --samplesheet [path/bam_samplesheet.csv] \
  --plot_depth_script [bin/plot_depth.R] \
  -profile stjude
```