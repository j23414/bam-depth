include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH } from './modules/nf-core/samtools/depth/main'

process SAMTOOLS_COVERAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

  input: tuple val(meta), path(bam)
  output:
  tuple val(meta), path("${bam.baseName}_coverage.tsv"), emit: coverage
  tuple val(meta), path("${bam.baseName}_tophit.tsv"), emit: tophit

  script:
  def segments = params.segments
  """
  samtools coverage ${bam} > ${bam.baseName}_coverage.tsv
  cat ${bam.baseName}_coverage.tsv | sort -k6,6nr > ${bam.baseName}_sorted.tsv

  # 1=reference; 4=numreads; 6=coverage; 7=meandepth
  for segment in ${segments}; do
    grep "|\$segment" ${bam.baseName}_sorted.tsv \
        | awk -F'\t' -v seg="\$segment" 'OFS="\t" {{print "${bam.baseName}", seg, \$1, \$4, \$6, \$7}}' \
        | head -n1 \
        >> ${bam.baseName}_tophit.tsv
  done
  """
}

process PLOT_DEPTH {
  tag "${meta.id}"
  label 'process_low'

  input: tuple val(meta), path(depth), path(tophit)
  output: tuple val(meta), path("${depth.baseName}_depth.png")
  script:
  def plot_depth_script=file(params.plot_depth_script)
  """
  Rscript ${plot_depth_script} ${depth} ${tophit} ${depth.baseName}_depth.png
  """
}

workflow {
    main:
    // ============ Load Files ==================== //
    // BAM: Should be run through bam-filtered already
    // REF: Reference to call variants against
    // ============================================ //
    if (params.samplesheet) {
      bam_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], file(row.bam))
          }
    } else if (params.bam) {
        bam_ch = channel.fromPath(params.bam, checkIfExists: true)
          | map { bamfile -> tuple([id: bamfile.baseName], bamfile) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --bam 'data/*.bam'"
    }

    // Calculate a single coverage value
    bam_ch
    | SAMTOOLS_COVERAGE

    // Get read depth across all positions
    bam_ch
    | SAMTOOLS_INDEX

    SAMTOOLS_DEPTH(
      bam_ch | join(SAMTOOLS_INDEX.out.index),
      [[],[]] // All positions, so pass an empty positional value, this might get long if there are many ~200 references
    )

    if (params.plot_depth_script){
      SAMTOOLS_DEPTH.out.tsv
      | join(SAMTOOLS_COVERAGE.out.tophit)
      | PLOT_DEPTH
    }
}
