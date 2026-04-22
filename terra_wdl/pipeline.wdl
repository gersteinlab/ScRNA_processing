## pipeline.wdl
##
## Single-cell / single-nucleus RNA-seq processing pipeline for Terra (Cromwell/WDL 1.0).
## Mirrors the three-stage local pipeline:
##   Stage 1  – CellRanger count  (one task per sample, scattered)
##   Stage 1b – CITE-seq-Count    (one task per hashed sample, conditional scatter)
##   Stage 2  – CellBender        (one task per sample, optional GPU, conditional scatter)
##   Stage 3  – Pegasus analysis  (one task for the entire batch)
##
## One Terra workflow submission = one sample_set (batch).
## Use generate_sample_table.py locally to produce sample_table.tsv and
## sample_set_table.tsv before uploading to Terra.
##
## Pilot / bypass mode (skip_cellranger = true):
##   Provide precomputed_h5_urls (HTTPS or gs:// URLs to CellBender-filtered .h5 files).
##   Stages 1, 1b, and 2 are skipped entirely; each URL is downloaded by LocalizeH5
##   and the resulting files are passed directly to Pegasus (Stage 3).
##
## Docker images:
##   cellranger_docker  – majidfarhadloo/scrna_processing_cellranger:latest
##                        (built by .github/workflows/build-cellranger-image.yml)
##                        contains: cellranger 10.0.0, gsutil, Python 3
##   cellbender_docker  – must contain: cellbender, CUDA drivers, Python 3
##   pegasus_docker     – majidfarhadloo/scrna_processing_pegasus:latest
##                        (built by .github/workflows/build-pegasus-image.yml)
##                        contains all pipeline scripts at /opt/pipeline/

version 1.0

# ---------------------------------------------------------------------------
# WORKFLOW
# ---------------------------------------------------------------------------

workflow SCRNAseqPipeline {

  input {
    # -----------------------------------------------------------------------
    # Per-sample arrays (parallel; one element per sample in this batch)
    # All arrays must have the same length.
    # -----------------------------------------------------------------------

    ## Sample name matching the FASTQ file prefix (CellRanger --sample / --id)
    ## Also used as metadata labels by Pegasus regardless of skip_cellranger mode.
    Array[String] sample_tags

    ## Whether each sample has cell hashing (parallel to sample arrays).
    ## Required even when skip_cellranger=true so Pegasus knows hashing status.
    Array[Boolean] hashing_flags

    # -----------------------------------------------------------------------
    # Stage 1 inputs – only required when skip_cellranger = false
    # -----------------------------------------------------------------------

    ## GCS paths to per-sample FASTQ directories (folders, not individual files)
    Array[String] fastq_dirs = []

    ## 10x chemistry: "v2", "v3", or "auto"
    Array[String] chemistry_flags = []

    ## "true" / "false" – include intronic reads (snRNA-seq = true)
    Array[String] include_introns_flags = []

    ## Hashing FASTQ R1 / R2 paths (empty string "" for non-hashed samples)
    Array[String] hashing_fastq_r1 = []
    Array[String] hashing_fastq_r2 = []

    ## Cell-hashing tag CSV paths (empty string "" for non-hashed samples)
    Array[String] hashing_tag_files = []

    ## GCS path to CellRanger reference transcriptome directory
    String transcriptome_gcs_path = ""

    # -----------------------------------------------------------------------
    # Bypass / pilot mode
    # -----------------------------------------------------------------------

    ## Set true to skip Stages 1, 1b, 2 entirely.
    ## Requires precomputed_h5_urls to be provided (same length as sample_tags).
    Boolean skip_cellranger = false

    ## HTTPS or gs:// URLs to pre-processed (CellBender-filtered) .h5 files,
    ## one per sample. Only used when skip_cellranger = true.
    Array[String] precomputed_h5_urls = []

    # -----------------------------------------------------------------------
    # Shared reference inputs
    # -----------------------------------------------------------------------

    ## GCS path to mitochondrial gene list CSV (used by Pegasus)
    File mito_file

    # -----------------------------------------------------------------------
    # Batch-level metadata (used by Pegasus)
    # -----------------------------------------------------------------------

    ## Human-readable batch / sample-set identifier (e.g. MS1023HH_prefrontal_cortex)
    String batch_name

    ## Overall run prefix for output naming
    String run_prefix

    # -----------------------------------------------------------------------
    # Pipeline flags
    # -----------------------------------------------------------------------

    ## Run CellBender ambient RNA removal (requires GPU worker).
    ## Ignored when skip_cellranger = true.
    Boolean run_cellbender = true

    ## GCP billing project for requester-pays GCS buckets (e.g. "terra-4ea165c8").
    ## Leave empty ("") if FASTQ buckets are not requester-pays.
    String billing_project = ""

    # -----------------------------------------------------------------------
    # Docker image URIs (pass as workflow inputs so they can be overridden
    # without editing the WDL)
    # -----------------------------------------------------------------------

    String cellranger_docker  = "majidfarhadloo/scrna_processing_cellranger:latest"
    String cellbender_docker  = "us.gcr.io/broad-dsde-methods/cellbender:0.3.0"
    String pegasus_docker     = "majidfarhadloo/scrna_processing_pegasus:latest"
    String localize_h5_docker = "ubuntu:22.04"

    # -----------------------------------------------------------------------
    # Resource tuning (optional overrides)
    # -----------------------------------------------------------------------

    Int cellranger_cpu     = 16
    Int cellranger_mem_gb  = 64
    Int cellranger_disk_gb = 500

    Int cellbender_gpu     = 1
    Int cellbender_cpu     = 4
    Int cellbender_mem_gb  = 32
    Int cellbender_disk_gb = 100

    Int pegasus_cpu     = 16
    Int pegasus_mem_gb  = 128
    Int pegasus_disk_gb = 200

    Int localize_h5_disk_gb = 150

    # -----------------------------------------------------------------------
    # QC / algorithm parameters (Pegasus defaults match original pipeline)
    # -----------------------------------------------------------------------
    Int    qc_min_umis           = 500
    Int    qc_percent_mito       = 10
    Int    qc_min_genes          = 200
    Int    min_umi_check         = 1000   # CellRanger low-quality flag threshold
    Int    cellbender_epochs     = 150
    Float  cellbender_fpr        = 0.01
    Int    dd_bst_n_iters        = 25
    Float  dd_pred_pthresh       = 0.0000000000000001
    Float  dd_pred_voterthresh   = 0.3
    Int    hvg_n_top             = 5000
    Int    n_jobs_pg             = 10
  }

  # -------------------------------------------------------------------------
  # Bypass mode: localize pre-computed h5 files from HTTPS / GCS URLs
  # -------------------------------------------------------------------------

  scatter (url in precomputed_h5_urls) {
    call LocalizeH5 {
      input:
        url     = url,
        docker  = localize_h5_docker,
        disk_gb = localize_h5_disk_gb
    }
  }

  # -------------------------------------------------------------------------
  # Stage 1 / 1b / 2: CellRanger → CITE-seq-Count → CellBender
  # Skipped entirely when skip_cellranger = true.
  # -------------------------------------------------------------------------

  scatter (i in range(length(sample_tags))) {

    if (!skip_cellranger) {

      # ---------------------------------------------------------------------
      # Stage 1: CellRanger count
      # ---------------------------------------------------------------------
      call CellRangerCount {
        input:
          fastq_dir       = fastq_dirs[i],
          sample_tag      = sample_tags[i],
          transcriptome   = transcriptome_gcs_path,
          chemistry       = chemistry_flags[i],
          include_introns = include_introns_flags[i],
          numproc         = cellranger_cpu,
          localmem        = cellranger_mem_gb,
          min_umi_check   = min_umi_check,
          billing_project = billing_project,
          docker          = cellranger_docker,
          cpu             = cellranger_cpu,
          mem_gb          = cellranger_mem_gb,
          disk_gb         = cellranger_disk_gb
      }

      # ---------------------------------------------------------------------
      # Stage 1b: CITE-seq-Count – only for hashed samples
      # ---------------------------------------------------------------------
      if (hashing_flags[i]) {
        call CiteSeqCount {
          input:
            sample_tag      = sample_tags[i],
            fastq_r1        = hashing_fastq_r1[i],
            fastq_r2        = hashing_fastq_r2[i],
            tag_file        = hashing_tag_files[i],
            chemistry       = chemistry_flags[i],
            estimated_cells = select_first([CellRangerCount.estimated_cells]),
            numproc         = cellranger_cpu,
            billing_project = billing_project,
            docker          = pegasus_docker,
            cpu             = cellranger_cpu,
            mem_gb          = cellranger_mem_gb,
            disk_gb         = cellranger_disk_gb
        }
      }

      # ---------------------------------------------------------------------
      # Stage 2: CellBender – optional, one task per sample
      # ---------------------------------------------------------------------
      if (run_cellbender) {
        call CellBender {
          input:
            sample_tag      = sample_tags[i],
            input_h5        = select_first([CellRangerCount.raw_h5]),
            expected_cells  = select_first([CellRangerCount.estimated_cells]),
            epochs          = cellbender_epochs,
            fpr             = cellbender_fpr,
            docker          = cellbender_docker,
            gpu             = cellbender_gpu,
            cpu             = cellbender_cpu,
            mem_gb          = cellbender_mem_gb,
            disk_gb         = cellbender_disk_gb
        }
      }
    }
  }

  # -------------------------------------------------------------------------
  # Resolve the h5 list passed to Pegasus:
  #   skip_cellranger=true  → use files localized by LocalizeH5
  #   skip_cellranger=false, run_cellbender=true  → CellBender filtered h5
  #   skip_cellranger=false, run_cellbender=false → CellRanger filtered h5
  # -------------------------------------------------------------------------

  Array[File] h5_for_pegasus =
    if skip_cellranger
    then LocalizeH5.h5_file
    else if run_cellbender
      then select_all(CellBender.cb_filtered_h5)
      else select_all(CellRangerCount.filtered_h5)

  # -------------------------------------------------------------------------
  # Stage 3: Pegasus – one task for the whole batch
  # -------------------------------------------------------------------------

  call PegasusPipeline {
    input:
      batch_name          = batch_name,
      run_prefix          = run_prefix,
      sample_tags         = sample_tags,
      h5_files            = h5_for_pegasus,
      hto_count_dirs      = select_all(CiteSeqCount.hto_count_dir),
      hashing_flags       = hashing_flags,
      mito_file           = mito_file,
      n_jobs_pg           = n_jobs_pg,
      qc_min_umis         = qc_min_umis,
      qc_percent_mito     = qc_percent_mito,
      qc_min_genes        = qc_min_genes,
      dd_bst_n_iters      = dd_bst_n_iters,
      dd_pred_pthresh     = dd_pred_pthresh,
      dd_pred_voterthresh = dd_pred_voterthresh,
      hvg_n_top           = hvg_n_top,
      docker              = pegasus_docker,
      cpu                 = pegasus_cpu,
      mem_gb              = pegasus_mem_gb,
      disk_gb             = pegasus_disk_gb
  }

  # -------------------------------------------------------------------------
  # Workflow-level outputs
  # -------------------------------------------------------------------------

  output {
    # CellRanger per-sample outputs (empty when skip_cellranger=true)
    Array[File?]    cr_raw_h5            = CellRangerCount.raw_h5
    Array[File?]    cr_filtered_h5       = CellRangerCount.filtered_h5
    Array[File?]    cr_molecule_info     = CellRangerCount.molecule_info
    Array[File?]    cr_metrics_csv       = CellRangerCount.metrics_csv
    Array[Int?]     cr_estimated_cells   = CellRangerCount.estimated_cells
    Array[Boolean?] cr_low_quality_flags = CellRangerCount.is_low_quality

    # CellBender per-sample outputs (empty when skip_cellranger=true or run_cellbender=false)
    Array[File?]    cb_filtered_h5_list  = CellBender.cb_filtered_h5

    # Pegasus batch outputs
    File            pegasus_output_tar   = PegasusPipeline.output_tar
    File            pegasus_summary      = PegasusPipeline.summary_stats
  }
}


# ---------------------------------------------------------------------------
# TASK: LocalizeH5
# Downloads a single .h5 file from an HTTPS or gs:// URL.
# Used in bypass mode (skip_cellranger = true) to stage pre-processed h5 files
# before passing them to Pegasus.
# ---------------------------------------------------------------------------

task LocalizeH5 {
  input {
    String url

    String docker  = "ubuntu:22.04"
    Int    disk_gb = 150
  }

  command <<<
    set -euo pipefail
    apt-get update -qq && apt-get install -y -qq curl ca-certificates

    echo "Downloading: ~{url}"
    curl -L --retry 5 --retry-delay 10 --retry-max-time 300 \
      -o downloaded.h5 "~{url}"

    echo "Download complete. File size: $(du -sh downloaded.h5 | cut -f1)"
  >>>

  output {
    File h5_file = "downloaded.h5"
  }

  runtime {
    docker:      docker
    cpu:         2
    memory:      "4 GB"
    disks:       "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}


# ---------------------------------------------------------------------------
# TASK: CellRangerCount
# ---------------------------------------------------------------------------

task CellRangerCount {
  input {
    String  fastq_dir
    String  sample_tag
    String  transcriptome
    String  chemistry       = "auto"
    String  include_introns = "true"
    Int     numproc         = 16
    Int     localmem        = 64
    Int     min_umi_check   = 1000
    String  billing_project = ""

    String  docker
    Int     cpu     = 16
    Int     mem_gb  = 64
    Int     disk_gb = 500
  }

  command <<<
    set -euo pipefail

    python3 /opt/pipeline/cellranger_task.py \
      --fastq_dir        "~{fastq_dir}" \
      --sample_tag       "~{sample_tag}" \
      --transcriptome    "~{transcriptome}" \
      --chemistry        "~{chemistry}" \
      --include_introns  "~{include_introns}" \
      --numproc          ~{numproc} \
      --localmem         ~{localmem} \
      --min_umi_check    ~{min_umi_check} \
      --billing_project  "~{billing_project}"
  >>>

  output {
    File    raw_h5          = read_string("raw_h5_path.txt")
    File    filtered_h5     = read_string("filtered_h5_path.txt")
    File    molecule_info   = read_string("molecule_info_path.txt")
    File    metrics_csv     = read_string("metrics_csv_path.txt")
    Int     estimated_cells = read_int("estimated_cells.txt")
    Boolean is_low_quality  = read_boolean("is_low_quality.txt")
    File    metrics_json    = "cellranger_metrics.json"
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} HDD"
    preemptible: 1
  }
}


# ---------------------------------------------------------------------------
# TASK: CiteSeqCount
# ---------------------------------------------------------------------------

task CiteSeqCount {
  input {
    String sample_tag
    String fastq_r1
    String fastq_r2
    String tag_file
    String chemistry       = "v3"
    Int    estimated_cells
    Int    numproc         = 16
    String billing_project = ""

    String docker
    Int    cpu     = 16
    Int    mem_gb  = 64
    Int    disk_gb = 200
  }

  # UMI length: v2 = 26 bp, v3 = 28 bp (from original pipeline)
  Int umil = if chemistry == "v2" then 26 else 28

  # Build gsutil prefix: add -u <project> for requester-pays buckets
  String gsutil_prefix = if billing_project != "" then "gsutil -u ~{billing_project}" else "gsutil"

  command <<<
    set -euo pipefail

    # Localize hashing FASTQs from GCS
    ~{gsutil_prefix} cp "~{fastq_r1}" ./hashing_R1.fastq.gz
    ~{gsutil_prefix} cp "~{fastq_r2}" ./hashing_R2.fastq.gz
    ~{gsutil_prefix} cp "~{tag_file}" ./tags.csv

    mkdir -p CITE-seq-Results

    CITE-seq-Count \
      -R1 ./hashing_R1.fastq.gz \
      -R2 ./hashing_R2.fastq.gz \
      -t  ./tags.csv \
      -cbf 1 -cbl 16 \
      -umif 17 -umil ~{umil} \
      -cells ~{estimated_cells} \
      -T ~{numproc} \
      -o CITE-seq-Results

    # Write the umi_count directory path for downstream use
    echo "$(pwd)/CITE-seq-Results/umi_count" > hto_count_dir.txt
  >>>

  output {
    File   hto_count_dir      = read_string("hto_count_dir.txt")
    String hto_count_dir_path = read_string("hto_count_dir.txt")
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} HDD"
    preemptible: 1
  }
}


# ---------------------------------------------------------------------------
# TASK: CellBender
# ---------------------------------------------------------------------------

task CellBender {
  input {
    String sample_tag
    File   input_h5
    Int    expected_cells
    Int    epochs          = 150
    Float  fpr             = 0.01

    String docker
    Int    gpu     = 1
    Int    cpu     = 4
    Int    mem_gb  = 32
    Int    disk_gb = 100
  }

  command <<<
    set -euo pipefail

    python /opt/pipeline/cellbender_task.py \
      --input_h5        "~{input_h5}" \
      --output_dir      "./~{sample_tag}_cb_outputs" \
      --expected_cells  ~{expected_cells} \
      --epochs          ~{epochs} \
      --fpr             ~{fpr}
  >>>

  output {
    File cb_filtered_h5 = read_string("cb_filtered_h5_path.txt")
  }

  runtime {
    docker:      docker
    cpu:         cpu
    memory:      "~{mem_gb} GB"
    disks:       "local-disk ~{disk_gb} HDD"
    gpuType:     "nvidia-tesla-t4"
    gpuCount:    gpu
    nvidiaDriverVersion: "418.87.00"
    preemptible: 0   # GPU tasks are not preemptible by default
  }
}


# ---------------------------------------------------------------------------
# TASK: PegasusPipeline
# ---------------------------------------------------------------------------

task PegasusPipeline {
  input {
    String         batch_name
    String         run_prefix
    Array[String]  sample_tags
    Array[File]    h5_files              # localized h5 files (one per sample)
    Array[String]  hto_count_dirs        # empty array when no hashing in batch
    Array[Boolean] hashing_flags
    File           mito_file
    Int            n_jobs_pg           = 10

    # QC / algorithm parameters
    Int   qc_min_umis           = 500
    Int   qc_percent_mito       = 10
    Int   qc_min_genes          = 200
    Int   dd_bst_n_iters        = 25
    Float dd_pred_pthresh       = 0.0000000000000001
    Float dd_pred_voterthresh   = 0.3
    Int   hvg_n_top             = 5000

    String docker
    Int    cpu     = 16
    Int    mem_gb  = 128
    Int    disk_gb = 200
  }

  # Build the matrix_directory JSON array: [[sample_tag, h5_path], ...]
  # and the per-batch JSON config, then call Pegasus-Pipeline.py.

  command <<<
    set -euo pipefail

    # ------------------------------------------------------------------
    # 1. Build Pegasus JSON config from WDL-provided arrays
    # ------------------------------------------------------------------
    python3 - <<'PYEOF'
import json, sys, os

sample_tags    = '~{sep="," sample_tags}'.split(",")
h5_files       = '~{sep="," h5_files}'.split(",")
hashing_flags  = '~{sep="," hashing_flags}'.split(",")    # "true"/"false" strings
hto_dirs       = '~{sep="," hto_count_dirs}'.split(",")

# matrix_directory: list of [sample_tag, h5_path]
# h5_files are Cromwell-localized local paths
matrix_directory = [[t.strip(), p.strip()] for t, p in zip(sample_tags, h5_files)]

# Batch-level hashing: True if ANY sample in batch is hashed
batch_hashing = any(f.strip().lower() == "true" for f in hashing_flags)

config = {
    "currdir":          os.getcwd(),
    "matrix_directory": matrix_directory,
    "mito_file":        "~{mito_file}",
    "n_jobs":           ~{n_jobs_pg},
    "hashing":          str(batch_hashing),
    "qc_min_umis":      ~{qc_min_umis},
    "qc_percent_mito":  ~{qc_percent_mito},
    "qc_min_genes":     ~{qc_min_genes},
    "dd_bst_n_iters":   ~{dd_bst_n_iters},
    "dd_pred_pthresh":  ~{dd_pred_pthresh},
    "dd_pred_voterthresh": ~{dd_pred_voterthresh},
    "hvg_n_top":        ~{hvg_n_top},
}

# If any sample has hashing, attach the first non-empty HTO directory
if batch_hashing:
    non_empty = [d.strip() for d in hto_dirs if d.strip()]
    if non_empty:
        config["hto_file"] = non_empty[0]

jsonfile = "~{run_prefix}_~{batch_name}_pg_input.json"
with open(jsonfile, "w") as f:
    json.dump(config, f, indent=2)

print(f"Wrote Pegasus config to: {jsonfile}")
print(json.dumps(config, indent=2))
PYEOF

    # ------------------------------------------------------------------
    # 2. Run Pegasus pipeline
    # ------------------------------------------------------------------
    JSONFILE="~{run_prefix}_~{batch_name}_pg_input.json"

    python3 /opt/pipeline/Pegasus-Pipeline.py \
      -J "${JSONFILE}" \
      -S "~{batch_name}"

    # ------------------------------------------------------------------
    # 3. Bundle all outputs into a single tar for WDL File output
    # ------------------------------------------------------------------
    tar -czf "~{batch_name}_pegasus_outputs.tar.gz" "~{batch_name}/"

    # Summary stats file path
    SUMMARY="~{batch_name}/~{batch_name}_summary_stats.txt"
    if [ ! -f "${SUMMARY}" ]; then
      # Fallback: look for any summary stats file
      SUMMARY=$(find "~{batch_name}/" -name "*summary_stats.txt" | head -1)
    fi
    echo "${SUMMARY}" > summary_path.txt
  >>>

  output {
    File output_tar    = "~{batch_name}_pegasus_outputs.tar.gz"
    File summary_stats = read_string("summary_path.txt")
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} SSD"
    preemptible: 0
  }
}
