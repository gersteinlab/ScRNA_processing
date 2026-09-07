version 1.0

# ============================================================================
# downstream.wdl
#
# Downstream Azimuth reference-mapping + label reconciliation + h5ad
# finalization workflow. Operates on the per-batch h5ad emitted by the
# Pegasus stage of pipeline.wdl (output: pegasus_hybrid_filtered_h5ad).
#
# Triggered as a SEPARATE Terra submission after the main pipeline
# completes. See PIPELINE_RUN_SUMMARY.md section 8 for the full design.
#
# Stages:
#   HybridAzimuth    (R + Seurat + Azimuth)   ~2-4 hr on n1-highmem-64
#   ReconcileLabels  (Python)                 ~5 min on n1-standard-2
#   FinalizeH5ad     (Python + Pegasus)       ~15-30 min on n1-highmem-16
# ============================================================================

workflow ScrnaDownstream {
  input {
    # ----------------------------------------------------------------------
    # Inputs from the main pipeline
    # ----------------------------------------------------------------------
    File pegasus_h5ad                  # pipeline.wdl: pegasus_hybrid_filtered_h5ad
    String batch_name                  # matches the Pegasus batch name

    # ----------------------------------------------------------------------
    # Azimuth references (PRE-PROCESSED once by prep_azimuth_references.wdl:
    # SCTransform + RunPCA already applied, saved as Seurat .RDS). Stage the
    # two prep outputs in GCS and point here.
    # ----------------------------------------------------------------------
    File ahba_ref_rds                  # AHBA_sct_pca.rds
    File ma_sestan_ref_rds             # Ma_Sestan_sct_pca.rds

    # ----------------------------------------------------------------------
    # Per-batch sample -> individual map
    # (from generate_sample_table.py --write_individual_map)
    # ----------------------------------------------------------------------
    File sample_to_individual_tsv

    # ----------------------------------------------------------------------
    # Knobs
    # ----------------------------------------------------------------------
    Float leiden_resolution = 4.5      # Use 4.0 for SZBD-Kellis-like data
    Int azimuth_dims = 30

    # ----------------------------------------------------------------------
    # Docker images
    # ----------------------------------------------------------------------
    String azimuth_docker              # majidfarhadloo/scrna_processing_azimuth:sha-...
    String pegasus_docker              # majidfarhadloo/scrna_processing_pegasus:sha-...

    # ----------------------------------------------------------------------
    # Compute (see PIPELINE_RUN_SUMMARY.md section 8.7)
    # ----------------------------------------------------------------------
    Int azimuth_cpu      = 64
    Int azimuth_mem_gb   = 416
    Int azimuth_disk_gb  = 500

    Int reconcile_cpu      = 2
    Int reconcile_mem_gb   = 8
    Int reconcile_disk_gb  = 50

    Int finalize_cpu      = 16
    Int finalize_mem_gb   = 104
    Int finalize_disk_gb  = 200
  }

  call HybridAzimuth {
    input:
      pegasus_h5ad       = pegasus_h5ad,
      batch_name         = batch_name,
      ahba_ref_rds       = ahba_ref_rds,
      ma_sestan_ref_rds  = ma_sestan_ref_rds,
      dims               = azimuth_dims,
      docker             = azimuth_docker,
      cpu                = azimuth_cpu,
      mem_gb             = azimuth_mem_gb,
      disk_gb            = azimuth_disk_gb
  }

  call ReconcileLabels {
    input:
      ahba_csv      = HybridAzimuth.ahba_predictions_csv,
      ma_sestan_csv = HybridAzimuth.ma_sestan_predictions_csv,
      batch_name    = batch_name,
      docker        = pegasus_docker,
      cpu           = reconcile_cpu,
      mem_gb        = reconcile_mem_gb,
      disk_gb       = reconcile_disk_gb
  }

  call FinalizeH5ad {
    input:
      pegasus_h5ad             = pegasus_h5ad,
      reconciled_csv           = ReconcileLabels.reconciled_csv,
      sample_to_individual_tsv = sample_to_individual_tsv,
      batch_name               = batch_name,
      resolution               = leiden_resolution,
      docker                   = pegasus_docker,
      cpu                      = finalize_cpu,
      mem_gb                   = finalize_mem_gb,
      disk_gb                  = finalize_disk_gb
  }

  output {
    File annotated_h5ad             = FinalizeH5ad.annotated_h5ad
    File reconciled_labels_csv      = ReconcileLabels.reconciled_csv
    File inconsistent_labels_csv    = ReconcileLabels.inconsistent_csv
    File ahba_predictions_csv       = HybridAzimuth.ahba_predictions_csv
    File ma_sestan_predictions_csv  = HybridAzimuth.ma_sestan_predictions_csv
    File ahba_umap_png              = HybridAzimuth.ahba_umap_png
    File ma_sestan_umap_png         = HybridAzimuth.ma_sestan_umap_png
  }
}

# ============================================================================
# TASK: HybridAzimuth
# Reference-mapping against AHBA + Ma-Sestan atlases using Seurat / Azimuth.
# ============================================================================
task HybridAzimuth {
  input {
    File pegasus_h5ad
    String batch_name
    File ahba_ref_rds
    File ma_sestan_ref_rds
    Int dims
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    Rscript /opt/pipeline/hybrid_azimuth_task.R \
      --h5ad             "~{pegasus_h5ad}" \
      --ahba_ref_rds     "~{ahba_ref_rds}" \
      --ma_sestan_ref_rds "~{ma_sestan_ref_rds}" \
      --batch_name       "~{batch_name}" \
      --dims             ~{dims}
  >>>

  output {
    File ahba_predictions_csv      = "~{batch_name}_Azimuth_predictions_AHBA.csv"
    File ma_sestan_predictions_csv = "~{batch_name}_Azimuth_predictions_Ma_Sestan.csv"
    File ahba_umap_png             = "~{batch_name}_Azimuth_Transferred_UMAP_AHBA.png"
    File ma_sestan_umap_png        = "~{batch_name}_Azimuth_Transferred_UMAP_Ma_Sestan.png"
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} SSD"
    preemptible: 0
  }
}

# ============================================================================
# TASK: ReconcileLabels
# Merge AHBA + Ma-Sestan predictions; keep cells where broad class agrees;
# use AHBA labels for neurons, Ma-Sestan for glia.
# ============================================================================
task ReconcileLabels {
  input {
    File ahba_csv
    File ma_sestan_csv
    String batch_name
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    python /opt/pipeline/reconcile_labels.py \
      --ahba_csv       "~{ahba_csv}" \
      --ma_sestan_csv  "~{ma_sestan_csv}" \
      --output_prefix  "~{batch_name}"
  >>>

  output {
    File reconciled_csv   = "~{batch_name}_Azimuth_predictions_Ma_Sestan_AHBA_reconcile.csv"
    File inconsistent_csv = "~{batch_name}_Azimuth_predictions_Ma_Sestan_AHBA_inconsistent.csv"
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} HDD"
    preemptible: 1
  }
}

# ============================================================================
# TASK: FinalizeH5ad
# Merge reconciled Azimuth labels into the Pegasus h5ad, attach individualID,
# re-run HVG/PCA/Harmony/Leiden/UMAP, emit final annotated h5ad.
# ============================================================================
task FinalizeH5ad {
  input {
    File pegasus_h5ad
    File reconciled_csv
    File sample_to_individual_tsv
    String batch_name
    Float resolution
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    python /opt/pipeline/finalize_h5ad_task.py \
      --h5ad                      "~{pegasus_h5ad}" \
      --azimuth_csv               "~{reconciled_csv}" \
      --sample_to_individual_tsv  "~{sample_to_individual_tsv}" \
      --batch_name                "~{batch_name}" \
      --resolution                ~{resolution} \
      --n_jobs                    ~{cpu}
  >>>

  output {
    File annotated_h5ad = "~{batch_name}_annotated.h5ad"
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} SSD"
    preemptible: 0
  }
}
