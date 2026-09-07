version 1.0

# ============================================================================
# prep_azimuth_references.wdl
#
# ONE-TIME precomputation of the two Azimuth reference atlases used by the
# downstream workflow (terra_wdl/downstream.wdl / hybrid_azimuth_task.R).
#
# Runs SCTransform + RunPCA on the AHBA and Ma-Sestan references ONCE and
# writes the processed Seurat objects as .RDS. Stage the two outputs into
# gs://.../references/azimuth/ and point downstream.wdl's HybridAzimuth task
# at them (ahba_ref_rds / ma_sestan_ref_rds). After that, each per-batch
# HybridAzimuth run skips the (identical, expensive) reference SCTransform.
#
# Run as a standalone Terra submission with no root entity (all inputs are
# literal gs:// paths). Needs the big machine ONCE because it does the
# SCTransform fitting that dominates memory (~416 GB, same as a single
# legacy HybridAzimuth run).
# ============================================================================

workflow PrepAzimuthReferences {
  input {
    File ahba_mat_rds
    File ahba_meta_rds
    File ma_sestan_mat_rds

    String azimuth_docker            # majidfarhadloo/scrna_processing_azimuth:sha-...

    Int cpu     = 64
    Int mem_gb  = 416
    Int disk_gb = 200
  }

  call PrepReferences {
    input:
      ahba_mat_rds      = ahba_mat_rds,
      ahba_meta_rds     = ahba_meta_rds,
      ma_sestan_mat_rds = ma_sestan_mat_rds,
      docker            = azimuth_docker,
      cpu               = cpu,
      mem_gb            = mem_gb,
      disk_gb           = disk_gb
  }

  output {
    File ahba_sct_pca_rds      = PrepReferences.ahba_sct_pca_rds
    File ma_sestan_sct_pca_rds = PrepReferences.ma_sestan_sct_pca_rds
  }
}

task PrepReferences {
  input {
    File ahba_mat_rds
    File ahba_meta_rds
    File ma_sestan_mat_rds
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    Rscript /opt/pipeline/prep_azimuth_references.R \
      --ahba_mat_rds      "~{ahba_mat_rds}" \
      --ahba_meta_rds     "~{ahba_meta_rds}" \
      --ma_sestan_mat_rds "~{ma_sestan_mat_rds}" \
      --out_ahba          "AHBA_sct_pca.rds" \
      --out_ma            "Ma_Sestan_sct_pca.rds"
  >>>

  output {
    File ahba_sct_pca_rds      = "AHBA_sct_pca.rds"
    File ma_sestan_sct_pca_rds = "Ma_Sestan_sct_pca.rds"
  }

  runtime {
    docker: docker
    cpu:    cpu
    memory: "~{mem_gb} GB"
    disks:  "local-disk ~{disk_gb} SSD"
    preemptible: 0
  }
}
