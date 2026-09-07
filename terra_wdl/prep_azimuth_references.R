#!/usr/bin/env Rscript
# prep_azimuth_references.R
#
# One-time precomputation of the two Azimuth reference atlases for the
# downstream Terra WDL workflow (terra_wdl/downstream.wdl).
#
# The per-batch task (hybrid_azimuth_task.R) previously rebuilt BOTH
# references from raw matrices and ran SCTransform + RunPCA on them on
# every single batch (~66 min and the dominant memory spike, repeated
# identically 133x across the cohort). This script performs that work
# ONCE and saves the processed Seurat objects as .RDS. hybrid_azimuth_task.R
# then readRDS()s them directly and only transforms the (small) query.
#
# The reference-processing steps below are byte-for-byte identical to the
# reference branch of hybrid_azimuth_task.R (load -> Update/CreateSeuratObject
# -> SCTransform -> RunPCA), so downstream label transfer is unchanged.
#
# Outputs (in cwd):
#   <out_ahba>   processed AHBA Seurat object    (default: AHBA_sct_pca.rds)
#   <out_ma>     processed Ma-Sestan Seurat obj  (default: Ma_Sestan_sct_pca.rds)
#
# Run once (see terra_wdl/prep_azimuth_references.wdl), then stage the two
# outputs into gs://.../references/azimuth/ and point downstream.wdl at them.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(anndata)
  library(SparseM)
  library(tidyverse)
})

option_list <- list(
  make_option(c("--ahba_mat_rds"), type = "character",
              help = "Path to AHBA_mat.RDS (raw count matrix)."),
  make_option(c("--ahba_meta_rds"), type = "character",
              help = "Path to AHBA_meta_share.RDS (cell metadata)."),
  make_option(c("--ma_sestan_mat_rds"), type = "character",
              help = "Path to Ma_Sestan_mat.rds (pre-built Seurat object)."),
  make_option(c("--out_ahba"), type = "character", default = "AHBA_sct_pca.rds",
              help = "Output path for processed AHBA object (default: %default)."),
  make_option(c("--out_ma"), type = "character", default = "Ma_Sestan_sct_pca.rds",
              help = "Output path for processed Ma-Sestan object (default: %default).")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$ahba_mat_rds))
stopifnot(!is.null(opt$ahba_meta_rds))
stopifnot(!is.null(opt$ma_sestan_mat_rds))

# Match hybrid_azimuth_task.R's future globals size ceiling for SCTransform.
options(future.globals.maxSize = 1000000 * 1024^2)

# ---------------------------------------------------------------------------
# Ma-Sestan reference: load -> UpdateSeuratObject -> SCTransform -> RunPCA
# ---------------------------------------------------------------------------
cat(sprintf("[1/4] Loading Ma-Sestan reference: %s\n", opt$ma_sestan_mat_rds))
seurat_obj_ref_Ma <- readRDS(opt$ma_sestan_mat_rds)
seurat_obj_ref_Ma <- UpdateSeuratObject(seurat_obj_ref_Ma)

# ---------------------------------------------------------------------------
# AHBA reference: build from raw matrix + metadata
# ---------------------------------------------------------------------------
cat("[2/4] Loading AHBA reference matrix + metadata\n")
rawRDS_AHBA <- readRDS(opt$ahba_mat_rds)
metaRDS_AHBA <- as.data.frame(readRDS(opt$ahba_meta_rds))
rownames(metaRDS_AHBA) <- metaRDS_AHBA$sample_id
seurat_obj_ref_AHBA <- CreateSeuratObject(
  counts = rawRDS_AHBA,
  project = "dlpfc",
  meta.data = metaRDS_AHBA,
  min.cells = 1,
  min.features = 0
)
Idents(seurat_obj_ref_AHBA) <- seurat_obj_ref_AHBA$method

# ---------------------------------------------------------------------------
# SCTransform + PCA on both references (the expensive, cacheable work)
# ---------------------------------------------------------------------------
cat("[3/4] SCTransform + PCA on AHBA reference\n")
seurat_obj_ref_AHBA <- SCTransform(seurat_obj_ref_AHBA)
seurat_obj_ref_AHBA <- RunPCA(seurat_obj_ref_AHBA)

cat("[3/4] SCTransform + PCA on Ma-Sestan reference\n")
seurat_obj_ref_Ma <- SCTransform(seurat_obj_ref_Ma)
seurat_obj_ref_Ma <- RunPCA(seurat_obj_ref_Ma)

# ---------------------------------------------------------------------------
# Persist processed references
# ---------------------------------------------------------------------------
cat(sprintf("[4/4] Writing processed references:\n  %s\n  %s\n",
            opt$out_ahba, opt$out_ma))
saveRDS(seurat_obj_ref_AHBA, opt$out_ahba)
saveRDS(seurat_obj_ref_Ma, opt$out_ma)

cat("prep_azimuth_references.R: done\n")
