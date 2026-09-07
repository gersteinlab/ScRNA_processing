#!/usr/bin/env Rscript
# hybrid_azimuth_task.R
#
# Runs Azimuth-style reference mapping against TWO atlases (AHBA + Ma-Sestan)
# for the downstream Terra WDL workflow (terra_wdl/downstream.wdl).
#
# Refactored from downstream/Hybrid_Azimuth.R:
#   - CLI args via optparse replace hardcoded `data_prefix` literal.
#   - Removes `use_condaenv("test2")` (system Python in container).
#   - Reads .RDS reference paths from CLI args, not relative paths.
#   - Writes outputs to current working directory; the WDL `output {}`
#     block captures them via globs.
#
# Reference caching (see prep_azimuth_references.R / .wdl):
#   The two reference atlases are now PRE-processed (SCTransform + RunPCA)
#   once by prep_azimuth_references.R and passed in here as ready-to-use
#   Seurat .RDS objects (--ahba_ref_rds / --ma_sestan_ref_rds). This task
#   no longer rebuilds/SCTransforms the references on every batch; it only
#   loads them and transforms the (small) query. The objects handed to
#   FindTransferAnchors are identical to the previous in-task build, so
#   label-transfer output is unchanged.
#
# Outputs (in cwd):
#   {batch_name}_Azimuth_predictions_AHBA.csv
#   {batch_name}_Azimuth_predictions_Ma_Sestan.csv
#   {batch_name}_Azimuth_Transferred_UMAP_AHBA.png
#   {batch_name}_Azimuth_Transferred_UMAP_Ma_Sestan.png

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
  make_option(c("--h5ad"), type = "character",
              help = "Path to the Pegasus Hybrid_filtered h5ad."),
  make_option(c("--ahba_ref_rds"), type = "character",
              help = "Path to the PRE-processed AHBA reference (SCTransform + PCA) .RDS from prep_azimuth_references.R."),
  make_option(c("--ma_sestan_ref_rds"), type = "character",
              help = "Path to the PRE-processed Ma-Sestan reference (SCTransform + PCA) .RDS from prep_azimuth_references.R."),
  make_option(c("--batch_name"), type = "character",
              help = "Batch name used as output filename prefix."),
  make_option(c("--dims"), type = "integer", default = 30,
              help = "Number of PCs for FindTransferAnchors (default: 30).")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$h5ad))
stopifnot(!is.null(opt$ahba_ref_rds))
stopifnot(!is.null(opt$ma_sestan_ref_rds))
stopifnot(!is.null(opt$batch_name))

batch <- opt$batch_name

cat(sprintf("[1/6] Loading pre-processed references:\n  AHBA: %s\n  Ma-Sestan: %s\n",
            opt$ahba_ref_rds, opt$ma_sestan_ref_rds))
seurat_obj_ref_AHBA <- readRDS(opt$ahba_ref_rds)
seurat_obj_ref_Ma <- readRDS(opt$ma_sestan_ref_rds)

cat(sprintf("[2/6] Loading query h5ad: %s\n", opt$h5ad))
seurat_obj_query <- read_h5ad(opt$h5ad)
cat(sprintf("  Query raw dim: %d x %d\n",
            dim(seurat_obj_query$raw$X)[1], dim(seurat_obj_query$raw$X)[2]))

count.data <- t(seurat_obj_query$raw$X)
count.data <- as(count.data, "matrix.csr")
count.data <- as(count.data, "dgCMatrix")
colnames(count.data) <- row.names(seurat_obj_query$obs)
row.names(count.data) <- row.names(seurat_obj_query$var)
seurat_obj_query <- SeuratObject::CreateSeuratObject(
  counts = count.data,
  meta.data = seurat_obj_query$obs
)

cat("[3/6] SCTransform on query\n")
options(future.globals.maxSize = 1000000 * 1024^2)
seurat_obj_query <- SCTransform(seurat_obj_query)

cat("[4/6] PCA on query\n")
seurat_obj_query <- RunPCA(seurat_obj_query)

# ---------------------------------------------------------------------------
# Ma-Sestan label transfer
# ---------------------------------------------------------------------------
cat("[5/6] Transferring Ma-Sestan labels\n")
transfer_anchors_Ma <- FindTransferAnchors(
  reference = seurat_obj_ref_Ma,
  query = seurat_obj_query,
  dims = 1:opt$dims,
  reduction = "pcaproject"
)
predictions_Ma <- TransferData(
  anchorset = transfer_anchors_Ma,
  refdata = seurat_obj_ref_Ma$subclass,
  dims = 1:opt$dims
)
seurat_obj_query <- AddMetaData(seurat_obj_query, metadata = predictions_Ma)

write.csv(
  seurat_obj_query$predicted.id,
  sprintf("./%s_Azimuth_predictions_Ma_Sestan.csv", batch),
  quote = FALSE, row.names = TRUE
)

seurat_obj_query <- RunUMAP(seurat_obj_query, dims = 1:opt$dims)
Idents(seurat_obj_query) <- seurat_obj_query$predicted.id
png(file = sprintf("./%s_Azimuth_Transferred_UMAP_Ma_Sestan.png", batch),
    width = 1200, height = 1000)
print(DimPlot(seurat_obj_query, reduction = "umap",
              label = TRUE, pt.size = 0.5) + NoLegend())
dev.off()

# ---------------------------------------------------------------------------
# AHBA label transfer
# ---------------------------------------------------------------------------
cat("[6/6] Transferring AHBA labels\n")
transfer_anchors_AHBA <- FindTransferAnchors(
  reference = seurat_obj_ref_AHBA,
  query = seurat_obj_query,
  dims = 1:opt$dims,
  reduction = "pcaproject"
)
predictions_AHBA <- TransferData(
  anchorset = transfer_anchors_AHBA,
  refdata = seurat_obj_ref_AHBA$within_area_subclass,
  dims = 1:opt$dims
)
seurat_obj_query <- AddMetaData(seurat_obj_query, metadata = predictions_AHBA)

write.csv(
  seurat_obj_query$predicted.id,
  sprintf("./%s_Azimuth_predictions_AHBA.csv", batch),
  quote = FALSE, row.names = TRUE
)

seurat_obj_query <- RunUMAP(seurat_obj_query, dims = 1:opt$dims)
Idents(seurat_obj_query) <- seurat_obj_query$predicted.id
png(file = sprintf("./%s_Azimuth_Transferred_UMAP_AHBA.png", batch),
    width = 1200, height = 1000)
print(DimPlot(seurat_obj_query, reduction = "umap",
              label = TRUE, pt.size = 0.5) + NoLegend())
dev.off()

cat("hybrid_azimuth_task.R: done\n")
