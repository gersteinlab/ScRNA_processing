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
  make_option(c("--ahba_mat_rds"), type = "character",
              help = "Path to AHBA_mat.RDS (raw count matrix)."),
  make_option(c("--ahba_meta_rds"), type = "character",
              help = "Path to AHBA_meta_share.RDS (cell metadata)."),
  make_option(c("--ma_sestan_mat_rds"), type = "character",
              help = "Path to Ma_Sestan_mat.rds (pre-built Seurat object)."),
  make_option(c("--batch_name"), type = "character",
              help = "Batch name used as output filename prefix."),
  make_option(c("--dims"), type = "integer", default = 30,
              help = "Number of PCs for FindTransferAnchors (default: 30).")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$h5ad))
stopifnot(!is.null(opt$ahba_mat_rds))
stopifnot(!is.null(opt$ahba_meta_rds))
stopifnot(!is.null(opt$ma_sestan_mat_rds))
stopifnot(!is.null(opt$batch_name))

batch <- opt$batch_name

cat(sprintf("[1/7] Loading Ma-Sestan reference: %s\n", opt$ma_sestan_mat_rds))
seurat_obj_ref_Ma <- readRDS(opt$ma_sestan_mat_rds)
seurat_obj_ref_Ma <- UpdateSeuratObject(seurat_obj_ref_Ma)

cat(sprintf("[2/7] Loading AHBA reference matrix + metadata\n"))
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

cat(sprintf("[3/7] Loading query h5ad: %s\n", opt$h5ad))
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

cat("[4/7] SCTransform on all three objects\n")
options(future.globals.maxSize = 1000000 * 1024^2)
seurat_obj_ref_AHBA <- SCTransform(seurat_obj_ref_AHBA)
seurat_obj_ref_Ma <- SCTransform(seurat_obj_ref_Ma)
seurat_obj_query <- SCTransform(seurat_obj_query)

cat("[5/7] PCA on all three objects\n")
seurat_obj_ref_AHBA <- RunPCA(seurat_obj_ref_AHBA)
seurat_obj_ref_Ma <- RunPCA(seurat_obj_ref_Ma)
seurat_obj_query <- RunPCA(seurat_obj_query)

# ---------------------------------------------------------------------------
# Ma-Sestan label transfer
# ---------------------------------------------------------------------------
cat("[6/7] Transferring Ma-Sestan labels\n")
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
cat("[7/7] Transferring AHBA labels\n")
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
