# PsychENCODE / SCORCH scRNA-seq Processing Pipeline (Terra/WDL)

This repository contains the cloud-based implementation of the PsychENCODE /
SCORCH single-nucleus RNA-seq processing pipeline, running on
[Terra](https://app.terra.bio) via WDL/Cromwell on Google Cloud.

For a detailed description of the scientific workflow and individual script
logic, see `Detailed_PsychENCODE_snRNAseq_pipeline.pdf`.

---

## Pipeline Overview

The pipeline is split into two Terra workflows:

**Upstream** (`terra_wdl/pipeline.wdl`) — raw FASTQ through cell-type-agnostic
filtering:

```
Stage 1  CellRanger count        one task per sample (scattered)
Stage 1b CITE-seq-Count          one task per hashed sample (conditional)
Stage 2  CellBender              one task per sample, GPU (optional)
Stage 3  Pegasus analysis        one task per batch, produces Hybrid_filtered.h5ad
```

**Downstream** (`terra_wdl/downstream.wdl`) — Azimuth label transfer +
reconciliation + final re-clustering:

```
Stage A  HybridAzimuth           AHBA + Ma-Sestan Azimuth label transfer (R)
Stage B  ReconcileLabels         Cross-reference AHBA vs Ma-Sestan calls
Stage C  FinalizeH5ad            HVG -> PCA -> Harmony -> Leiden -> UMAP
```

**One Terra workflow submission = one batch (sample_set).**
Multiple batches are run as separate submissions.

A pilot/bypass mode is supported: set `skip_cellranger=true` and provide
`precomputed_h5_urls` (HTTPS or `gs://`) to run only the Pegasus stage on
already-CellBender-filtered `.h5` files.

---

## Repository Structure

```
ScRNA_processing/
├── .github/workflows/
│   ├── build-pegasus-image.yml          CI: Pegasus/downstream image
│   ├── build-cellranger-image.yml       CI: CellRanger image
│   └── build-azimuth-image.yml          CI: Azimuth (R) image
│
├── terra_wdl/
│   ├── pipeline.wdl                     Upstream WDL (upload to Terra)
│   ├── downstream.wdl                   Downstream WDL (upload to Terra)
│   ├── Dockerfile.pegasus               Image for Pegasus + downstream Python tasks
│   ├── Dockerfile.cellranger            Image for CellRanger stage
│   ├── Dockerfile.azimuth               Image for HybridAzimuth (R) task
│   ├── requirements.txt                 Python deps baked into Pegasus image
│   ├── cellranger_task.py               CellRanger runner
│   ├── cellbender_task.py               CellBender runner
│   ├── hybrid_azimuth_task.R            AHBA + Ma-Sestan Azimuth label transfer
│   ├── reconcile_labels.py              Cross-reference AHBA vs Ma-Sestan calls
│   ├── finalize_h5ad_task.py            Final re-cluster + write annotated h5ad
│   ├── generate_sample_table.py         Local: generate Terra data tables
│   └── requirements_local.txt           Local pip deps
│
├── Pegasus-Pipeline.py                  Pegasus analysis (baked into image)
├── AHBA_PFC_filtered.json               Marker file (baked into image)
├── Hybrid_subclass_markers.json         Marker file (baked into image)
├── mito_genes.csv                       Mitochondrial gene list (upload to GCS)
├── dataset/
│   ├── pilot_workflow_inputs.json               Example upstream inputs
│   ├── nemo_public_pilot_workflow_inputs.json   Example bypass-mode inputs
│   ├── downstream_pilot_inputs.json             Example downstream inputs
│   └── sample_to_individual_pilot.tsv           sample_tag -> individualID map
├── Detailed_PsychENCODE_snRNAseq_pipeline.pdf
└── README.md
```

---

## Prerequisites

- A [Terra](https://app.terra.bio) workspace linked to a GCP project
- NeMo FASTQ files already staged in a GCS bucket
- Locally: Python 3.9+, `pip`, `gsutil`

---

## Step 1 — Generate Terra data tables (local)

Run once per cohort to produce the two TSV files that Terra uses as data tables.

```bash
pip install -r terra_wdl/requirements_local.txt

python terra_wdl/generate_sample_table.py \
  --prelim_data   Prelim_Data.tsv \
  --nemo_metadata nemo_metadata.tsv \
  --output_dir    ./terra_tables \
  --batch_key     donor_region \
  --write_individual_map
```

**Inputs:**

| Argument | Description |
|---|---|
| `--prelim_data` | Terra `Prelim_Data` TSV exported from the NeMo workspace |
| `--nemo_metadata` | NeMo metadata TSV (mixed library- and sample-level rows) |
| `--output_dir` | Directory where output TSVs are written |
| `--batch_key` | Batch grouping: `donor_region` (default), `donor`, or `project` |
| `--write_individual_map` | Also emit `sample_to_individual.tsv` for downstream |

**Outputs:**

| File | Description |
|---|---|
| `sample_table.tsv` | Upload as the `sample` data table in Terra |
| `sample_set_table.tsv` | Upload as the `sample_set` data table in Terra |
| `sample_to_individual.tsv` | (optional) sample_tag -> individualID map for downstream |

---

## Step 2 — Build and publish Docker images

Three images are built automatically by GitHub Actions on push to `main` or
`dev` whenever their inputs change:

| Image | Trigger paths |
|---|---|
| `majidfarhadloo/scrna_processing_pegasus` | `Dockerfile.pegasus`, `requirements.txt`, `Pegasus-Pipeline.py`, `cellranger_task.py`, `reconcile_labels.py`, `finalize_h5ad_task.py`, marker JSONs |
| `majidfarhadloo/scrna_processing_cellranger` | `Dockerfile.cellranger`, `cellranger_task.py` |
| `majidfarhadloo/scrna_processing_azimuth` | `Dockerfile.azimuth`, `hybrid_azimuth_task.R` |

Each push tags `sha-<short>`, `latest`, and any release/manual tag. **Use the
`sha-<short>` tag** (printed in the Actions step summary) as the `*_docker`
input in Terra, not `latest`, so runs are reproducible.

After the first successful build, make each repository public on Docker Hub so
Terra can pull without credentials.

CellBender uses the external public Broad image
`us.gcr.io/broad-dsde-methods/cellbender:0.3.0` — no build required.

To manually build/tag a release version:

> **GitHub → Actions → select workflow → Run workflow**

---

## Step 3 — Set up Terra

1. **Upload the WDLs** — Workflows → Find a Workflow → upload
   `terra_wdl/pipeline.wdl` and `terra_wdl/downstream.wdl` as separate methods.

2. **Stage references in GCS** (once):

   ```bash
   gsutil cp mito_genes.csv                    gs://YOUR_BUCKET/refs/
   gsutil cp sample_to_individual.tsv          gs://YOUR_BUCKET/refs/
   gsutil cp AHBA_mat.RDS AHBA_meta_share.RDS  gs://YOUR_BUCKET/refs/azimuth/
   gsutil cp Ma_Sestan_mat.rds                 gs://YOUR_BUCKET/refs/azimuth/
   ```

   The CellRanger reference (`refdata-gex-GRCh38-*`) must also live in GCS.

3. **Import data tables** — Data → import `sample_table.tsv` as the `sample`
   table and `sample_set_table.tsv` as the `sample_set` table.

---

## Step 4 — Run the upstream workflow

In the Terra workflow configuration UI, set these inputs:

| Input | Value |
|---|---|
| `pegasus_docker` | `majidfarhadloo/scrna_processing_pegasus:sha-XXXXXXX` |
| `cellranger_docker` | `majidfarhadloo/scrna_processing_cellranger:sha-XXXXXXX` |
| `cellbender_docker` | `us.gcr.io/broad-dsde-methods/cellbender:0.3.0` |
| `transcriptome_gcs_path` | e.g. `gs://YOUR_BUCKET/refs/refdata-gex-GRCh38-2024-A` |
| `mito_file` | `gs://YOUR_BUCKET/refs/mito_genes.csv` |
| `run_cellbender` | `true` or `false` |
| `billing_project` | Your GCP project (needed for requester-pays buckets) |
| `skip_cellranger` | `false` for full pipeline; `true` for pilot bypass |
| `precomputed_h5_urls` | (bypass mode only) list of HTTPS / `gs://` URLs |

Select a **`sample_set`** row as the root entity, then click **Run Analysis**.
One submission per batch.

Example input JSONs are in `dataset/` — use `pilot_workflow_inputs.json` as a
starting point.

---

## Step 5 — Upstream outputs

| Output | Description |
|---|---|
| `cr_raw_h5` | Per-sample CellRanger `raw_feature_bc_matrix.h5` |
| `cr_filtered_h5` | Per-sample CellRanger `filtered_feature_bc_matrix.h5` |
| `cr_metrics_csv` | Per-sample CellRanger `metrics_summary.csv` |
| `cr_estimated_cells` | Per-sample estimated cell count |
| `cr_low_quality_flags` | Per-sample low-quality flag (median UMI < `min_umi_check`) |
| `cb_filtered_h5_list` | Per-sample CellBender filtered H5 (if `run_cellbender=true`) |
| `pegasus_output_tar` | Batch-level Pegasus results tarball |
| `pegasus_summary` | Batch-level QC summary stats text file |
| `pegasus_hybrid_filtered_h5ad` | Filtered h5ad — feed into the downstream workflow |

---

## Step 6 — Run the downstream workflow

Upload `terra_wdl/downstream.wdl` as a separate Terra method and set inputs:

| Input | Value |
|---|---|
| `pegasus_h5ad` | GCS URL of `pegasus_hybrid_filtered_h5ad` from Step 5 |
| `batch_name` | Batch identifier (matches the `sample_set` name) |
| `ahba_mat_rds` | `gs://YOUR_BUCKET/refs/azimuth/AHBA_mat.RDS` |
| `ahba_meta_rds` | `gs://YOUR_BUCKET/refs/azimuth/AHBA_meta_share.RDS` |
| `ma_sestan_mat_rds` | `gs://YOUR_BUCKET/refs/azimuth/Ma_Sestan_mat.rds` |
| `sample_to_individual_tsv` | `gs://YOUR_BUCKET/refs/sample_to_individual.tsv` |
| `azimuth_docker` | `majidfarhadloo/scrna_processing_azimuth:sha-XXXXXXX` |
| `pegasus_docker` | Same tag as Step 4 |
| `leiden_resolution` | `4.5` (use `4.0` for SZBD-Kellis batches) |
| `azimuth_dims` | `30` |

**Downstream outputs:**

| Output | Description |
|---|---|
| `{batch}_Azimuth_predictions_AHBA.csv` | Per-cell AHBA label + score |
| `{batch}_Azimuth_predictions_Ma_Sestan.csv` | Per-cell Ma-Sestan label + score |
| `{batch}_Azimuth_Transferred_UMAP_{AHBA,Ma_Sestan}.png` | UMAP figures |
| `{batch}_Azimuth_predictions_Ma_Sestan_AHBA_reconcile.csv` | Reconciled calls |
| `{batch}_Azimuth_predictions_Ma_Sestan_AHBA_inconsistent.csv` | Disagreements |
| `{batch}_annotated.h5ad` | Final h5ad with `obs.azimuth`, `obs.subclass`, `obs.individualID` |

`dataset/downstream_pilot_inputs.json` is an example input set.

---

## Docker images

```
majidfarhadloo/scrna_processing_pegasus:{sha-tag}       Pegasus + downstream Python tasks
majidfarhadloo/scrna_processing_cellranger:{sha-tag}    CellRanger stage
majidfarhadloo/scrna_processing_azimuth:{sha-tag}       HybridAzimuth (R + Seurat 5 + Azimuth)
us.gcr.io/broad-dsde-methods/cellbender:0.3.0           External Broad public image
```

The Pegasus image contains Python 3.9 with all deps from
`terra_wdl/requirements.txt`, plus `Pegasus-Pipeline.py`, `cellranger_task.py`,
`cellbender_task.py`, `reconcile_labels.py`, `finalize_h5ad_task.py`, and the
two marker JSONs at `/opt/pipeline/`.

The Azimuth image is based on `bioconductor/bioconductor_docker:RELEASE_3_18`
with Seurat 5, Azimuth, SeuratDisk, and anndata preinstalled.
