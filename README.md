# PsychENCODE scRNA-seq Processing Pipeline (Terra/WDL)

This repository contains the cloud-based implementation of the PsychENCODE
single-nucleus RNA-seq processing pipeline, running on
[Terra](https://app.terra.bio) via WDL/Cromwell on Google Cloud.

For a detailed description of the scientific workflow and individual script
logic, see `Detailed_PsychENCODE_snRNAseq_pipeline.pdf`.

---

## Pipeline Overview

The pipeline processes snRNA-seq data from the NeMo portal (PsychENCODE /
SCORCH) through three stages:

```
Stage 1  CellRanger count        one task per sample (scattered)
Stage 1b CITE-seq-Count          one task per hashed sample (conditional)
Stage 2  CellBender              one task per sample, GPU (optional)
Stage 3  Pegasus analysis        one task per batch
```

**One Terra workflow submission = one batch (sample_set).**
Multiple batches are run as separate submissions.

---

## Repository Structure

```
ScRNA_processing/
├── .github/
│   └── workflows/
│       └── build-pegasus-image.yml   CI/CD: auto-builds Docker image on push
│
├── terra_wdl/
│   ├── pipeline.wdl                  WDL workflow — upload this to Terra
│   ├── Dockerfile.pegasus            Docker image for Pegasus/CellRanger tasks
│   ├── requirements.txt              Python deps baked into Docker image
│   ├── cellranger_task.py            CellRanger runner (baked into image)
│   ├── cellbender_task.py            CellBender runner (baked into image)
│   ├── generate_sample_table.py      Local script: generates Terra data tables
│   └── requirements_local.txt        Local pip deps for generate_sample_table.py
│
├── Pegasus-Pipeline.py               Pegasus analysis script (baked into image)
├── AHBA_PFC_filtered.json            Marker file (baked into image)
├── Hybrid_subclass_markers.json      Marker file (baked into image)
├── mito_genes.csv                    Mitochondrial gene list (upload to GCS)
├── Detailed_PsychENCODE_snRNAseq_pipeline.pdf   Full pipeline documentation
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
  --batch_key     donor_region
```

**Inputs:**

| Argument | Description |
|---|---|
| `--prelim_data` | Terra `Prelim_Data` TSV exported from the NeMo workspace (one row per FASTQ file) |
| `--nemo_metadata` | NeMo metadata TSV (mixed library- and sample-level rows) |
| `--output_dir` | Directory where output TSVs are written |
| `--batch_key` | Batch grouping strategy: `donor_region` (default), `donor`, or `project` |

**Outputs:**

| File | Description |
|---|---|
| `sample_table.tsv` | One row per sample — upload as the `sample` data table in Terra |
| `sample_set_table.tsv` | One row per batch — upload as the `sample_set` data table in Terra |

---

## Step 2 — Build and publish the Docker image

The Docker image is built automatically by GitHub Actions whenever any of the
following files change on `main`:

```
terra_wdl/Dockerfile.pegasus
terra_wdl/requirements.txt
terra_wdl/cellranger_task.py
terra_wdl/cellbender_task.py
Pegasus-Pipeline.py
AHBA_PFC_filtered.json
Hybrid_subclass_markers.json
```

After the **first** successful build, make the package public so Terra can pull
it without credentials:

> **GitHub → gersteinlab org → Packages → `scrna_processing/pegasus` →
> Package settings → Change visibility → Public**

The Actions step summary prints the exact image URI to use in Terra, e.g.:

```
ghcr.io/gersteinlab/scrna_processing/pegasus:sha-a3f2c1b
```

To manually trigger a build or tag a release version:

> **GitHub → Actions → Build and Push Pegasus Docker Image → Run workflow**

Enter an optional version tag (e.g. `1.1`) in the input field.

---

## Step 3 — Set up Terra

1. **Upload the WDL** — in your Terra workspace go to
   Workflows → Find a Workflow → upload `terra_wdl/pipeline.wdl`

2. **Upload `mito_genes.csv` to GCS**

   ```bash
   gsutil cp mito_genes.csv gs://YOUR_BUCKET/refs/mito_genes.csv
   ```

3. **Import data tables** — go to Data in your Terra workspace and import:
   - `sample_table.tsv` as the `sample` table
   - `sample_set_table.tsv` as the `sample_set` table

---

## Step 4 — Configure and run

In the Terra workflow configuration UI, set these inputs:

| Input | Value |
|---|---|
| `pegasus_docker` | `ghcr.io/gersteinlab/scrna_processing/pegasus:sha-XXXXXXX` (from Actions summary) |
| `transcriptome_gcs_path` | GCS path to CellRanger reference, e.g. `gs://YOUR_BUCKET/refs/refdata-gex-GRCh38-2020-A` |
| `mito_file` | `gs://YOUR_BUCKET/refs/mito_genes.csv` |
| `run_cellbender` | `true` or `false` |

Select a **`sample_set`** row as the root entity, then click **Run Analysis**.
Repeat one submission per batch.

---

## Step 5 — Outputs

| Output | Description |
|---|---|
| `cr_raw_h5` | Per-sample CellRanger `raw_feature_bc_matrix.h5` |
| `cr_filtered_h5` | Per-sample CellRanger `filtered_feature_bc_matrix.h5` |
| `cr_metrics_csv` | Per-sample CellRanger `metrics_summary.csv` |
| `cr_estimated_cells` | Per-sample estimated cell count |
| `cr_low_quality_flags` | Per-sample low-quality flag (median UMI < 1000) |
| `cb_filtered_h5_list` | Per-sample CellBender filtered H5 (when `run_cellbender=true`) |
| `pegasus_output_tar` | Batch-level Pegasus results tarball |
| `pegasus_summary` | Batch-level QC summary stats text file |

---

## Docker image

```
ghcr.io/gersteinlab/scrna_processing/pegasus:{tag}
```

The image contains:
- Python 3.9 with all dependencies from `terra_wdl/requirements.txt`
- `Pegasus-Pipeline.py`, `cellranger_task.py`, `cellbender_task.py` at `/opt/pipeline/`
- `AHBA_PFC_filtered.json`, `Hybrid_subclass_markers.json` at `/opt/pipeline/`
