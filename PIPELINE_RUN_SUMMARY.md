# scRNA-seq Pipeline — Run Summary & Input Reference

A snapshot of what we've run, what inputs each run mode requires, and how
to set up new cohorts. Companion to `CLAUDE.md` (which holds architecture
and debugging history).

---

## 1. Latest Runs

### Full Pipeline Run (CellRanger -> CellBender -> Pegasus)
- **Batch**: `insular_cortex_pilot`
- **Mode**: `skip_cellranger: false`, `run_cellbender: true`
- **Samples**: 2

  | Sample tag | Condition | FASTQ dir | Lane |
  |---|---|---|---|
  | `HCTBA_CTR_INS1_HNT_S2_L003` | control | `gs://scorch-tmp-jwbzm77-nsu/.../HCTBA_CTR_INS1_HNT_S2_L003/` | L003 |
  | `HCCLC_OUD_INS1_HNT_S11_L002` | OUD | `gs://scorch-tmp-jwbzm77-nsu/.../HCCLC_OUD_INS1_HNT_S11_L002/` | L002 |

- **FASTQ files**: 6 total (3 per sample: R1 + R2 + I1, single-lane)

### Pegasus-Only Run (skip CellRanger + CellBender)
- **Batch**: `pfc_pilot_batch`
- **Mode**: `skip_cellranger: true`, `run_cellbender: false`
- **Samples**: 13 (precomputed CellBender-filtered `.h5` from NeMo)
  - Kellis cohort (9): `HCcFW_07`, `HCcHO_05`, `HCcLS_01`, `HCcOE_06`,
    `HCcOV_03`, `HCcPN_02`, `HccUD_10`, `HCcWK_12`, `HCtPA_16`
    (all `_BA9-FC_21A`)
  - Spudich cohort (4): `HCCCX_PFC`, `HCCDE_PFC`, `NIH967_PFC`, `NIH1055_PFC`
- **FASTQ files**: 0 (h5 inputs only)

---

## 2. FASTQ File Anatomy

### What R1, R2, I1 mean
Not "regions" — they are distinct **Illumina sequencing reads** produced per
DNA cluster on the flow cell. For 10x 3' GEX v3 / v3.1 (the pilot's chemistry):

| File | Name | Length | Contents | Used for |
|---|---|---|---|---|
| **R1** | Read 1 | 28 bp | 16 bp **cell barcode** + 12 bp **UMI** | Demultiplex into cells; PCR-dedup via UMI |
| **R2** | Read 2 | 90 bp | **cDNA sequence** | Aligned to transcriptome (STAR inside CellRanger) |
| **I1** | Index 1 | 8 bp | **Sample index (i7)** | Demultiplex multiplexed lane into per-sample FASTQs (used pre-CellRanger) |

R1, R2, and I1 from the same line number all originate from the **same cluster**
— that's how CellRanger ties cell barcode (R1), transcript (R2), and sample (I1)
together. CellRanger strictly needs only R1 + R2; I1 is informational by the
time the FASTQs reach the pipeline.

Other variants you may encounter:
- **I2** — i5 index, for dual-indexed libraries (NovaSeq)
- **R3** — appears in 10x Multiome / ATAC chemistries (R1+R3 paired, R2 = barcode)

### What "single-lane" means
A flow cell has multiple physical **lanes** (NovaSeq S4 = 4, NextSeq = 4,
MiSeq = 1). A library can be loaded onto one or more lanes for additional
sequencing depth. The lane number appears in the FASTQ filename:

```
{Sample}_S{SampleNumber}_L{LaneNumber}_{R1|R2|I1}_001.fastq.gz
```

Both pilot samples appear on a single lane (`L003`, `L002` respectively),
so each yields **3 FASTQs** (R1 + R2 + I1).

A 4-lane sample would yield **12 FASTQs**:

```
HCTBA_CTR_INS1_HNT_S2_L001_R1_001.fastq.gz
HCTBA_CTR_INS1_HNT_S2_L001_R2_001.fastq.gz
HCTBA_CTR_INS1_HNT_S2_L001_I1_001.fastq.gz
HCTBA_CTR_INS1_HNT_S2_L002_R1_001.fastq.gz
... (and so on for L003, L004)
```

CellRanger automatically merges reads across lanes for the same sample by
matching the `--sample` prefix against filenames in the FASTQ directory.
**No WDL changes are required for multi-lane samples** — they go in the
same GCS folder, and the `fastq_sample_prefixes` input identifies them.

---

## 3. Naming Convention Pitfall (encountered in pilot)

CellRanger 10.0 enforces strict 10x FASTQ naming:

```
{Sample}_S{N}_L{NNN}_R[12I]\d*_001.fastq.gz
```

The `--sample` flag passed to `cellranger count` must match the `{Sample}`
portion of the FASTQ filenames **exactly** (CellRanger does a literal prefix
match against files in the `--fastqs` directory). Mismatches produce errors
like:

```
error: Could not find FASTQs matching --sample=...
```

### What we hit in the pilot
Our `sample_tags` include the lane / sample-number suffix, e.g.:

```
HCTBA_CTR_INS1_HNT_S2_L003
```

But the actual FASTQ files are named:

```
HCTBA_CTR_INS1_HNT_S2_L003_R1_001.fastq.gz   <- {Sample}=HCTBA_CTR_INS1_HNT
HCTBA_CTR_INS1_HNT_S2_L003_R2_001.fastq.gz
HCTBA_CTR_INS1_HNT_S2_L003_I1_001.fastq.gz
```

The `_S2_L003_R1_001.fastq.gz` segment is the standard Illumina suffix, and
the actual `{Sample}` prefix CellRanger expects is `HCTBA_CTR_INS1_HNT` —
**not** `HCTBA_CTR_INS1_HNT_S2_L003`.

### Resolution: `fastq_sample_prefixes` input
A separate Array input was added to the WDL to **decouple sample identity
from CellRanger's `--sample` argument**:

| Input | Used as | Example |
|---|---|---|
| `sample_tags[i]` | CellRanger `--id`, Pegasus metadata key, output dir name | `HCTBA_CTR_INS1_HNT_S2_L003` |
| `fastq_sample_prefixes[i]` | CellRanger `--sample` (must match FASTQ filename prefix) | `HCTBA_CTR_INS1_HNT` |

`cellranger_task.py` also auto-detects the prefix by regex-scanning the
FASTQ directory as a robust fallback if `fastq_sample_prefixes[i]` is empty.

### Rule of thumb when adding new samples
1. Inspect filenames in the FASTQ dir: `gsutil ls gs://.../sample_dir/`
2. Strip the standard Illumina suffix `_S\d+_L\d+_R[12I]\d*_001.fastq.gz`
3. Whatever remains = `fastq_sample_prefixes[i]`
4. Choose `sample_tags[i]` independently — must be globally unique across
   the batch and stable across re-sequencing (i.e., **avoid baking in
   `_S{N}_L{NNN}` if you may add lanes later**)
5. Multi-lane samples: same prefix appears with multiple `_L{NNN}_`;
   CellRanger merges them automatically — no special handling

> **Caveat noted in pilot**: Current `sample_tags` (`_S2_L003`, `_S11_L002`)
> bake the lane number into the sample identifier. If a biological sample is
> later re-sequenced on additional lanes, the tag will no longer uniquely
> identify the *biological* sample — only that specific lane's run. Drop
> `_S{N}_L{NNN}` from sample tags when scaling beyond the pilot.

---

## 4. Required Inputs by Run Mode

### 4a. Full Pipeline Mode (`skip_cellranger: false`)

#### Per-sample arrays (all index-aligned, same length)
| Input | Type | Purpose |
|---|---|---|
| `fastq_dirs` | `Array[String]` | GCS path to each sample's FASTQ directory |
| `sample_tags` | `Array[String]` | Unique sample ID (CellRanger `--id`, Pegasus metadata) |
| `fastq_sample_prefixes` | `Array[String]` | FASTQ filename prefix for CellRanger `--sample` |
| `chemistry_flags` | `Array[String]` | CellRanger chemistry (`"auto"` or e.g. `"SC3Pv3"`) |
| `include_introns_flags` | `Array[String]` | `"true"` for snRNA-seq (nuclei), `"false"` for scRNA-seq |
| `hashing_flags` | `Array[Boolean]` | Whether sample is HTO/CITE-seq hashed |
| `hashing_fastq_r1` | `Array[String]` | HTO R1 FASTQ (or `""`) |
| `hashing_fastq_r2` | `Array[String]` | HTO R2 FASTQ (or `""`) |
| `hashing_tag_files` | `Array[String]` | HTO tag list CSV (or `""`) |

#### FASTQ file requirements
Each `fastq_dirs[i]` must contain 10x-naming-compliant files:
- `{prefix}_S{N}_L{NNN}_R1_001.fastq.gz`
- `{prefix}_S{N}_L{NNN}_R2_001.fastq.gz`
- `{prefix}_S{N}_L{NNN}_I1_001.fastq.gz` (optional but typical)

Single-lane -> 3 FASTQs/sample. Multi-lane -> 3 x N_lanes FASTQs/sample.

#### Shared references
| Input | Description |
|---|---|
| `transcriptome_gcs_path` | CellRanger reference dir (e.g., `refdata-gex-GRCh38-2024-A`) |
| `mito_file` | `mito_genes.csv` for Pegasus QC |

#### Batch / run config
| Input | Description |
|---|---|
| `batch_name` | Pegasus batch identifier |
| `run_prefix` | Output filename prefix |
| `run_cellbender` | `true` to run CellBender between CellRanger and Pegasus |
| `billing_project` | Required for requester-pays buckets (e.g., `terra-4ea165c8`) |

#### Docker images (pin by SHA, not `:latest`)
- `cellranger_docker`, `cellbender_docker`, `pegasus_docker`

#### Optional QC / Pegasus knobs (with defaults)
`qc_min_umis`, `qc_percent_mito`, `qc_min_genes`, `n_jobs_pg`, `hvg_n_top`,
`dd_bst_n_iters`, `dd_pred_pthresh`, `dd_pred_voterthresh`,
`pegasus_cpu/mem_gb/disk_gb`

#### Bypass-mode keys (must be set to "off")
- `skip_cellranger: false`
- `precomputed_h5_urls: []`

---

### 4b. Pegasus-Only Mode (`skip_cellranger: true`)

#### Required
| Input | Type | Purpose |
|---|---|---|
| `skip_cellranger` | `Boolean` | Must be `true` |
| `run_cellbender` | `Boolean` | Typically `false` (h5 already CellBender-filtered) |
| `sample_tags` | `Array[String]` | Sample IDs (still needed for Pegasus metadata) |
| `hashing_flags` | `Array[Boolean]` | Per-sample hashing flag |
| `precomputed_h5_urls` | `Array[String]` | HTTPS or `gs://` URL per sample (CellBender-filtered `.h5`) |
| `mito_file` | `String` | Pegasus QC reference |
| `batch_name`, `run_prefix` | `String` | Output identifiers |
| `pegasus_docker` | `String` | Pinned image SHA |

#### Not needed
- `fastq_dirs`, `fastq_sample_prefixes`, `chemistry_flags`,
  `include_introns_flags` — no FASTQ stage runs
- `transcriptome_gcs_path` — CellRanger skipped
- `cellranger_docker`, `cellbender_docker` — those tasks are skipped
- `billing_project` — only required if `precomputed_h5_urls` point to
  requester-pays GCS buckets

---

## 5. Reference: Input Files in This Repo

### Full Pipeline (insular cortex pilot)
- `dataset/pilot_workflow_inputs.json`
- `dataset/pilot_sample_table.tsv` — Terra `sample` entity table
- `dataset/pilot_sample_set_membership.tsv` — Terra `sample_set` membership

### Pegasus-Only (NeMo PFC pilot)
- `dataset/nemo_public_pilot_workflow_inputs.json`
- `dataset/nemo_public_pilot_sample_table.tsv`
- `dataset/nemo_public_pilot_sample_set_membership.tsv`

### Generating new tables
`terra_wdl/generate_sample_table.py` builds the Terra entity TSVs from
upstream metadata (`Prelim_Data.tsv` + `nemo_metadata.tsv`). See `CLAUDE.md`
for the `--batch_key` strategy options.

#### Outputs (per cohort):
- `{prefix}_sample_table.tsv` — Terra `sample` entity table
- `{prefix}_sample_set_table.tsv` — Terra `sample_set` membership
- `{prefix}_sample_to_individual.tsv` — **NEW**: 2-column map
  (`sample_tag` -> `individualID`) emitted when `--write-individual-map`
  is passed. Required input to the downstream Azimuth workflow (§8).

#### New flag:
- `--write-individual-map` (bool): emits the `sample_to_individual.tsv`
  by reading `subject_id` / `individualID` from `Prelim_Data.tsv` or
  `nemo_metadata.tsv` and aligning to `sample_tag`.

---

## 6. Currently Pinned Docker Images

| Stage | Image | Tag |
|---|---|---|
| CellRanger | `majidfarhadloo/scrna_processing_cellranger` | `sha-d47ded4` |
| CellBender | `us.gcr.io/broad-dsde-methods/cellbender` | `0.3.0` |
| Pegasus (full pilot + NeMo pilot + downstream) | `majidfarhadloo/scrna_processing_pegasus` | `sha-4617b9a` |
| Azimuth (downstream) | `majidfarhadloo/scrna_processing_azimuth` | *first build in progress* |

> Pegasus pins unified to `sha-4617b9a` (post filter-step rebuild). Both
> the full-pipeline (`pilot_workflow_inputs.json`) and NeMo-only pilot
> (`nemo_public_pilot_workflow_inputs.json`) now use this tag.

---

## 7. Quick Start: Adding a New Cohort

1. Stage FASTQs in GCS — one folder per sample, 10x-named files
2. Build index-aligned arrays (`fastq_dirs`, `sample_tags`,
   `fastq_sample_prefixes`, chemistry/intron/hashing flags) — see §3 for
   how to derive `fastq_sample_prefixes` from filenames
3. Stage references in GCS (`transcriptome_gcs_path`, `mito_file`) — one-time
4. Update workflow JSON (or generate Terra tables via
   `terra_wdl/generate_sample_table.py --write-individual-map`) with
   pinned docker SHAs and `batch_name`
5. Submit on Terra with **call caching unchecked**, **delete intermediates
   unchecked** (until pipeline is fully validated)
6. Verify CellRanger h5 -> CellBender filtered h5 -> Pegasus tarball appear
   in workspace bucket
7. **(Optional)** Run downstream Azimuth workflow as a separate Terra
   submission against the Pegasus task's `hybrid_filtered_h5ad` output
   (see §8). Requires `sample_to_individual.tsv` from step 4.

---

## 8. Downstream Analysis Workflow

A **separate Terra workflow** (`terra_wdl/downstream.wdl`) that runs
Azimuth reference-mapping, label reconciliation, and h5ad finalization
on the per-batch h5ad emitted by the Pegasus stage of the main pipeline.

### 8.1 Why a Separate WDL (Decoupling Rationale)

- **Independent compute profile**: Azimuth needs ~400 GB RAM (R/Seurat
  heavy); fundamentally different from Pegasus stage's GPU+Python profile
- **Re-runnable on existing outputs**: can iterate on Leiden resolution,
  atlas reference, or label reconciliation logic without re-running
  CellRanger or CellBender
- **Image isolation**: R + Seurat + Azimuth (~5 GB) kept separate from
  the Python-only Pegasus image
- **Matches the original Slurm pipeline structure**: each step was its
  own job — no scientific reason to fuse them
- **Promotion path**: if/when cohort cadence stabilizes, downstream can
  be inlined into `pipeline.wdl` as a final WDL `call` (see §9 follow-ups)

### 8.2 Locked Design Decisions

These were locked during planning — do not re-litigate without
explicit reason.

| # | Decision |
|---|---|
| 0.1 | Stage Azimuth references in **`gs://fc-dde46beb-1862-45cf-92bc-06c767e05710/references/azimuth/`** (Terra workspace bucket). **Caveat**: this bucket is auto-provisioned by Terra (`fc-<UUID>` prefix = "Firecloud") and is **deleted if the workspace is deleted**. For long-term persistence, plan to migrate to a dedicated bucket (see §9). |
| 0.2 | Use **generic `sample_to_individual.tsv`** (2-column: `sample_tag`, `individualID`) instead of hardcoded study branches in `finalize_h5ad.py`. Generated by `generate_sample_table.py --write-individual-map`. |
| 0.3 | **Add filter step to `Pegasus-Pipeline.py`** so it emits `{batch}_Hybrid_filtered.h5ad` (drops clusters labeled only as `"Inhibitory neuron"`, `"Excitatory neuron"`, or bare integer cluster numbers). Downstream WDL starts directly at Azimuth — no clustering duplication. |
| 0.4 | **Keep both atlases** (AHBA + Ma-Sestan hybrid scheme). Accept the ~400 GB RAM Azimuth runtime cost. |
| 0.5 | **Option (ii)**: explicit `File hybrid_filtered_h5ad` added to Pegasus task's WDL `output {}` via `glob()`. Downstream WDL takes `File pegasus_h5ad` directly — no untar logic needed. |
| F3 | **Auto-generate** `sample_to_individual.tsv` from `generate_sample_table.py` via new `--write-individual-map` flag (documented in §5). |
| F4 | Azimuth machine type: **`n1-highmem-64`** (64 vCPU / 416 GB RAM), **non-preemptible**. Matches Slurm's 400 GB request; ~2-4 hr runtime; ~$8-15 per Azimuth run. |

### 8.3 Workflow Stages

```
[Pegasus stage (main pipeline)]
    | emits hybrid_filtered_h5ad
    v
[downstream.wdl]
    Task 1: HybridAzimuth      (R + Seurat + Azimuth, n1-highmem-64, ~2-4 hr)
            -> AHBA_predictions.csv, Ma_Sestan_predictions.csv, UMAP PNGs
    Task 2: ReconcileLabels    (Python, n1-standard-2, ~5 min)
            -> reconciled.csv, inconsistent.csv
    Task 3: FinalizeH5ad       (Python + Pegasus, n1-highmem-16, ~15-30 min)
            -> annotated.h5ad with obs.azimuth, obs.subclass, obs.individualID
```

### 8.4 Required Inputs

| Input | Type | Source |
|---|---|---|
| `pegasus_h5ad` | `File` | Pegasus task's `hybrid_filtered_h5ad` output |
| `batch_name` | `String` | Matches the Pegasus batch name |
| `ahba_mat_rds` | `File` | `gs://.../references/azimuth/AHBA_mat.RDS` |
| `ahba_meta_rds` | `File` | `gs://.../references/azimuth/AHBA_meta_share.RDS` |
| `ma_sestan_mat_rds` | `File` | `gs://.../references/azimuth/Ma_Sestan_mat.rds` |
| `sample_to_individual_tsv` | `File` | From `generate_sample_table.py --write-individual-map` |
| `azimuth_docker` | `String` | Pinned `majidfarhadloo/scrna_processing_azimuth:sha-...` |
| `pegasus_docker` | `String` | Reused for Python tasks (reconcile, finalize) |
| `leiden_resolution` | `Float` | Re-clustering resolution in FinalizeH5ad (default 4.5; 4.0 for SZBD-Kellis-like data) |

### 8.5 Reference Data Location

Stage once (one-time `gsutil cp`):

```
gs://fc-dde46beb-1862-45cf-92bc-06c767e05710/references/azimuth/
    AHBA_mat.RDS              (~1.7 GB)
    AHBA_meta_share.RDS       (~0.7 MB)
    Ma_Sestan_mat.rds         (~2.0 GB)
```

> **Terra workspace bucket caveat**: The `fc-<UUID>` bucket is owned by
> Terra's Firecloud infrastructure and is **automatically deleted when
> the workspace is deleted**. For pilot work this is fine. For production
> / long-lived references, migrate to a dedicated non-Terra GCS bucket
> in `terra-4ea165c8` (or another GCP project) — see §9 follow-ups.

### 8.6 Docker Images

| Stage | Image | Notes |
|---|---|---|
| HybridAzimuth | `majidfarhadloo/scrna_processing_azimuth` (NEW) | rocker/r-ver:4.3.2 + Seurat 5 + Azimuth + SeuratDisk + anndata + reticulate + Python helpers (~5 GB) |
| ReconcileLabels | `majidfarhadloo/scrna_processing_pegasus` (reused) | Lightweight Python; reuses existing image |
| FinalizeH5ad | `majidfarhadloo/scrna_processing_pegasus` (reused) | Needs Pegasus + scipy |

### 8.7 Compute Requirements

| Task | Machine | Disk | Preemptible | Est. cost (pilot) |
|---|---|---|---|---|
| HybridAzimuth | `n1-highmem-64` (64 vCPU, 416 GB) | 500 GB | No | $4-8 |
| ReconcileLabels | `n1-standard-2` (2 vCPU, 7.5 GB) | 50 GB | Yes | <$0.10 |
| FinalizeH5ad | `n1-highmem-16` (16 vCPU, 104 GB) | 200 GB | No | ~$0.50 |
| **Total per batch** | — | — | — | **~$5-9 (pilot), ~$10-20 (typical cohort)** |

### 8.8 Trigger Model

**Manual second submission** per batch on Terra:
1. Main pipeline (`pipeline.wdl`) completes Pegasus stage
2. Note the `hybrid_filtered_h5ad` URL from Pegasus task output
3. Submit `downstream.wdl` against that h5ad (or import as a separate
   workflow in the same workspace and select the `sample_set` row whose
   Pegasus output is ready)

### 8.9 Output Schema

| Output | Content |
|---|---|
| `annotated_h5ad` | Final h5ad with `obs.azimuth` (reconciled label), `obs.subclass` (derived), `obs.individualID` (from sample_to_individual map), re-clustered UMAP |
| `reconciled_csv` | Per-cell barcode -> final label (cells where AHBA + Ma-Sestan classes agree) |
| `inconsistent_csv` | Per-cell barcode -> both labels (cells where AHBA + Ma-Sestan disagreed; typically <10% drop rate) |
| `azimuth_umap_pngs` | UMAP visualizations colored by Azimuth labels (AHBA, Ma-Sestan) |

### 8.10 Code Refactoring (from existing `downstream/` scripts)

Scripts in `downstream/` are not WDL-friendly (hardcoded paths, Slurm
env vars, `data_prefix = "..."` literals). Refactored versions go in
`terra_wdl/`:

| Original (Slurm) | Refactored (Terra) | Key changes |
|---|---|---|
| `downstream/PEC_generate_matrices_hybrid.py` | (folded into `Pegasus-Pipeline.py`) | Only the filter sub-block survives; rest is duplicate of existing Pegasus stage |
| `downstream/Hybrid_Azimuth.R` | `terra_wdl/hybrid_azimuth_task.R` | CLI via `optparse`; drop `use_condaenv("test2")`; drop hardcoded `data_prefix` |
| `downstream/ahba_ma_map_v2.py` | `terra_wdl/reconcile_labels.py` | CLI args; drop `dataset_list` hardcoding; keep `ahba_map`/`ma_map` schema dicts as module constants |
| `downstream/finalize_h5ad.py` | `terra_wdl/finalize_h5ad_task.py` | Strip all `elif study=="CMC"` / `"UCLA-ASD"` / etc. branches; accept generic `--sample_to_individual_tsv` |

### 8.11 Implementation Status

| Item | Status |
|---|---|
| `terra_wdl/downstream.wdl` | Authored (3 tasks, 7 outputs) |
| `terra_wdl/hybrid_azimuth_task.R` | Authored (optparse CLI) |
| `terra_wdl/reconcile_labels.py` | Authored (CLI: `--ahba_csv`, `--ma_sestan_csv`, `--output_prefix`) |
| `terra_wdl/finalize_h5ad_task.py` | Authored; review-pass fix applied (only drop cells whose freshly-merged `azimuth` col is null, not any obs NaN) |
| `terra_wdl/Dockerfile.azimuth` | Authored (rocker/r-ver:4.3.2 + Seurat 5 + Azimuth) |
| `.github/workflows/build-azimuth-image.yml` | Authored |
| `.github/workflows/build-pegasus-image.yml` | `paths:` extended with `reconcile_labels.py`, `finalize_h5ad_task.py` |
| `terra_wdl/Dockerfile.pegasus` | `COPY` added for the two new downstream Python scripts |
| `Pegasus-Pipeline.py` filter step (0.3) | Implemented; writes `{batch}_Hybrid_filtered.h5ad` |
| `terra_wdl/pipeline.wdl` `hybrid_filtered_h5ad` output | Added via `glob()` |
| `terra_wdl/generate_sample_table.py --write_individual_map` | Implemented |
| `dataset/sample_to_individual_pilot.tsv` | Created (HCTBA, HCCLC) |
| `dataset/downstream_pilot_inputs.json` | Created (docker tags = `:latest` placeholder, h5ad URL placeholder) |
| **Commit + push** new files / changes | Pending (uncommitted in working tree) |
| GitHub Actions builds (Pegasus rebuild, Azimuth first build) | Pending — gated on commit/push |
| Pin SHA tags in pilot + downstream JSONs | Pending — gated on builds |
| Stage `.RDS` references to `gs://...references/azimuth/` | Pending (user `gsutil cp`) |
| Stage `sample_to_individual_pilot.tsv` to GCS | Pending (user `gsutil cp` to `gs://.../references/`) |
| Upload `pipeline.wdl` Snapshot 5 to Terra | Pending |
| Upload `downstream.wdl` as new Terra method | Pending |
| End-to-end pilot validation | Pending |

---

## 9. Open Follow-ups (Future Work)

- **Migrate Azimuth references** to a dedicated long-lived GCS bucket
  (outside the Terra workspace `fc-<UUID>` bucket) once cohort cadence
  is established. Avoids losing ~4 GB of reference data if workspace
  is recreated.
- **Promote downstream to Option A** (inline `call DownstreamAnnotation`
  at the end of `pipeline.wdl`) once the standalone workflow is
  validated against several real cohorts. Reduces operator burden of
  managing two submissions per batch.
- **Preemptible Azimuth variant** for non-time-sensitive batch
  processing of many cohorts — could cut ~70% off Azimuth compute cost
  with `preemptible: 1` and a non-preemptible fallback on second attempt.
- **Right-size Azimuth machine** after profiling actual RAM usage on
  the first few real cohorts — `n1-highmem-32` (208 GB) may suffice
  for smaller batches.
- **Unify Pegasus image SHAs** across full-pipeline and NeMo pilots
  once `sha-c5d5dbb` (or its post-filter-step successor) is validated.
- **`sample_to_individual.tsv` automation**: extend
  `generate_sample_table.py --write-individual-map` to handle batches
  that pool data from multiple metadata sources (currently only
  single-source `Prelim_Data.tsv` or `nemo_metadata.tsv`).
