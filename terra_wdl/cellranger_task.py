"""
cellranger_task.py

WDL-compatible CellRanger count runner.
Called by the CellRangerCount task in pipeline.wdl.

Replaces the run_cellranger() function from PEC_scRNA_pipeline_CellRanger.py.

Key differences from the original:
  - Accepts individual CLI arguments instead of a JSON config file.
  - Runs CellRanger in --jobmode=local (no Slurm).
  - No multiprocessing.Pool: WDL scatter handles per-sample parallelism.
  - Localizes the FASTQ directory from GCS to the local VM disk before running.
  - Writes outputs to structured paths that WDL captures as task outputs.
  - Prints expected_cells and low_quality flag to stdout for WDL read_int/read_boolean.

Usage (called from WDL command block):
    python /opt/pipeline/cellranger_task.py \
        --fastq_dir        gs://bucket/path/sample_folder/ \
        --sample_tag       MS1023HH_010720 \
        --transcriptome    gs://bucket/references/refdata-gex-GRCh38-2024-A \
        --chemistry        v3 \
        --include_introns  true \
        --numproc          16 \
        --min_umi_check    1000 \
        --billing_project  terra-4ea165c8
"""

import argparse
import csv
import json
import os
import subprocess
import sys


def localize_fastq_dir(gcs_path: str, local_dir: str,
                       billing_project: str = "") -> str:
    """
    Copy a GCS FASTQ directory to the local VM disk using gsutil.
    Returns the local path to the downloaded folder.
    CellRanger requires FASTQs to be on local disk.

    If billing_project is set, passes -u <project> for requester-pays buckets.
    """
    os.makedirs(local_dir, exist_ok=True)
    print(f"[cellranger_task] Localizing FASTQs from {gcs_path} → {local_dir}")
    cmd = ["gsutil"]
    if billing_project:
        cmd += ["-u", billing_project]
    cmd += ["-m", "cp", "-r", gcs_path.rstrip("/"), local_dir]
    result = subprocess.run(cmd, check=True)
    # gsutil cp -r copies the folder itself into local_dir
    folder_name = gcs_path.rstrip("/").split("/")[-1]
    local_fastq_path = os.path.join(local_dir, folder_name)
    print(f"[cellranger_task] FASTQs localized to: {local_fastq_path}")
    return local_fastq_path


def localize_reference(gcs_path: str, local_dir: str = "./reference") -> str:
    """
    Copy CellRanger reference directory from GCS to local disk.
    Returns the local path to the reference directory.
    Skips download if gcs_path is already a local path (not gs://).
    """
    if not gcs_path.startswith("gs://"):
        print(f"[cellranger_task] Transcriptome path is already local: {gcs_path}")
        return gcs_path
    os.makedirs(local_dir, exist_ok=True)
    print(f"[cellranger_task] Localizing reference from {gcs_path} -> {local_dir}")
    cmd = ["gsutil", "-m", "cp", "-r", gcs_path.rstrip("/"), local_dir]
    subprocess.run(cmd, check=True)
    folder_name = gcs_path.rstrip("/").split("/")[-1]
    local_ref_path = os.path.join(local_dir, folder_name)
    print(f"[cellranger_task] Reference localized to: {local_ref_path}")
    return local_ref_path


def run_cellranger_count(args) -> str:
    """
    Run cellranger count and return the path to the output folder.
    """
    include_flag = "--include-introns" if args.include_introns.lower() == "true" else ""

    cmd = (
        f"cellranger count"
        f"  --id={args.sample_tag}"
        f"  --fastqs={args.local_fastq_dir}"
        f"  --sample={args.sample_tag}"
        f"  --transcriptome={args.transcriptome}"
        f"  --chemistry={args.chemistry}"
        f"  --localcores={args.numproc}"
        f"  --localmem={args.localmem}"
        f"  --jobmode=local"
        f"  {include_flag}"
    )
    print(f"[cellranger_task] Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    output_dir = os.path.join(os.getcwd(), args.sample_tag)
    print(f"[cellranger_task] CellRanger count complete. Output dir: {output_dir}")
    return output_dir


def parse_metrics(metrics_csv_path: str) -> dict:
    """
    Parse CellRanger metrics_summary.csv and return a dict of key metrics.
    Handles comma-formatted numbers (e.g. "12,345" → 12345).
    """
    with open(metrics_csv_path, newline="") as f:
        reader = csv.DictReader(f)
        row = next(reader)

    def clean_int(val: str) -> int:
        return int(str(val).replace(",", "").strip())

    metrics = {}
    # Estimated Number of Cells is the first column
    first_key = list(row.keys())[0]
    metrics["estimated_cells"] = clean_int(row[first_key])

    if "Median UMI Counts per Cell" in row:
        metrics["median_umis_per_cell"] = clean_int(row["Median UMI Counts per Cell"])
    else:
        # Fallback key name variants across CellRanger versions
        for key in row:
            if "median" in key.lower() and "umi" in key.lower():
                metrics["median_umis_per_cell"] = clean_int(row[key])
                break
        else:
            metrics["median_umis_per_cell"] = 0

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="WDL-compatible CellRanger count runner for Terra."
    )
    parser.add_argument(
        "--fastq_dir", required=True,
        help="GCS path to the FASTQ folder (gs://bucket/path/sample_folder/)."
    )
    parser.add_argument(
        "--sample_tag", required=True,
        help="Sample name matching the FASTQ filename prefix (CellRanger --sample)."
    )
    parser.add_argument(
        "--transcriptome", required=True,
        help="GCS path or local path to the CellRanger reference transcriptome "
             "directory. If a gs:// path is given, it will be downloaded first."
    )
    parser.add_argument(
        "--chemistry", required=False, default="auto",
        help="10x chemistry version: v2, v3, or auto (default: auto)."
    )
    parser.add_argument(
        "--include_introns", required=False, default="true",
        help="Include intronic reads for snRNA-seq: true or false (default: true)."
    )
    parser.add_argument(
        "--numproc", required=False, type=int, default=16,
        help="Number of local cores for CellRanger (default: 16)."
    )
    parser.add_argument(
        "--localmem", required=False, type=int, default=60,
        help="Memory in GB available to CellRanger (default: 60)."
    )
    parser.add_argument(
        "--min_umi_check", required=False, type=int, default=1000,
        help="Median UMI threshold below which sample is flagged low quality (default: 1000)."
    )
    parser.add_argument(
        "--billing_project", required=False, default="",
        help="GCP billing project for requester-pays buckets (e.g. terra-4ea165c8). "
             "If empty, gsutil runs without -u flag."
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Localize FASTQ directory from GCS to local disk
    # ------------------------------------------------------------------
    local_fastqs_root = "./fastqs"
    args.local_fastq_dir = localize_fastq_dir(
        args.fastq_dir, local_fastqs_root,
        billing_project=args.billing_project
    )

    # ------------------------------------------------------------------
    # 1b. Localize reference transcriptome from GCS to local disk
    # ------------------------------------------------------------------
    args.transcriptome = localize_reference(args.transcriptome)

    # ------------------------------------------------------------------
    # 2. Run CellRanger count
    # ------------------------------------------------------------------
    output_dir = run_cellranger_count(args)

    # ------------------------------------------------------------------
    # 3. Validate outputs exist
    # ------------------------------------------------------------------
    outs_dir             = os.path.join(output_dir, "outs")
    raw_h5               = os.path.join(outs_dir, "raw_feature_bc_matrix.h5")
    filtered_h5          = os.path.join(outs_dir, "filtered_feature_bc_matrix.h5")
    metrics_csv          = os.path.join(outs_dir, "metrics_summary.csv")
    molecule_info_h5     = os.path.join(outs_dir, "molecule_info.h5")

    for path in [raw_h5, filtered_h5, metrics_csv, molecule_info_h5]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: Expected CellRanger output not found: {path}")

    # ------------------------------------------------------------------
    # 4. Parse QC metrics
    # ------------------------------------------------------------------
    metrics = parse_metrics(metrics_csv)
    estimated_cells    = metrics["estimated_cells"]
    median_umis        = metrics["median_umis_per_cell"]
    is_low_quality     = median_umis < args.min_umi_check

    print(f"[cellranger_task] Estimated cells:      {estimated_cells}")
    print(f"[cellranger_task] Median UMIs per cell: {median_umis}")
    print(f"[cellranger_task] Low quality flag:     {is_low_quality}")

    if is_low_quality:
        print(
            f"[cellranger_task] WARNING: Sample {args.sample_tag} has only "
            f"{median_umis} median UMIs per cell (threshold: {args.min_umi_check}). "
            f"Flagged as potentially low quality."
        )

    # ------------------------------------------------------------------
    # 5. Write output files consumed by WDL read_* functions
    #    WDL captures these via output { } declarations.
    # ------------------------------------------------------------------
    # Absolute paths written to text files for WDL read_string()
    with open("raw_h5_path.txt", "w") as f:
        f.write(raw_h5)
    with open("filtered_h5_path.txt", "w") as f:
        f.write(filtered_h5)
    with open("molecule_info_path.txt", "w") as f:
        f.write(molecule_info_h5)
    with open("metrics_csv_path.txt", "w") as f:
        f.write(metrics_csv)

    # Integer outputs for WDL read_int()
    with open("estimated_cells.txt", "w") as f:
        f.write(str(estimated_cells))

    # Boolean output for WDL read_boolean()
    with open("is_low_quality.txt", "w") as f:
        f.write(str(is_low_quality).lower())  # "true" or "false"

    # Full metrics JSON for downstream reference
    metrics["sample_tag"]   = args.sample_tag
    metrics["fastq_dir"]    = args.fastq_dir
    metrics["chemistry"]    = args.chemistry
    metrics["include_introns"] = args.include_introns
    with open("cellranger_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[cellranger_task] Done. Outputs:")
    print(f"  raw_h5:           {raw_h5}")
    print(f"  filtered_h5:      {filtered_h5}")
    print(f"  molecule_info_h5: {molecule_info_h5}")
    print(f"  metrics_csv:      {metrics_csv}")
    print(f"  estimated_cells:  {estimated_cells}")
    print(f"  is_low_quality:   {is_low_quality}")


if __name__ == "__main__":
    main()
