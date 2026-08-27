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
import re
import subprocess
import sys


def localize_fastq_dirs(gcs_paths, local_root: str,
                        billing_project: str = "",
                        sample_prefix: str = ""):
    """
    Copy one or more GCS FASTQ directories to the local VM disk using gsutil.
    Each source directory is copied into its own local subfolder to avoid
    filename collisions between lanes stored in separate directories.
    Returns a list of local directory paths (one per source).

    CellRanger requires FASTQs to be on local disk. `--fastqs` accepts a
    comma-separated list of directories, so multi-lane samples whose lanes
    live in different GCS folders are all provided.

    If `sample_prefix` is given, only files matching `{sample_prefix}_S*` are
    copied. This is REQUIRED for flat shared buckets (e.g. gs://.../raw/) that
    hold many samples — and ATAC libraries — in one directory: copying `dir/*`
    there would download the whole bucket per sample and drag in the wrong
    (e.g. *_ATAC_*) files. When empty, the whole directory is copied (correct
    for per-sample directories).

    If billing_project is set, passes -u <project> for requester-pays buckets.
    """
    os.makedirs(local_root, exist_ok=True)
    local_paths = []
    for idx, gcs_path in enumerate(gcs_paths):
        gcs_path = gcs_path.strip()
        if not gcs_path:
            continue
        # Unique local subfolder per source dir (index prefix guards against
        # two source dirs sharing the same basename).
        subdir = os.path.join(local_root, f"d{idx:03d}")
        os.makedirs(subdir, exist_ok=True)
        # Restrict to this sample's files when a prefix is known, else grab all.
        if sample_prefix:
            src = gcs_path.rstrip("/") + f"/{sample_prefix}_S*"
        else:
            src = gcs_path.rstrip("/") + "/*"
        print(f"[cellranger_task] Localizing FASTQs from {src} → {subdir}")
        cmd = ["gsutil"]
        if billing_project:
            cmd += ["-u", billing_project]
        # Copy the matching files into the local subdir.
        cmd += ["-m", "cp", "-r", src, subdir]
        subprocess.run(cmd, check=True)
        local_paths.append(subdir)
    if not local_paths:
        sys.exit("ERROR: no FASTQ directories were provided/localized.")
    print(f"[cellranger_task] Localized {len(local_paths)} FASTQ dir(s).")
    return local_paths


def normalize_fastq_names(local_dirs, convention: str) -> None:
    """
    Rename localized FASTQ files in-place so CellRanger can parse them.

      standard : no change (files already `_S{n}_L{lane}_R{k}_001.fastq.gz`
                 or `_S{n}_R{k}_001.fastq.gz`).
      dotlane  : `..._L001.R1_001.fastq.gz` → `..._L001_R1_001.fastq.gz`
                 (replace the dot before the read token with an underscore).
      r3       : 4-file 10x dump `R1/R2/R3/I1` where (confirmed by read-length
                 inspection on GCS):
                   R1 = 28 bp  cell barcode + UMI   → CellRanger R1 (keep)
                   R2 = 10 bp  sample index         → drop (not needed)
                   R3 = 90 bp  cDNA                 → CellRanger R2 (rename)
                   I1 = 10 bp  sample index         → drop (not needed)
                 `cellranger count` on already-demultiplexed FASTQs needs only
                 the barcode (R1) and cDNA (R2) reads. The separator before the
                 read token may be '.' (Spudich) or '_' (raw); both handled.
    """
    if convention == "standard":
        return

    if convention == "dotlane":
        dot_read = re.compile(r"\.([RI][123])_(\d+\.fastq\.gz)$")
        n_renamed = 0
        for d in local_dirs:
            for fn in os.listdir(d):
                new = dot_read.sub(r"_\1_\2", fn)
                if new != fn:
                    os.rename(os.path.join(d, fn), os.path.join(d, new))
                    n_renamed += 1
        print(f"[cellranger_task] dotlane normalisation: renamed {n_renamed} file(s).")
        return

    if convention == "r3":
        # Read token preceded by '.' or '_'; group(1)=token, group(2)=suffix.
        read_re = re.compile(r"[._]([RI][123])_(\d+\.fastq\.gz)$")
        n_dropped = n_renamed = 0
        for d in local_dirs:
            # Pass 1: remove the index reads (original R2, I1/I2) FIRST so the
            # R3→R2 rename below cannot collide with an existing R2 file.
            for fn in list(os.listdir(d)):
                m = read_re.search(fn)
                if m and m.group(1) in ("R2", "I1", "I2"):
                    os.remove(os.path.join(d, fn))
                    n_dropped += 1
            # Pass 2: rename R3→R2 (cDNA) and normalise the R1 separator to '_'.
            for fn in list(os.listdir(d)):
                m = read_re.search(fn)
                if not m:
                    continue
                token = m.group(1)
                if token == "R3":
                    new = read_re.sub(r"_R2_\2", fn)
                elif token == "R1":
                    new = read_re.sub(r"_R1_\2", fn)
                else:
                    continue
                if new != fn:
                    os.rename(os.path.join(d, fn), os.path.join(d, new))
                    n_renamed += 1
        print(f"[cellranger_task] r3 remap: dropped {n_dropped} index file(s), "
              f"renamed {n_renamed} file(s) (R3→R2, R1 separator normalised).")
        return

    sys.exit(f"ERROR: unknown fastq_convention '{convention}'.")


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


def detect_fastq_sample_prefix(fastq_dirs) -> str:
    """
    Auto-detect the sample prefix from FASTQ filenames across one or more
    directories. CellRanger expects filenames like
    {SampleName}_S{N}_L{lane}_R{read}_001.fastq.gz (or without the lane token).
    Returns the {SampleName} portion from the first matching file.
    """
    pattern = re.compile(
        r'^(.+?)_S\d+(?:_L\d+)?_[RI][123]_001\.fastq\.gz$')
    for d in fastq_dirs:
        for fname in sorted(os.listdir(d)):
            m = pattern.match(fname)
            if m:
                prefix = m.group(1)
                print(f"[cellranger_task] Auto-detected FASTQ sample prefix: {prefix}")
                return prefix
    return ""


def run_cellranger_count(args) -> str:
    """
    Run cellranger count and return the path to the output folder.
    """
    include_flag = f"--include-introns {args.include_introns.lower()}"

    # Determine sample prefix: explicit arg > auto-detect > fallback to sample_tag
    sample_prefix = args.fastq_sample_prefix if args.fastq_sample_prefix else ""
    if not sample_prefix:
        sample_prefix = detect_fastq_sample_prefix(args.local_fastq_dirs)
    if not sample_prefix:
        sample_prefix = args.sample_tag

    fastqs_arg = ",".join(args.local_fastq_dirs)

    cmd = (
        f"cellranger count"
        f"  --id={args.sample_tag}"
        f"  --fastqs={fastqs_arg}"
        f"  --sample={sample_prefix}"
        f"  --transcriptome={args.transcriptome}"
        f"  --chemistry={args.chemistry}"
        f"  --localcores={args.numproc}"
        f"  --localmem={args.localmem}"
        f"  --jobmode=local"
        f"  --create-bam false"
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
        "--fastq_dirs", required=True,
        help="Comma-separated GCS path(s) to the FASTQ folder(s) for this "
             "sample. Multiple entries are used when a sample's lanes live in "
             "separate directories (gs://.../L001/,gs://.../L002/)."
    )
    parser.add_argument(
        "--fastq_convention", required=False, default="standard",
        choices=["standard", "dotlane", "r3"],
        help="FASTQ naming/read-layout convention, used to normalise filenames "
             "before CellRanger (default: standard)."
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
        "--fastq_sample_prefix", required=False, default="",
        help="FASTQ sample name prefix (the part before _S{N}_L{lane} in filenames). "
             "If empty, falls back to --sample_tag."
    )
    parser.add_argument(
        "--billing_project", required=False, default="",
        help="GCP billing project for requester-pays buckets (e.g. terra-4ea165c8). "
             "If empty, gsutil runs without -u flag."
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Localize FASTQ directory(ies) from GCS to local disk, then normalise
    #    filenames so CellRanger can parse them.
    #    A known sample prefix restricts the copy to this sample's files —
    #    essential for flat shared buckets (gs://.../raw/) that also hold ATAC
    #    and other subjects. Falls back to the whole directory when unknown.
    # ------------------------------------------------------------------
    gcs_dirs = [d for d in args.fastq_dirs.split(",") if d.strip()]
    local_fastqs_root = "./fastqs"
    args.local_fastq_dirs = localize_fastq_dirs(
        gcs_dirs, local_fastqs_root,
        billing_project=args.billing_project,
        sample_prefix=args.fastq_sample_prefix,
    )
    normalize_fastq_names(args.local_fastq_dirs, args.fastq_convention)

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
    # 5. Write small output files consumed by WDL read_int/read_boolean.
    #    File outputs (h5, csv) are captured by WDL via glob() patterns
    #    in the output { } block; no path-marker files needed.
    # ------------------------------------------------------------------
    # Integer outputs for WDL read_int()
    with open("estimated_cells.txt", "w") as f:
        f.write(str(estimated_cells))

    # Boolean output for WDL read_boolean()
    with open("is_low_quality.txt", "w") as f:
        f.write(str(is_low_quality).lower())  # "true" or "false"

    # Full metrics JSON for downstream reference
    metrics["sample_tag"]   = args.sample_tag
    metrics["fastq_dirs"]   = args.fastq_dirs
    metrics["fastq_convention"] = args.fastq_convention
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
