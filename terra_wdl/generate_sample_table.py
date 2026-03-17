"""
generate_sample_table.py

Transforms an exported Terra Prelim_Data TSV and a NeMo metadata TSV
into two TSV files ready for direct import into Terra:
  - {output_prefix}_sample_table.tsv     → import as "sample" entity
  - {output_prefix}_sample_set_table.tsv → import as "sample_set" entity

All fields are derived automatically from the two input files.
No manual editing of the output TSVs is required.

Usage:
    python generate_sample_table.py \
        --prelim_data   Prelim_Data.tsv \
        --nemo_metadata nemo_metadata.tsv \
        --batch_key     donor_region \
        --output_dir    ./output/ \
        --output_prefix my_study

Batch key options:
    donor_region  Group by subject_id + sample_anatomical_region (default)
    donor         Group by subject_id only
    project       Group by project_id only
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_fastq_dir(url: str) -> str:
    """
    Return the GCS directory path (with trailing slash) from a full GCS file URL.
    e.g. gs://bucket/path/sample_folder/sample_R1_001.fastq.gz
         → gs://bucket/path/sample_folder/
    """
    return url.rsplit("/", 1)[0] + "/"


def parse_sample_tag(filename: str) -> str:
    """
    Extract the CellRanger --sample tag from a FASTQ filename.
    Strips the Illumina lane/read suffix: _S{n}_L{n}_R{1,2}_{index}.fastq.gz

    e.g. MS1023HH_010720_S5_L004_R1_001.fastq.gz → MS1023HH_010720
    Falls back to stripping at _R1 / _R2 if the standard pattern is not found.
    """
    # Standard Illumina pattern
    match = re.match(r"^(.+?)_S\d+_L\d+_R[12]_\d+\.fastq\.gz$", filename)
    if match:
        return match.group(1)
    # Fallback: strip at _R1 or _R2
    for suffix in ("_R1", "_R2"):
        idx = filename.find(suffix)
        if idx != -1:
            return filename[:idx]
    # Last resort: strip extension only
    return filename.replace(".fastq.gz", "").replace(".fastq", "")


def parse_chemistry(technique: str) -> str:
    """
    Derive 10x chemistry version from the NeMo sample_technique string.
    e.g. "10x chromium 3' v3 sequencing" → "v3"
         "10x chromium 3' v3.1 sequencing" → "v3"
         "10x chromium 3' v2 sequencing" → "v2"
    Returns "v3" as a safe default if the version cannot be parsed.
    """
    if pd.isna(technique) or technique is None:
        return "v3"
    match = re.search(r"v(\d+)", str(technique))
    return f"v{match.group(1)}" if match else "v3"


def parse_include_introns(specimen_type: str) -> str:
    """
    Return "true" if specimen type is nuclei (snRNA-seq requires --include-introns),
    "false" for cells (scRNA-seq).
    """
    if pd.isna(specimen_type) or specimen_type is None:
        return "false"
    return "true" if str(specimen_type).strip().lower() == "nuclei" else "false"


def parse_hashing(modality: str, hashing_keywords: list) -> str:
    """
    Return "true" if the NeMo sample_modality matches any hashing keyword.
    """
    if pd.isna(modality) or modality is None:
        return "false"
    modality_lower = str(modality).strip().lower()
    return "true" if any(kw.lower() in modality_lower for kw in hashing_keywords) else "false"


# ---------------------------------------------------------------------------
# Secondary join: get anatomical region from smp-level metadata rows
# ---------------------------------------------------------------------------

def build_smp_region_map(smp_df: pd.DataFrame) -> dict:
    """
    Build a mapping: subject_id → list of (sample_name, anatomical_region)
    from sample-level (nemo:smp-*) metadata rows.
    """
    region_map = {}
    for _, row in smp_df.iterrows():
        sid = row.get("subject_id", None)
        if pd.isna(sid) or sid is None:
            continue
        region = row.get("sample_anatomical_region", None)
        name = row.get("sample_name", "")
        if sid not in region_map:
            region_map[sid] = []
        region_map[sid].append((str(name) if not pd.isna(name) else "", region))
    return region_map


def lookup_anatomical_region(subject_id: str, lib_sample_name: str,
                              region_map: dict) -> str:
    """
    Look up the anatomical region for a library row using the smp-level map.

    Priority:
      1. Single smp row for this subject → use its region directly.
      2. Multiple smp rows → try prefix matching on sample_name.
      3. Still ambiguous → use the first match and emit a warning.
      4. Subject not found in smp rows → return "unknown".
    """
    entries = region_map.get(subject_id, [])
    if not entries:
        return "unknown"
    if len(entries) == 1:
        region = entries[0][1]
        return str(region).strip().lower().replace(" ", "_") if not pd.isna(region) else "unknown"

    # Multiple entries: try prefix matching
    lib_name = str(lib_sample_name) if not pd.isna(lib_sample_name) else ""
    for smp_name, region in entries:
        if lib_name.startswith(smp_name) or smp_name.startswith(lib_name):
            return str(region).strip().lower().replace(" ", "_") if not pd.isna(region) else "unknown"

    # Ambiguous: warn and use first entry
    warnings.warn(
        f"Multiple anatomical regions found for subject {subject_id} and no "
        f"unambiguous prefix match for library sample_name='{lib_sample_name}'. "
        f"Using first match: '{entries[0][1]}'. "
        f"Consider using --batch_key donor to avoid ambiguity.",
        UserWarning
    )
    region = entries[0][1]
    return str(region).strip().lower().replace(" ", "_") if not pd.isna(region) else "unknown"


# ---------------------------------------------------------------------------
# Batch number assignment
# ---------------------------------------------------------------------------

def assign_batch_numbers(df: pd.DataFrame, batch_key: str,
                          region_map: dict, label_col: str) -> pd.DataFrame:
    """
    Add batch_num (int) and sample_set_id (human-readable string) columns to df.

    batch_key values:
      donor_region  → (subject_id, anatomical_region) compound key
      donor         → subject_id only
      project       → project_id only
    """
    if batch_key == "donor_region":
        # Ensure anatomical_region column is populated via secondary join
        if "anatomical_region" not in df.columns:
            df["anatomical_region"] = df.apply(
                lambda r: lookup_anatomical_region(
                    r["subject_id"], r.get("sample_name", ""), region_map
                ),
                axis=1
            )
        group_cols = ["subject_id", "anatomical_region"]
        label_fn = lambda r: (
            f"{r[label_col]}_{r['anatomical_region']}"
            .replace(" ", "_").replace("/", "_")
        )

    elif batch_key == "donor":
        group_cols = ["subject_id"]
        label_fn = lambda r: str(r[label_col]).replace(" ", "_")

    elif batch_key == "project":
        group_cols = ["project_id"]
        label_fn = lambda r: str(r["project_id"]).replace(":", "_")

    else:
        raise ValueError(
            f"Unknown --batch_key '{batch_key}'. "
            "Choose from: donor_region, donor, project"
        )

    # Create a sorted mapping of group key → batch number
    group_keys = (
        df[group_cols]
        .drop_duplicates()
        .sort_values(group_cols)
        .reset_index(drop=True)
    )
    group_keys["batch_num"] = group_keys.index + 1

    df = df.merge(group_keys, on=group_cols, how="left")

    # Human-readable sample_set_id
    df["sample_set_id"] = df.apply(label_fn, axis=1)

    return df


# ---------------------------------------------------------------------------
# Hashing FASTQ pairing
# ---------------------------------------------------------------------------

def pair_hashing_fastqs(df: pd.DataFrame) -> pd.DataFrame:
    """
    For samples where hashing=true, attempt to automatically populate
    hashing_fastq_r1 and hashing_fastq_r2 by finding the corresponding
    hashing library rows in the same dataframe (matched by subject_id).

    Hashing libraries are identified by their sample_modality being in the
    hashing keywords (already parsed into the hashing column = "true").
    """
    # Rows where hashing library itself is the source of HTO FASTQs
    hashing_libs = df[df["hashing"] == "true"].copy()
    rna_libs = df[df["hashing"] == "false"].copy()

    if hashing_libs.empty:
        df["hashing_fastq_r1"] = ""
        df["hashing_fastq_r2"] = ""
        return df

    # Build subject → hashing FASTQ dir map
    hash_fastq_map = {}
    for _, row in hashing_libs.iterrows():
        subj = row["subject_id"]
        fdir = row["fastq_dir"]
        if subj not in hash_fastq_map:
            hash_fastq_map[subj] = fdir

    def get_hash_r1(row):
        if row["hashing"] != "true":
            return ""
        fdir = hash_fastq_map.get(row["subject_id"], "")
        # Reconstruct R1 path: fastq_dir + sample_tag + _R1 suffix placeholder
        # Terra/CellRanger will resolve the actual file; store dir-level path
        return fdir

    def get_hash_r2(row):
        if row["hashing"] != "true":
            return ""
        return hash_fastq_map.get(row["subject_id"], "")

    df["hashing_fastq_r1"] = df.apply(get_hash_r1, axis=1)
    df["hashing_fastq_r2"] = df.apply(get_hash_r2, axis=1)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate Terra-ready sample and sample_set TSV tables from "
            "a Terra Prelim_Data export and a NeMo metadata file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--prelim_data", required=True,
        help="Path to the Prelim_Data TSV exported from Terra."
    )
    parser.add_argument(
        "--nemo_metadata", required=True,
        help="Path to the NeMo metadata TSV downloaded from the NeMo portal."
    )
    parser.add_argument(
        "--batch_key", default="donor_region",
        choices=["donor_region", "donor", "project"],
        help=(
            "Strategy for grouping samples into batches. "
            "donor_region: subject_id + anatomical_region (default). "
            "donor: subject_id only. "
            "project: project_id only."
        )
    )
    parser.add_argument(
        "--output_dir", default=".",
        help="Directory to write output TSV files (default: current directory)."
    )
    parser.add_argument(
        "--output_prefix", default="sample",
        help="Prefix for output filenames (default: 'sample')."
    )
    parser.add_argument(
        "--batch_label_col", default="subject_name",
        help=(
            "Metadata column to use for human-readable sample_set_id labels "
            "(default: subject_name)."
        )
    )
    parser.add_argument(
        "--hashing_keywords",
        default="hashing,HTO,antibody_capture,multiseq,citeseq",
        help=(
            "Comma-separated keywords in sample_modality that indicate a "
            "hashing library (default: hashing,HTO,antibody_capture,multiseq,citeseq)."
        )
    )
    parser.add_argument(
        "--min_umi_check", type=int, default=1000,
        help=(
            "Flag samples with estimated cell counts below this UMI threshold "
            "as potentially low quality in the output table (default: 1000). "
            "This mirrors the quality check in CellRanger.py and adds an "
            "'expected_low_quality' column to the sample table for review."
        )
    )
    args = parser.parse_args()

    hashing_keywords = [k.strip() for k in args.hashing_keywords.split(",")]

    # ------------------------------------------------------------------
    # 1. Load input files
    # ------------------------------------------------------------------
    print(f"[1/7] Loading Prelim_Data from: {args.prelim_data}")
    prelim_df = pd.read_csv(args.prelim_data, sep="\t", dtype=str)
    # Normalise column name — Terra exports with 'entity:Prelim_Data_id'
    prelim_df.columns = [c.replace("entity:", "") for c in prelim_df.columns]
    required_prelim = {"Prelim_Data_id", "sample_id", "urls"}
    missing = required_prelim - set(prelim_df.columns)
    if missing:
        sys.exit(
            f"ERROR: Prelim_Data TSV is missing required columns: {missing}\n"
            f"Found columns: {list(prelim_df.columns)}"
        )

    print(f"[1/7] Loading NeMo metadata from: {args.nemo_metadata}")
    meta_df = pd.read_csv(args.nemo_metadata, sep="\t", dtype=str)
    required_meta = {"sample_id", "subject_id", "subject_name",
                     "sample_technique", "sample_specimentype",
                     "sample_modality", "project_id"}
    missing = required_meta - set(meta_df.columns)
    if missing:
        sys.exit(
            f"ERROR: NeMo metadata TSV is missing required columns: {missing}\n"
            f"Found columns: {list(meta_df.columns)}"
        )

    # ------------------------------------------------------------------
    # 2. Split metadata into library-level and sample-level rows
    # ------------------------------------------------------------------
    print("[2/7] Separating library-level and sample-level metadata rows ...")
    lib_df = meta_df[meta_df["sample_id"].str.startswith("nemo:lib-")].copy()
    smp_df = meta_df[meta_df["sample_id"].str.startswith("nemo:smp-")].copy()

    print(f"      Library-level rows (nemo:lib-*): {len(lib_df)}")
    print(f"      Sample-level rows  (nemo:smp-*): {len(smp_df)}")

    if lib_df.empty:
        sys.exit(
            "ERROR: No nemo:lib-* rows found in the metadata file. "
            "Ensure the metadata TSV contains library-level entries."
        )

    # ------------------------------------------------------------------
    # 3. Join Prelim_Data ← lib-level metadata on sample_id
    # ------------------------------------------------------------------
    print("[3/7] Joining Prelim_Data to library-level metadata ...")
    merged = prelim_df.merge(lib_df, on="sample_id", how="left", suffixes=("", "_meta"))

    unmatched = merged["subject_id"].isna().sum()
    if unmatched > 0:
        warnings.warn(
            f"{unmatched} rows in Prelim_Data had no matching library entry "
            f"in the NeMo metadata (sample_id not found in nemo:lib-* rows). "
            f"These rows will be skipped.",
            UserWarning
        )
    merged = merged[merged["subject_id"].notna()].copy()

    # ------------------------------------------------------------------
    # 4. Derive per-row fields
    # ------------------------------------------------------------------
    print("[4/7] Deriving fastq_dir, sample_tag, chemistry, include_introns, hashing ...")
    merged["fastq_dir"]       = merged["urls"].apply(parse_fastq_dir)
    merged["sample_tag"]      = merged["Prelim_Data_id"].apply(parse_sample_tag)
    merged["chemistry"]       = merged["sample_technique"].apply(parse_chemistry)
    merged["include_introns"] = merged["sample_specimentype"].apply(parse_include_introns)
    merged["hashing"]         = merged["sample_modality"].apply(
        lambda m: parse_hashing(m, hashing_keywords)
    )

    # ------------------------------------------------------------------
    # 5. Deduplicate by sample_tag (collapse R1 + R2 rows → one per sample)
    #    Keep the first occurrence; fastq_dir is the same for both reads.
    # ------------------------------------------------------------------
    print("[5/7] Deduplicating R1/R2 rows → one row per sample ...")
    before = len(merged)
    merged = merged.drop_duplicates(subset="sample_tag", keep="first").copy()
    after = len(merged)
    print(f"      Collapsed {before} file rows → {after} sample rows.")

    # ------------------------------------------------------------------
    # 6. Assign batch numbers
    # ------------------------------------------------------------------
    print(f"[6/7] Assigning batch numbers using batch_key='{args.batch_key}' ...")
    region_map = build_smp_region_map(smp_df)
    merged = assign_batch_numbers(merged, args.batch_key, region_map, args.batch_label_col)
    print(f"      Found {merged['batch_num'].nunique()} unique batch(es).")

    # ------------------------------------------------------------------
    # 6b. Hashing FASTQ pairing
    # ------------------------------------------------------------------
    merged = pair_hashing_fastqs(merged)

    # ------------------------------------------------------------------
    # 6c. Low-quality flag
    # ------------------------------------------------------------------
    # We don't have UMI counts at this stage (pre-CellRanger), so this column
    # is a placeholder. It will be populated post-CellRanger in Terra.
    merged["expected_low_quality"] = "false"

    # ------------------------------------------------------------------
    # 7. Build output tables
    # ------------------------------------------------------------------
    print("[7/7] Writing output TSV files ...")
    os.makedirs(args.output_dir, exist_ok=True)

    # --- sample_table ---
    sample_cols = [
        "sample_tag",           # used as entity:sample_id
        "fastq_dir",
        "sample_tag",           # CellRanger --sample= argument (same value)
        "chemistry",
        "include_introns",
        "batch_num",
        "sample_set_id",
        "hashing",
        "hashing_fastq_r1",
        "hashing_fastq_r2",
        "expected_low_quality",
        "subject_id",
        "subject_name",
        "project_id",
        "sample_id",            # NeMo library ID (nemo:lib-*)
    ]
    # Add anatomical_region if it was computed
    if "anatomical_region" in merged.columns:
        sample_cols.append("anatomical_region")

    # Remove duplicates from column list (sample_tag appears twice above intentionally
    # as both entity id and the sample_tag value column — keep distinct names)
    sample_out = merged.copy()
    sample_out.rename(columns={"sample_id": "nemo_library_id"}, inplace=True)

    final_sample_cols = [
        "sample_tag", "fastq_dir", "chemistry", "include_introns",
        "batch_num", "sample_set_id", "hashing",
        "hashing_fastq_r1", "hashing_fastq_r2",
        "expected_low_quality", "subject_id", "subject_name",
        "project_id", "nemo_library_id",
    ]
    if "anatomical_region" in sample_out.columns:
        final_sample_cols.append("anatomical_region")

    sample_out = sample_out[final_sample_cols].copy()
    sample_out.insert(0, "entity:sample_id", sample_out["sample_tag"])

    sample_path = os.path.join(args.output_dir, f"{args.output_prefix}_sample_table.tsv")
    sample_out.to_csv(sample_path, sep="\t", index=False)
    print(f"      Wrote sample table ({len(sample_out)} rows): {sample_path}")

    # --- sample_set_table ---
    set_rows = []
    for set_id, group in sample_out.groupby("sample_set_id"):
        member_ids = ",".join(group["entity:sample_id"].tolist())
        batch_num  = group["batch_num"].iloc[0]
        set_rows.append({
            "entity:sample_set_id": set_id,
            "samples":              member_ids,
            "batch_num":            batch_num,
        })
    set_df = pd.DataFrame(set_rows).sort_values("batch_num").reset_index(drop=True)
    set_path = os.path.join(args.output_dir, f"{args.output_prefix}_sample_set_table.tsv")
    set_df.to_csv(set_path, sep="\t", index=False)
    print(f"      Wrote sample_set table ({len(set_df)} batch(es)): {set_path}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n--- Summary ---")
    print(f"  Total samples:        {len(sample_out)}")
    print(f"  Total batches:        {len(set_df)}")
    print(f"  Hashing samples:      {(sample_out['hashing'] == 'true').sum()}")
    print(f"  snRNA (nuclei):       {(sample_out['include_introns'] == 'true').sum()}")
    print(f"  scRNA (cells):        {(sample_out['include_introns'] == 'false').sum()}")
    chemistry_counts = sample_out["chemistry"].value_counts().to_dict()
    for chem, count in chemistry_counts.items():
        print(f"  Chemistry {chem}:         {count}")
    print(f"\nNext steps:")
    print(f"  1. Import '{sample_path}' into Terra as a 'sample' entity table.")
    print(f"  2. Import '{set_path}' into Terra as a 'sample_set' entity table.")
    print(f"  3. Add workspace attributes: transcriptome_ref, mito_file,")
    print(f"     ahba_markers_json, hybrid_markers_json.")
    print(f"  4. Submit pipeline.wdl on the desired sample_set(s).")


if __name__ == "__main__":
    main()
