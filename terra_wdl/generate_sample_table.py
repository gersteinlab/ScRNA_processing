"""
generate_sample_table.py

Transforms the SCORCH deduped file/subject metadata manifest (or any
cohort subset of it) into two TSV files ready for direct import into Terra:
  - {output_prefix}_sample_table.tsv     → import as "sample" entity
  - {output_prefix}_sample_set_table.tsv → import as "sample_set" entity

This generator is **manifest-driven**. It reads the tab-separated manifest
whose library-type column is `technique` (e.g. "10xMultiome_RNAseq"), groups
the per-lane library records into CellRanger samples, and emits one sample
row per CellRanger run.

Key behaviours
--------------
* Chemistry is derived from `technique` via TECHNIQUE_CHEMISTRY. Multiome RNA
  (`10xMultiome_RNAseq`) is mapped to CellRanger chemistry `ARC-v1` so it can
  be processed by the standard `cellranger count` tool (Option A). Techniques
  not in the map (cell-hashing, ATAC, Flex, VDJ, nanopore, pacbio, ...) are
  excluded from the run with a per-technique summary.
* The manifest `library_id` is a per-*lane* record. Multiple lanes of the same
  10x sample share a FASTQ sample prefix (the token before `_S{N}`) but live in
  separate GCS directories. Rows are grouped by (subject_id, sample_prefix) into
  one CellRanger sample, and *all* lane directories are collected into a
  comma-separated `fastq_dirs` field (CellRanger `--fastqs` accepts a list).
* FASTQ filenames use different conventions across labs. Each sample is tagged
  with a `fastq_convention` (standard | dotlane | r3) that the CellRanger task
  uses to normalise filenames / read roles before running.

Usage
-----
    python generate_sample_table.py \
        --manifest    Data/scorch_pfc_cohort_subset_2026-08-26.tsv \
        --batch_key   donor \
        --output_dir  ./terra_tables \
        --output_prefix scorch_pfc
"""

import argparse
import os
import re
import sys
import warnings

import pandas as pd


# ---------------------------------------------------------------------------
# Technique → CellRanger chemistry map. Keys double as the include allow-list;
# any technique not present here is excluded from the run.
# ---------------------------------------------------------------------------
TECHNIQUE_CHEMISTRY = {
    "10xMultiome_RNAseq": "ARC-v1",   # Option A: standard CellRanger, Multiome chemistry
    "10x_v2":   "v2",
    "10x_v3":   "v3",
    "10x_v3.1": "v3",
    "10x_v4":   "auto",               # GEM-X; let CellRanger 8+ auto-detect
}

DEFAULT_ANATOMICAL_SITES = [
    "Brodmann (1909) area 9",
    "Brodmann (1909) area 10",
    "prefrontal cortex",
]

# Manifest columns the generator relies on.
REQUIRED_MANIFEST_COLS = {
    "library_id", "subject_id", "subject_name", "project_short_name",
    "anatomical_site", "species", "technique", "specimen_type",
    "data_type", "format", "url", "file_name",
}


# ---------------------------------------------------------------------------
# FASTQ filename parsing
# ---------------------------------------------------------------------------
#   standard : {prefix}_S{n}_L{lane}_R1_001.fastq.gz      (underscore before read)
#   dotlane  : {prefix}_S{n}_L{lane}.R1_001.fastq.gz      (dot before read)
#   nolane   : {prefix}_S{n}_R1_001.fastq.gz              (no lane token)
FASTQ_PATTERNS = [
    ("standard", re.compile(
        r"^(?P<prefix>.+?)_S\d+_L(?P<lane>\d+)_(?P<read>[RI][123])_\d+\.fastq\.gz$")),
    ("dotlane", re.compile(
        r"^(?P<prefix>.+?)_S\d+_L(?P<lane>\d+)\.(?P<read>[RI][123])_\d+\.fastq\.gz$")),
    ("nolane", re.compile(
        r"^(?P<prefix>.+?)_S\d+_(?P<read>[RI][123])_\d+\.fastq\.gz$")),
]


def parse_fastq_dir(url: str) -> str:
    """Return the GCS directory (with trailing slash) from a full file URL."""
    return url.rsplit("/", 1)[0] + "/"


def parse_fastq_meta(filename: str):
    """
    Parse a 10x FASTQ filename into its sample prefix, read token, and naming
    style. Returns a dict {prefix, read, style, lane} or None if unrecognised.
    """
    for style, pat in FASTQ_PATTERNS:
        m = pat.match(filename)
        if m:
            g = m.groupdict()
            return {
                "prefix": g["prefix"],
                "read": g["read"],
                "style": style,
                "lane": g.get("lane"),
            }
    return None


def classify_convention(styles: set, reads: set) -> str:
    """
    Decide the normalisation convention for a CellRanger sample.

      r3       : an R3 read is present → non-standard read layout, the
                 CellRanger task must remap read roles before running.
      dotlane  : filenames use `.R1` (dot) → task must rename to `_R1`.
      standard : usable by CellRanger as-is (covers standard + nolane naming).
    """
    if "R3" in reads:
        return "r3"
    if "dotlane" in styles:
        return "dotlane"
    return "standard"


def parse_include_introns(specimen_type: str) -> str:
    """'true' for nuclei (snRNA-seq needs --include-introns), else 'false'."""
    if pd.isna(specimen_type) or specimen_type is None:
        return "false"
    return "true" if str(specimen_type).strip().lower() == "nuclei" else "false"


# ---------------------------------------------------------------------------
# Manifest loading + filtering
# ---------------------------------------------------------------------------

def load_manifest(path: str, species: str, sites: list) -> pd.DataFrame:
    """
    Load the manifest, filter to the requested cohort and to demultiplexed
    FASTQ rows, and return the surviving rows with parsed FASTQ metadata.
    """
    print(f"[1/5] Loading manifest: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)

    missing = REQUIRED_MANIFEST_COLS - set(df.columns)
    if missing:
        sys.exit(
            f"ERROR: manifest is missing required columns: {sorted(missing)}\n"
            f"Found columns: {list(df.columns)}"
        )

    n0 = len(df)

    # The manifest stores subject_id as a float-like string (e.g. "3647.0").
    # Strip the trailing ".0" so batch entity IDs / labels are clean integers.
    df["subject_id"] = df["subject_id"].str.replace(r"\.0$", "", regex=True)

    # Cohort filters
    if species and species != "*":
        df = df[df["species"] == species]
    if sites:
        df = df[df["anatomical_site"].isin(sites)]

    print(f"      Cohort filter (species={species!r}, sites={len(sites)}): "
          f"{n0} → {len(df)} rows")

    # Keep only demultiplexed FASTQ files (the runnable CellRanger input).
    before_fastq = len(df)
    df = df[(df["data_type"] == "demultiplexed_fastq") & (df["format"] == "fastq")]
    print(f"      FASTQ rows only: {before_fastq} → {len(df)} rows")

    if df.empty:
        sys.exit("ERROR: no demultiplexed FASTQ rows remain after filtering.")

    # Derive chemistry from technique; drop excluded techniques.
    df = df.copy()
    df["chemistry"] = df["technique"].map(TECHNIQUE_CHEMISTRY)
    excluded = df[df["chemistry"].isna()]
    if not excluded.empty:
        print("      Excluding techniques not in the run allow-list:")
        for tech, n in excluded["technique"].value_counts().items():
            print(f"        - {tech}: {n} file rows dropped")
    df = df[df["chemistry"].notna()].copy()

    if df.empty:
        sys.exit("ERROR: no rows left after technique allow-list filtering.")

    # Parse FASTQ filenames. NOTE: the manifest `file_name` column is unreliable
    # (~16% of rows drop infixes like `_HIV_` / `_HIVOUD_` that the real object
    # carries), so we parse the FASTQ metadata from the `url` basename, which is
    # the source of truth for what actually exists on GCS.
    df["url_base"] = df["url"].str.rsplit("/", n=1).str[-1]
    meta = df["url_base"].apply(parse_fastq_meta)
    unparsed = df[meta.isna()]
    if not unparsed.empty:
        warnings.warn(
            f"{len(unparsed)} FASTQ file(s) had unrecognised names and were "
            f"skipped. Example: {unparsed['url_base'].iloc[0]}",
            UserWarning,
        )
    df = df[meta.notna()].copy()
    meta = meta[meta.notna()]

    df["fastq_prefix"] = meta.apply(lambda d: d["prefix"])
    df["fastq_read"]   = meta.apply(lambda d: d["read"])
    df["fastq_style"]  = meta.apply(lambda d: d["style"])
    df["fastq_dir"]    = df["url"].apply(parse_fastq_dir)
    df["include_introns"] = df["specimen_type"].apply(parse_include_introns)

    return df


# ---------------------------------------------------------------------------
# Group per-lane FASTQ rows into CellRanger samples
# ---------------------------------------------------------------------------

def group_into_samples(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse per-lane FASTQ file rows into one row per CellRanger sample,
    grouped by (subject_id, fastq_prefix). Collect all lane directories into a
    comma-separated `fastq_dirs` field and classify the naming convention.
    """
    print("[2/5] Grouping per-lane FASTQ rows into CellRanger samples ...")
    records = []
    for (subject_id, prefix), grp in df.groupby(["subject_id", "fastq_prefix"], sort=True):
        dirs = sorted(grp["fastq_dir"].unique())
        styles = set(grp["fastq_style"])
        reads = set(grp["fastq_read"])
        convention = classify_convention(styles, reads)

        def first(col):
            return grp[col].iloc[0]

        # Chemistry / include_introns should be consistent within a sample;
        # warn if a sample mixes techniques (should not happen in practice).
        if grp["chemistry"].nunique() > 1:
            warnings.warn(
                f"Sample {prefix} (subject {subject_id}) mixes chemistries "
                f"{sorted(grp['chemistry'].unique())}; using the first.",
                UserWarning,
            )

        records.append({
            "sample_prefix":     prefix,
            "subject_id":        subject_id,
            "subject_name":      first("subject_name"),
            "project_id":        first("project_short_name"),
            "anatomical_region": str(first("anatomical_site")).strip().lower()
                                 .replace(" ", "_").replace("/", "_"),
            "technique":         first("technique"),
            "chemistry":         first("chemistry"),
            "include_introns":   first("include_introns"),
            "fastq_convention":  convention,
            "fastq_dirs":        ",".join(dirs),
            "n_lane_dirs":       len(dirs),
            "nemo_library_ids":  ",".join(sorted(grp["library_id"].unique())),
        })

    samples = pd.DataFrame(records)

    # Ensure a globally-unique sample_tag (entity id). Disambiguate the rare
    # case of the same FASTQ prefix appearing under two subjects.
    dup_mask = samples.duplicated("sample_prefix", keep=False)
    samples["sample_tag"] = samples["sample_prefix"]
    if dup_mask.any():
        samples.loc[dup_mask, "sample_tag"] = (
            samples.loc[dup_mask, "sample_prefix"] + "__" +
            samples.loc[dup_mask, "subject_id"].astype(str)
        )
        warnings.warn(
            f"{int(dup_mask.sum())} sample prefixes collided across subjects; "
            f"suffixed with subject_id to keep sample_tag unique.",
            UserWarning,
        )

    # CellRanger --sample prefix (the real FASTQ token, without disambiguation).
    samples["fastq_sample_prefix"] = samples["sample_prefix"]

    # No hashing on this manifest path (hashing techniques are excluded).
    samples["hashing"]          = "false"
    samples["hashing_fastq_r1"] = ""
    samples["hashing_fastq_r2"] = ""
    samples["hashing_tag_file"] = ""

    print(f"      {len(df)} FASTQ rows → {len(samples)} CellRanger samples "
          f"across {samples['subject_id'].nunique()} subjects.")

    # Flag samples that need FASTQ normalisation before CellRanger.
    need_norm = samples[samples["fastq_convention"] != "standard"]
    if not need_norm.empty:
        print("      Samples needing FASTQ normalisation:")
        for conv, n in need_norm["fastq_convention"].value_counts().items():
            print(f"        - {conv}: {n} sample(s)")

    return samples


# ---------------------------------------------------------------------------
# Batch assignment
# ---------------------------------------------------------------------------

def assign_batches(samples: pd.DataFrame, batch_key: str,
                   label_col: str) -> pd.DataFrame:
    """Add batch_num (int) and sample_set_id (string) columns."""
    print(f"[3/5] Assigning batches using batch_key='{batch_key}' ...")
    if batch_key == "donor":
        group_cols = ["subject_id"]
        label_fn = lambda r: str(r[label_col]).replace(" ", "_")
    elif batch_key == "donor_region":
        group_cols = ["subject_id", "anatomical_region"]
        label_fn = lambda r: (
            f"{r[label_col]}_{r['anatomical_region']}".replace(" ", "_").replace("/", "_")
        )
    elif batch_key == "project":
        group_cols = ["project_id"]
        label_fn = lambda r: str(r["project_id"]).replace(":", "_").replace(" ", "_")
    else:
        raise ValueError(
            f"Unknown --batch_key '{batch_key}'. Choose: donor, donor_region, project"
        )

    keys = (samples[group_cols].drop_duplicates()
            .sort_values(group_cols).reset_index(drop=True))
    keys["batch_num"] = keys.index + 1
    samples = samples.merge(keys, on=group_cols, how="left")
    samples["sample_set_id"] = samples.apply(label_fn, axis=1)

    # Guarantee sample_set_id is 1:1 with batch_num. Human-readable labels
    # (e.g. subject_name) can collide across distinct batch keys (e.g. two
    # different subject_ids sharing a subject_name); left unresolved that would
    # silently merge separate donors/batches into one sample_set. Disambiguate
    # any colliding label by suffixing the underlying group key(s).
    label_to_batches = samples.groupby("sample_set_id")["batch_num"].nunique()
    collided = set(label_to_batches[label_to_batches > 1].index)
    if collided:
        mask = samples["sample_set_id"].isin(collided)
        suffix = samples.loc[mask, group_cols].astype(str).agg("_".join, axis=1)
        samples.loc[mask, "sample_set_id"] = (
            samples.loc[mask, "sample_set_id"] + "__" + suffix
        )
        warnings.warn(
            f"{len(collided)} sample_set label(s) collided across distinct "
            f"'{batch_key}' batches and were suffixed with {group_cols} to keep "
            f"sample_set_id unique (one submission per batch).",
            UserWarning,
        )

    print(f"      {samples['batch_num'].nunique()} batch(es).")
    return samples


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

SAMPLE_OUTPUT_COLS = [
    "sample_tag", "fastq_dirs", "fastq_sample_prefix", "fastq_convention",
    "chemistry", "include_introns", "hashing",
    "hashing_fastq_r1", "hashing_fastq_r2", "hashing_tag_file",
    "batch_num", "sample_set_id",
    "subject_id", "subject_name", "project_id", "anatomical_region",
    "technique", "n_lane_dirs", "nemo_library_ids",
]


def write_tables(samples: pd.DataFrame, output_dir: str, prefix: str,
                 write_individual_map: bool) -> None:
    print("[4/5] Writing output tables ...")
    os.makedirs(output_dir, exist_ok=True)

    sample_out = samples[SAMPLE_OUTPUT_COLS].copy()
    sample_out.insert(0, "entity:sample_id", samples["sample_tag"])

    sample_path = os.path.join(output_dir, f"{prefix}_sample_table.tsv")
    sample_out.to_csv(sample_path, sep="\t", index=False)
    print(f"      Wrote sample table ({len(sample_out)} rows): {sample_path}")

    set_rows = []
    for set_id, grp in sample_out.groupby("sample_set_id"):
        set_rows.append({
            "entity:sample_set_id": set_id,
            "samples": ",".join(grp["entity:sample_id"].tolist()),
            "batch_num": grp["batch_num"].iloc[0],
        })
    set_df = (pd.DataFrame(set_rows).sort_values("batch_num").reset_index(drop=True))
    set_path = os.path.join(output_dir, f"{prefix}_sample_set_table.tsv")
    set_df.to_csv(set_path, sep="\t", index=False)
    print(f"      Wrote sample_set table ({len(set_df)} batch(es)): {set_path}")

    if write_individual_map:
        ind_path = os.path.join(output_dir, f"{prefix}_sample_to_individual.tsv")
        ind = (sample_out[["sample_tag", "subject_id"]]
               .rename(columns={"subject_id": "individualID"})
               .drop_duplicates("sample_tag"))
        ind.to_csv(ind_path, sep="\t", index=False)
        print(f"      Wrote sample_to_individual map ({len(ind)} rows): {ind_path}")

    # Summary
    print("\n[5/5] Summary")
    print(f"  CellRanger samples:  {len(sample_out)}")
    print(f"  Subjects:            {samples['subject_id'].nunique()}")
    print(f"  Batches:             {set_df.shape[0]}")
    print("  Chemistry breakdown:")
    for chem, n in sample_out["chemistry"].value_counts().items():
        print(f"    {chem}: {n}")
    print("  FASTQ convention breakdown:")
    for conv, n in sample_out["fastq_convention"].value_counts().items():
        print(f"    {conv}: {n}")
    print("\nNext steps:")
    print(f"  1. Import '{sample_path}' into Terra as a 'sample' entity table.")
    print(f"  2. Import '{set_path}' into Terra as a 'sample_set' entity table.")
    print("  3. Set workflow inputs (transcriptome, mito_file, docker SHAs) and run.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate Terra sample / sample_set tables from the SCORCH "
                    "manifest (or a cohort subset of it).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--manifest", required=True,
                        help="Path to the tab-separated manifest / cohort subset.")
    parser.add_argument("--batch_key", default="donor",
                        choices=["donor", "donor_region", "project"],
                        help="How to group samples into sample_sets (default: donor).")
    parser.add_argument("--species", default="human",
                        help="Species filter (default: human; '*' disables).")
    parser.add_argument("--anatomical_sites",
                        default=",".join(DEFAULT_ANATOMICAL_SITES),
                        help="Comma-separated anatomical_site values to keep "
                             "(empty string disables site filtering).")
    parser.add_argument("--output_dir", default=".",
                        help="Directory for output TSVs (default: current dir).")
    parser.add_argument("--output_prefix", default="scorch",
                        help="Prefix for output filenames (default: scorch).")
    parser.add_argument("--batch_label_col", default="subject_name",
                        help="Column for human-readable sample_set_id labels "
                             "(default: subject_name).")
    parser.add_argument("--write_individual_map", action="store_true",
                        help="Also emit a sample_tag → individualID TSV for the "
                             "downstream Azimuth workflow.")
    args = parser.parse_args()

    sites = [s.strip() for s in args.anatomical_sites.split(",") if s.strip()]

    df = load_manifest(args.manifest, args.species, sites)
    samples = group_into_samples(df)
    samples = assign_batches(samples, args.batch_key, args.batch_label_col)
    write_tables(samples, args.output_dir, args.output_prefix,
                 args.write_individual_map)


if __name__ == "__main__":
    main()
