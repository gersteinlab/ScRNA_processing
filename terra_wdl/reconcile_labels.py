"""
reconcile_labels.py

Reconciles AHBA and Ma-Sestan Azimuth predictions for the same set of cells.
Refactored from downstream/ahba_ma_map_v2.py for use in the downstream Terra
WDL workflow (terra_wdl/downstream.wdl).

Changes from the original:
- CLI args replace hardcoded `dataset_list = ['ptsd_brainomics_new']`.
- Removes the `do_remove_kept_analysis` exploratory branch (kept commented in
  original; not needed in production WDL).
- Outputs are written to the current working directory; the WDL `output {}`
  block captures them via globs.

Behavior:
1. Read AHBA and Ma-Sestan predictions CSVs (one row per cell barcode).
2. Map each subclass label to a cell class (Exc / Inh / Glial) via the
   schema dictionaries `ahba_map` and `ma_map`.
3. Keep cells where the two schemes agree at the class level.
4. For agreed cells: use AHBA labels for neurons (Exc / Inh), Ma-Sestan
   labels for glia.
5. Write reconciled labels + dropped (inconsistent) cells to separate CSVs.

Schema dictionaries (`ahba_map`, `ma_map`) are kept as module-level
constants since they encode cell-type taxonomy, not data-specific config.
"""

import argparse
import os

import pandas as pd


# ---------------------------------------------------------------------------
# Cell-type class schemas (subclass label -> broad class)
# ---------------------------------------------------------------------------

ahba_map = {
    "L6 IT": "Exc", "L5/6 NP": "Exc", "L2/3 IT": "Exc", "L5 IT": "Exc",
    "L6b": "Exc", "L5 ET": "Exc", "L4 IT": "Exc", "L6 IT Car3": "Exc",
    "L6 CT": "Exc",
    "Chandelier": "Inh", "Vip": "Inh", "Lamp5": "Inh", "Lamp5 Lhx6": "Inh",
    "Sst Chodl": "Inh", "Sst": "Inh", "Pvalb": "Inh", "Sncg": "Inh",
    "Pax6": "Inh",
    "Astro": "Glial", "Oligo": "Glial", "Micro/PVM": "Glial",
    "Endo": "Glial", "OPC": "Glial", "VLMC": "Glial",
}

ma_map = {
    "L5-6 NP": "Exc", "L6 IT-2": "Exc", "L6B": "Exc", "L3-5 IT-3": "Exc",
    "L5 ET": "Exc", "L3-5 IT-2": "Exc", "L6 IT-1": "Exc",
    "L3-5 IT-1": "Exc", "L2-3 IT": "Exc", "L6 CT": "Exc",
    "ADARB2 KCNG1": "Inh", "SST NPY": "Inh", "VIP": "Inh",
    "PVALB ChC": "Inh", "SST": "Inh", "LAMP5 RELN": "Inh",
    "SST HGF": "Inh", "PVALB": "Inh", "LAMP5 LHX6": "Inh",
    "PC": "Glial", "RB": "Glial", "SMC": "Glial", "Micro": "Glial",
    "Immune": "Glial", "Astro": "Glial", "Oligo": "Glial",
    "Endo": "Glial", "OPC": "Glial", "VLMC": "Glial",
}


def reconcile(ahba_csv: str, ma_sestan_csv: str, output_prefix: str) -> None:
    """
    Reconcile two Azimuth prediction CSVs.

    The input CSVs are expected to follow the format produced by
    Hybrid_Azimuth.R / hybrid_azimuth_task.R: first column is the cell
    barcode (index), second column is the predicted subclass label
    (originally named 'x' by Seurat's write.csv).
    """
    # ------------------------------------------------------------------
    # 1. Read and rename
    # ------------------------------------------------------------------
    df_ma = pd.read_csv(ma_sestan_csv, header=0, index_col=0)
    df_ma = df_ma.rename(mapper={"x": "Ma_Sestan"}, axis=1)
    df_ahba = pd.read_csv(ahba_csv, header=0, index_col=0)
    df_ahba = df_ahba.rename(mapper={"x": "AHBA"}, axis=1)

    print(f"AHBA predictions: {len(df_ahba)} cells")
    print(f"Ma-Sestan predictions: {len(df_ma)} cells")

    # ------------------------------------------------------------------
    # 2. Merge and assign broad classes
    # ------------------------------------------------------------------
    merged = pd.concat([df_ma, df_ahba], axis=1)
    merged["Ma_Sestan_class"] = merged["Ma_Sestan"].apply(
        lambda x: ma_map.get(x, "Unknown")
    )
    merged["AHBA_class"] = merged["AHBA"].apply(
        lambda x: ahba_map.get(x, "Unknown")
    )

    # ------------------------------------------------------------------
    # 3. Keep cells where classes agree
    # ------------------------------------------------------------------
    kept = merged[merged["Ma_Sestan_class"] == merged["AHBA_class"]].copy()
    removed = merged[merged["Ma_Sestan_class"] != merged["AHBA_class"]].copy()

    print(f"Class agreement: {len(kept)} cells kept, {len(removed)} dropped")
    if len(merged) > 0:
        drop_pct = 100.0 * len(removed) / len(merged)
        print(f"  Drop rate: {drop_pct:.2f}%")

    # ------------------------------------------------------------------
    # 4. Use AHBA for neurons, Ma-Sestan for glia
    # ------------------------------------------------------------------
    kept_neurons = kept[kept["AHBA_class"].isin(["Exc", "Inh"])]
    kept_glia = kept[kept["Ma_Sestan_class"] == "Glial"]

    neuron_labels = kept_neurons[["AHBA"]].rename(mapper={"AHBA": "x"}, axis=1)
    glia_labels = kept_glia[["Ma_Sestan"]].rename(
        mapper={"Ma_Sestan": "x"}, axis=1
    )

    reconciled = pd.concat([neuron_labels, glia_labels], axis=0)

    # ------------------------------------------------------------------
    # 5. Write outputs
    # ------------------------------------------------------------------
    reconciled_path = f"{output_prefix}_Azimuth_predictions_Ma_Sestan_AHBA_reconcile.csv"
    inconsistent_path = f"{output_prefix}_Azimuth_predictions_Ma_Sestan_AHBA_inconsistent.csv"

    reconciled.to_csv(reconciled_path)
    removed.to_csv(inconsistent_path)

    print(f"Wrote reconciled labels: {reconciled_path}")
    print(f"Wrote inconsistent cells: {inconsistent_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile AHBA and Ma-Sestan Azimuth predictions. "
            "Keeps cells where the two atlases agree on broad class "
            "(Exc/Inh/Glial); picks AHBA labels for neurons and "
            "Ma-Sestan labels for glia."
        )
    )
    parser.add_argument(
        "--ahba_csv", required=True,
        help="Path to AHBA Azimuth predictions CSV (one row per cell)."
    )
    parser.add_argument(
        "--ma_sestan_csv", required=True,
        help="Path to Ma-Sestan Azimuth predictions CSV (one row per cell)."
    )
    parser.add_argument(
        "--output_prefix", required=True,
        help="Output filename prefix (e.g. batch name)."
    )
    args = parser.parse_args()

    reconcile(args.ahba_csv, args.ma_sestan_csv, args.output_prefix)


if __name__ == "__main__":
    main()
