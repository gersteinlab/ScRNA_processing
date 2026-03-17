"""
cellbender_task.py

WDL-compatible CellBender remove-background runner.
Called by the CellBender task in pipeline.wdl.

Replaces the sbatch submission loop in PEC_scRNA_pipeline_CellBender.py.

Key differences from the original:
  - Accepts individual CLI arguments instead of a JSON config file.
  - Runs cellbender remove-background directly (no Slurm sbatch).
  - GPU is made available by the WDL runtime block; this script just calls --cuda.
  - No multiprocessing.Pool: WDL scatter handles per-sample parallelism.
  - Writes the filtered H5 output path to a .txt file for WDL read_string().

Usage (called from WDL command block):
    python /opt/pipeline/cellbender_task.py \
        --input_h5        /cromwell_root/.../raw_feature_bc_matrix.h5 \
        --output_dir      ./cellbender_outputs \
        --expected_cells  5000 \
        --epochs          150 \
        --fpr             0.01
"""

import argparse
import os
import subprocess
import sys


def run_cellbender(args) -> str:
    """
    Run cellbender remove-background and return the path to the filtered H5.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    output_h5 = os.path.join(args.output_dir, "cellbender-output.h5")
    filtered_h5 = os.path.join(args.output_dir, "cellbender-output_filtered.h5")

    # total_droplets_included follows the original pipeline convention:
    # expected_cells + 20,000
    total_droplets = args.expected_cells + 20000

    cmd = [
        "cellbender", "remove-background",
        f"--input={args.input_h5}",
        f"--output={output_h5}",
        f"--expected-cells={args.expected_cells}",
        f"--total-droplets-included={total_droplets}",
        f"--fpr={args.fpr}",
        f"--epochs={args.epochs}",
        "--cuda",
    ]

    # MKL_THREADING_LAYER=GNU is required by CellBender on GPU (from original pipeline)
    env = os.environ.copy()
    env["MKL_THREADING_LAYER"] = "GNU"

    print(f"[cellbender_task] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, env=env)

    if not os.path.exists(filtered_h5):
        sys.exit(
            f"ERROR: CellBender did not produce the expected filtered output: {filtered_h5}"
        )

    print(f"[cellbender_task] CellBender complete. Filtered H5: {filtered_h5}")
    return filtered_h5


def main():
    parser = argparse.ArgumentParser(
        description="WDL-compatible CellBender remove-background runner for Terra."
    )
    parser.add_argument(
        "--input_h5", required=True,
        help="Path to the CellRanger raw_feature_bc_matrix.h5 input file."
    )
    parser.add_argument(
        "--output_dir", required=False, default="./cellbender_outputs",
        help="Directory for CellBender outputs (default: ./cellbender_outputs)."
    )
    parser.add_argument(
        "--expected_cells", required=True, type=int,
        help="Estimated number of real cells (from CellRanger estimated_cells)."
    )
    parser.add_argument(
        "--epochs", required=False, type=int, default=150,
        help="Number of training epochs for CellBender (default: 150)."
    )
    parser.add_argument(
        "--fpr", required=False, type=float, default=0.01,
        help="False positive rate for CellBender (default: 0.01)."
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Run CellBender
    # ------------------------------------------------------------------
    filtered_h5 = run_cellbender(args)

    # ------------------------------------------------------------------
    # Write output path to file for WDL read_string()
    # ------------------------------------------------------------------
    with open("cb_filtered_h5_path.txt", "w") as f:
        f.write(filtered_h5)

    print(f"[cellbender_task] Done.")
    print(f"  filtered_h5: {filtered_h5}")


if __name__ == "__main__":
    main()
