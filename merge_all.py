from pathlib import Path
import math

import pandas as pd
import yaml


# ============================================================
# CONFIGURATION
# ============================================================

BASE_DIR = Path(".")
OUTPUT_DIR = BASE_DIR / "output_cleaned"

SPECIES = [
    "pipos",
    "pineg",
    "kaonpos",
    "kaonneg",
]

# folder, filename tag, Q_min, Q_max
Q2_RANGES = [
    ("1to10",       "1to10",       1.0,   3.17),
    ("10to100",     "10to100",     3.16,  10.1),
    ("100to1000",   "100to1000",   10.0,  31.63),
    ("1000to10000", "1000to10000", 31.62, 100.0),
]

SUFFIX = "9x275_26_07"


# ============================================================
# YAML VALIDITY CHECK
# ============================================================

def contains_invalid_value(obj):
    """
    Recursively check a YAML object.

    Returns True if any value is:
      - None / null
      - NaN
      - +inf / -inf

    Numeric strings are converted to float when possible.
    Non-numeric strings are left untouched.
    """

    if isinstance(obj, dict):
        return any(contains_invalid_value(v) for v in obj.values())

    if isinstance(obj, list):
        return any(contains_invalid_value(v) for v in obj)

    if obj is None:
        return True

    try:
        value = float(obj)

        if not math.isfinite(value):
            return True

    except (TypeError, ValueError):
        # Ordinary non-numeric strings are fine
        pass

    return False


# ============================================================
# CSV MERGE
# ============================================================

def merge_csv(species):

    frames = []

    print(f"\n--- CSV: {species} ---")

    for folder, q2_tag, q_min, q_max in Q2_RANGES:

        filename = f"{species}_q2_{q2_tag}_{SUFFIX}.csv"
        filepath = BASE_DIR / folder / filename

        if not filepath.exists():
            raise FileNotFoundError(f"Missing file: {filepath}")

        df = pd.read_csv(
            filepath,
            sep=",",
            skipinitialspace=True
        )

        # Remove possible leading/trailing spaces from column names
        df.columns = df.columns.str.strip()

        # ----------------------------------------------------
        # Remove NaN
        # Equivalent validity requirement for CSV
        # ----------------------------------------------------

        df = df.dropna()

        # ----------------------------------------------------
        # Remove +/- inf from numeric columns
        # ----------------------------------------------------

        numeric_cols = df.select_dtypes(include="number").columns

        finite_mask = df[numeric_cols].apply(
            lambda col: col.map(math.isfinite)
        ).all(axis=1)

        df = df[finite_mask]

        # ----------------------------------------------------
        # Remove physically invalid cross-section bins
        # ----------------------------------------------------

        df = df[
            (df["xSec_diff[pb/GeV^3]"] > 0) &
            (df["xSec_uncorr_error"] > 0)
        ]

        n_before = len(df)

        # ----------------------------------------------------
        # Keep only bins completely inside nominal Q range
        # ----------------------------------------------------

        df = df[
            (df["Q_min[GeV]"] >= q_min) &
            (df["Q_max[GeV]"] <= q_max)
        ]

        n_after = len(df)

        print(
            f"{folder}: "
            f"{n_before} valid -> "
            f"{n_after} kept "
            f"(discarded {n_before - n_after})"
        )

        frames.append(df)

    # --------------------------------------------------------
    # Merge the four Q2 productions
    # --------------------------------------------------------

    merged = pd.concat(frames, ignore_index=True)

    # Renumber IDs after merge
    merged["id"] = range(1, len(merged) + 1)

    # Preserve compatibility with old code:
    # add initial space to column names
    merged.columns = [" " + c for c in merged.columns]

    output_file = (
        OUTPUT_DIR /
        f"{species}_merged_clean_{SUFFIX}.csv"
    )

    merged.to_csv(
        output_file,
        index=False,
        float_format="%.6e"
    )

    print(f"CSV written: {output_file}")
    print(f"Total bins: {len(merged)}")

    return len(merged)


# ============================================================
# YAML MERGE
# ============================================================

def merge_yaml(species):

    merged_bins = []

    print(f"\n--- YAML: {species} ---")

    for folder, q2_tag, q_min, q_max in Q2_RANGES:

        filename = f"{species}_q2_{q2_tag}_{SUFFIX}.yaml"
        filepath = BASE_DIR / folder / filename

        if not filepath.exists():
            raise FileNotFoundError(f"Missing file: {filepath}")

        with open(filepath, "r") as f:
            data = yaml.safe_load(f)

        if data is None:
            raise ValueError(f"Empty YAML file: {filepath}")

        n_before = 0
        n_after = 0

        for entry in data:

            # ------------------------------------------------
            # Basic YAML structure check
            # ------------------------------------------------

            if not isinstance(entry, dict):
                continue

            if "bin" not in entry:
                continue

            b = entry["bin"]

            # ------------------------------------------------
            # Equivalent to CSV dropna + finite check
            #
            # Reject bin if ANY field contains:
            #   null
            #   NaN
            #   +inf
            #   -inf
            # ------------------------------------------------

            if contains_invalid_value(b):
                continue

            # ------------------------------------------------
            # Explicit numeric conversion
            #
            # PyYAML may interpret scientific notation
            # as strings.
            # ------------------------------------------------

            try:
                xsec = float(b["xSec_diff"])
                xsec_err = float(b["xSec_uncorr_error"])

                Q_min = float(b["Q"]["Q_min"])
                Q_max = float(b["Q"]["Q_max"])

            except (TypeError, ValueError, KeyError):
                continue

            # ------------------------------------------------
            # Extra safety check
            # ------------------------------------------------

            if not math.isfinite(xsec):
                continue

            if not math.isfinite(xsec_err):
                continue

            if not math.isfinite(Q_min):
                continue

            if not math.isfinite(Q_max):
                continue

            # ------------------------------------------------
            # Same physical cuts as CSV
            # ------------------------------------------------

            if xsec <= 0:
                continue

            if xsec_err <= 0:
                continue

            n_before += 1

            # ------------------------------------------------
            # Keep only bins completely inside nominal Q range
            # ------------------------------------------------

            if not (
                Q_min >= q_min
                and Q_max <= q_max
            ):
                continue

            merged_bins.append(entry)
            n_after += 1

        print(
            f"{folder}: "
            f"{n_before} valid -> "
            f"{n_after} kept "
            f"(discarded {n_before - n_after})"
        )

    # --------------------------------------------------------
    # Renumber IDs after merge
    # --------------------------------------------------------

    for i, entry in enumerate(merged_bins, start=1):
        entry["bin"]["id"] = i

    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    output_file = (
        OUTPUT_DIR /
        f"{species}_merged_clean_{SUFFIX}.yaml"
    )

    with open(output_file, "w") as f:
        yaml.safe_dump(
            merged_bins,
            f,
            sort_keys=False,
            default_flow_style=False
        )

    print(f"YAML written: {output_file}")
    print(f"Total bins: {len(merged_bins)}")

    return len(merged_bins)


# ============================================================
# MAIN
# ============================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    print("========================================")
    print("   MERGING SIDIS EXTRACTION TABLES")
    print("========================================")

    results = {}

    for species in SPECIES:

        print("\n========================================")
        print(f" Processing {species}")
        print("========================================")

        n_csv = merge_csv(species)
        n_yaml = merge_yaml(species)

        results[species] = (n_csv, n_yaml)

        # ----------------------------------------------------
        # Immediate consistency check
        # ----------------------------------------------------

        print("\n--- CONSISTENCY CHECK ---")

        if n_csv == n_yaml:
            print(
                f"{species}: "
                f"CSV = {n_csv}, "
                f"YAML = {n_yaml}  [OK]"
            )
        else:
            print(
                f"{species}: "
                f"CSV = {n_csv}, "
                f"YAML = {n_yaml}  [MISMATCH]"
            )

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    print("\n========================================")
    print(" FINAL SUMMARY")
    print("========================================")

    all_ok = True

    for species, (n_csv, n_yaml) in results.items():

        status = "OK" if n_csv == n_yaml else "MISMATCH"

        if n_csv != n_yaml:
            all_ok = False

        print(
            f"{species:10s} "
            f"CSV: {n_csv:6d}   "
            f"YAML: {n_yaml:6d}   "
            f"[{status}]"
        )

    print("========================================")

    if all_ok:
        print("All CSV/YAML outputs are consistent.")
    else:
        print("WARNING: some CSV/YAML outputs differ.")

    print(f"\nOutput directory: {OUTPUT_DIR}")


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    main()
