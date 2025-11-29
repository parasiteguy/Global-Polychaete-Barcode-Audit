#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#!/usr/bin/env python3
"""
Compute overall voucher and identifier completeness for extant polychaete COI records.

Inputs:
  - polychaete_COI_metadata_matched_extant.csv

Outputs:
  - Printed percentages:
      * % of COI records with a specimen voucher
      * % of COI records with an identified_by entry
"""

from pathlib import Path
import pandas as pd

# ==== EDIT THIS PATH ====
METADATA_FILE = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/polychaete_COI_metadata_matched_extant.csv")
# ========================

def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    required = ["worms_match_type", "specimen_voucher", "identified_by"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Metadata file missing required columns: {missing}")

    return df


def main():
    df = load_metadata(METADATA_FILE)
    print(f"Total COI records in file: {len(df):,}")

    # Keep only records that match WoRMS (accepted or synonym)
    df_valid = df[df["worms_match_type"] != "no_worms_match"].copy()
    n_total = len(df_valid)
    print(f"Valid COI records with WoRMS match: {n_total:,}")

    # Helper to check non-empty strings
    def non_empty(series):
        return series.notna() & (series.str.strip() != "")

    # Specimen voucher
    has_voucher = non_empty(df_valid["specimen_voucher"])
    n_voucher = has_voucher.sum()
    pct_voucher = 100 * n_voucher / n_total if n_total > 0 else 0.0

    # Identified_by
    has_identifier = non_empty(df_valid["identified_by"])
    n_identifier = has_identifier.sum()
    pct_identifier = 100 * n_identifier / n_total if n_total > 0 else 0.0

    print("\n=== Metadata completeness (record-level, extant-only) ===")
    print(f"Records with specimen voucher: {n_voucher:,}/{n_total:,} ({pct_voucher:.1f}%)")
    print(f"Records with identified_by:    {n_identifier:,}/{n_total:,} ({pct_identifier:.1f}%)")


if __name__ == "__main__":
    main()
