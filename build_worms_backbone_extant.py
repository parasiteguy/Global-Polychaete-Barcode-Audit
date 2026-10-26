#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 19:40:51 2025

@author: w00451725
"""

#!/usr/bin/env python3
"""
Build an *extant-only* WoRMS Polychaeta taxonomy backbone from a CSV export.

Changes vs. original:
- Filters out fossil-only taxa using the 'Recent/Fossil' column.
  Keeps only rows that contain 'recent' (e.g., 'recent', 'recent + fossil').
- Then:
  1. Filters to species-level records.
  2. Splits into accepted vs. non-accepted species.
  3. Builds synonym mapping.
  4. Writes extant-only accepted species & synonyms to CSV and SQLite.
"""

import argparse
import sqlite3
from pathlib import Path

import pandas as pd


def load_worms_csv(path: Path) -> pd.DataFrame:
    """Load the WoRMS Polychaeta CSV."""
    df = pd.read_csv(path, dtype=str)
    # Strip whitespace
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def filter_extant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only extant (recent) taxa based on 'Recent/Fossil' column.

    Logic:
      - If 'Recent/Fossil' column exists:
          Keep rows where value contains 'recent' (case-insensitive),
          e.g. 'recent', 'recent + fossil'.
          This drops pure fossil taxa ('fossil').
      - If column missing, return df unchanged (but warn).
    """
    col = None
    # Try to find the column even if slightly mis-named
    for c in df.columns:
        if c.lower().replace(" ", "") in ("recent/fossil", "recent_fossil", "recentfossil"):
            col = c
            break

    if col is None:
        print("WARNING: No 'Recent/Fossil' column detected. No fossil filtering applied.")
        return df

    before = len(df)
    mask = df[col].astype(str).str.lower().str.contains("recent")
    df_extant = df[mask].copy()
    after = len(df_extant)

    print(f"Fossil filter: kept {after:,} of {before:,} rows "
          f"({before - after:,} fossil-only rows removed) based on '{col}'.")
    return df_extant


def filter_species(df: pd.DataFrame) -> pd.DataFrame:
    """Return only species-rank records."""
    if "RankName" not in df.columns:
        raise ValueError("Expected a 'RankName' column in the WoRMS CSV.")
    species_df = df[df["RankName"] == "Species"].copy()
    return species_df


def split_accepted_and_synonyms(species_df: pd.DataFrame):
    """Split species-level records into accepted and non-accepted."""
    if "Status" not in species_df.columns:
        raise ValueError("Expected a 'Status' column in the WoRMS CSV.")

    accepted = species_df[species_df["Status"] == "accepted"].copy()
    synonyms = species_df[species_df["Status"] != "accepted"].copy()
    return accepted, synonyms


def build_synonym_table(synonyms_df: pd.DataFrame) -> pd.DataFrame:
    """Build a clean synonym mapping table."""
    required_cols = [
        "ScientificName",
        "Authority",
        "AphiaID",
        "Status",
        "ValidName",
        "ValidName_AphiaID",
    ]
    missing = [c for c in required_cols if c not in synonyms_df.columns]
    if missing:
        raise ValueError(f"Synonyms DataFrame is missing required columns: {missing}")

    cols_to_keep = [
        "ScientificName",
        "Authority",
        "AphiaID",
        "Status",
        "ReasonSynonymy",
        "ValidName",
        "ValidName_AphiaID",
        "OriginalName",
        "OriginalName_AphiaID",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Marine",
        "Brackish",
        "Fresh",
        "Terrestrial",
        "Recent/Fossil",
        "qualitystatus_name",
    ]
    cols_to_keep = [c for c in cols_to_keep if c in synonyms_df.columns]

    syn = synonyms_df[cols_to_keep].copy().drop_duplicates()
    return syn


def build_accepted_table(accepted_df: pd.DataFrame) -> pd.DataFrame:
    """Build a cleaned table of accepted species (extant only)."""
    cols_to_keep = [
        "ScientificName",
        "Authority",
        "AphiaID",
        "Status",
        "ValidName",
        "ValidName_AphiaID",
        "OriginalName",
        "OriginalName_AphiaID",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Marine",
        "Brackish",
        "Fresh",
        "Terrestrial",
        "Recent/Fossil",
        "qualitystatus_name",
    ]
    cols_to_keep = [c for c in cols_to_keep if c in accepted_df.columns]

    acc = accepted_df[cols_to_keep].copy()
    acc = acc.drop_duplicates(subset=["AphiaID"])
    return acc


def write_csv(df: pd.DataFrame, path: Path):
    df.to_csv(path, index=False)
    print(f"Wrote {len(df):,} rows to {path}")


def write_sqlite(accepted_df: pd.DataFrame, synonyms_df: pd.DataFrame, path: Path):
    conn = sqlite3.connect(path)
    accepted_df.to_sql("accepted_species", conn, if_exists="replace", index=False)
    synonyms_df.to_sql("synonyms", conn, if_exists="replace", index=False)
    conn.close()
    print(f"Created SQLite DB at {path} with tables 'accepted_species' and 'synonyms'")


def main():
    parser = argparse.ArgumentParser(
        description="Build extant-only WoRMS Polychaeta taxonomy backbone from CSV export."
    )
    parser.add_argument(
        "worms_csv",
        type=str,
        help="Path to WoRMS Polychaeta CSV file (exported from Excel).",
    )
    parser.add_argument(
        "--accepted-csv",
        type=str,
        default="worms_accepted_species_extant.csv",
        help="Output CSV for accepted extant species "
             "(default: worms_accepted_species_extant.csv)",
    )
    parser.add_argument(
        "--synonyms-csv",
        type=str,
        default="worms_synonyms_extant.csv",
        help="Output CSV for extant synonyms "
             "(default: worms_synonyms_extant.csv)",
    )
    parser.add_argument(
        "--sqlite",
        type=str,
        default="worms_polychaeta_backbone_extant.sqlite",
        help="Output SQLite DB filename "
             "(default: worms_polychaeta_backbone_extant.sqlite)",
    )

    args = parser.parse_args()

    worms_path = Path(args.worms_csv)
    if not worms_path.exists():
        raise SystemExit(f"Input file not found: {worms_path}")

    print(f"Loading WoRMS CSV from {worms_path} ...")
    worms_df = load_worms_csv(worms_path)
    print(f"Loaded {len(worms_df):,} total rows.")

    # NEW: filter out fossil-only taxa
    worms_extant = filter_extant(worms_df)

    # Filter to species-level records
    species_df = filter_species(worms_extant)
    print(f"Filtered to {len(species_df):,} extant species-rank rows.")

    # Split accepted vs synonyms
    accepted_df, synonyms_df = split_accepted_and_synonyms(species_df)
    print(f"Accepted extant species: {len(accepted_df):,}")
    print(f"Non-accepted extant species (synonyms, etc.): {len(synonyms_df):,}")

    # Build cleaned tables
    accepted_clean = build_accepted_table(accepted_df)
    synonyms_clean = build_synonym_table(synonyms_df)

    # Write CSVs
    write_csv(accepted_clean, Path(args.accepted_csv))
    write_csv(synonyms_clean, Path(args.synonyms_csv))

    # Write SQLite DB
    write_sqlite(accepted_clean, synonyms_clean, Path(args.sqlite))


if __name__ == "__main__":
    main()
