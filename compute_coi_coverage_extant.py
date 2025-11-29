#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#!/usr/bin/env python3
"""
Compute COI sequence counts per species using the extant-only WoRMS backbone
and the matched COI metadata.

Inputs:
  - worms_accepted_species_extant.csv
  - polychaete_COI_metadata_matched_extant.csv

Output:
  - COI_counts_by_species_extant.csv

Columns in output:
  - ScientificName   (WoRMS accepted species name)
  - Family
  - Genus
  - n_COI_sequences  (number of COI records in GenBank matched to this species)
"""

import pandas as pd
from pathlib import Path

# ===== EDIT THESE PATHS IF NEEDED =====
ACCEPTED_EXTANT = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/worms_accepted_species_extant.csv")
COI_META = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/polychaete_COI_metadata_matched_extant.csv")
OUTFILE = "COI_counts_by_species_extant.csv"
# =====================================


def load_accepted(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Keep only accepted species (extant already filtered upstream)
    if "Status" in df.columns:
        df = df[df["Status"] == "accepted"].copy()

    required = ["ScientificName", "Family", "Genus"]
    for r in required:
        if r not in df.columns:
            raise ValueError(f"Accepted species table missing required column: {r}")

    print(f"Loaded {len(df):,} extant accepted species from WoRMS.")
    return df


def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Keep only records that match WoRMS (accepted or synonym)
    if "worms_match_type" not in df.columns or "worms_valid_name" not in df.columns:
        raise ValueError("Metadata must contain 'worms_match_type' and 'worms_valid_name'.")

    df = df[df["worms_match_type"] != "no_worms_match"].copy()

    print(f"Loaded {len(df):,} COI records with a WoRMS match.")
    return df


def compute_counts_by_species(accepted_df: pd.DataFrame, coi_df: pd.DataFrame) -> pd.DataFrame:
    # Count how many COI records per valid WoRMS species name
    counts = (
        coi_df["worms_valid_name"]
        .dropna()
        .value_counts()
        .rename_axis("ScientificName")
        .reset_index(name="n_COI_sequences")
    )

    print(f"Found COI data for {len(counts):,} species (with ≥1 sequence).")

    # Merge onto accepted species table so species without COI get n_COI_sequences = 0
    merged = accepted_df.merge(counts, on="ScientificName", how="left")
    merged["n_COI_sequences"] = merged["n_COI_sequences"].fillna(0).astype(int)

    return merged


def main():
    accepted_df = load_accepted(ACCEPTED_EXTANT)
    coi_df = load_metadata(COI_META)

    species_counts = compute_counts_by_species(accepted_df, coi_df)

    species_counts.to_csv(OUTFILE, index=False)
    print(f"\nWrote species-level COI counts to: {OUTFILE}")

    # Quick summary
    total_species = len(species_counts)
    with_coi = (species_counts["n_COI_sequences"] > 0).sum()
    print(f"Species with ≥1 COI sequence: {with_coi:,}/{total_species:,} "
          f"({100 * with_coi / total_species:.2f}%)")

    print("\nTop 10 most heavily sampled species:")
    print(
        species_counts.sort_values("n_COI_sequences", ascending=False)
        .head(10)[["ScientificName", "Family", "Genus", "n_COI_sequences"]]
        .to_string(index=False)
    )

    print("\nExamples of species with no COI sequences:")
    print(
        species_counts[species_counts["n_COI_sequences"] == 0]
        .head(10)[["ScientificName", "Family", "Genus"]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()

