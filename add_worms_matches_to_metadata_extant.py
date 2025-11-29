#!/usr/bin/env python3

#!/usr/bin/env python3
"""
Add WoRMS (extant-only) match info to the existing COI metadata file.

Inputs:
  - worms_accepted_species_extant.csv   (extant accepted species from WoRMS)
  - worms_synonyms_extant.csv           (extant synonyms from WoRMS)
  - polychaete_COI_metadata_cleaned.csv (your existing cleaned COI metadata)

Output:
  - polychaete_COI_metadata_matched_extant.csv
    with two new columns:
      - worms_match_type: accepted_match / synonym_match / no_worms_match
      - worms_valid_name: valid WoRMS species name (if any)
"""

from pathlib import Path
import pandas as pd

# ====== EDIT THESE PATHS ======
ACCEPTED_EXTANT = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/worms_accepted_species_extant.csv")
SYNONYMS_EXTANT = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/worms_synonyms_extant.csv")
COI_METADATA_IN = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/polychaete_COI_metadata_cleaned.csv")
COI_METADATA_OUT = Path("polychaete_COI_metadata_matched_extant.csv")
# ==============================


def load_accepted(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    if "ScientificName" not in df.columns:
        raise ValueError("Accepted table missing 'ScientificName' column.")

    # Restrict to accepted species just in case
    if "Status" in df.columns:
        df = df[df["Status"] == "accepted"].copy()

    return df


def load_synonyms(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    required = ["ScientificName", "ValidName"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Synonyms table missing columns: {missing}")

    return df


def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    if "clean_name" not in df.columns:
        raise ValueError("Metadata missing 'clean_name' column.")
    return df


def build_name_maps(accepted_df: pd.DataFrame, syn_df: pd.DataFrame):
    """
    Build:
      - accepted_names: set of valid WoRMS names
      - synonym_to_valid: map synonym -> valid name
    """
    accepted_names = set(
        n for n in accepted_df["ScientificName"].dropna().unique() if n.strip() != ""
    )

    synonym_to_valid = {}
    for _, row in syn_df.iterrows():
        syn_name = row.get("ScientificName", "")
        valid_name = row.get("ValidName", "")
        if syn_name and valid_name:
            synonym_to_valid[syn_name] = valid_name

    return accepted_names, synonym_to_valid


def classify_records(meta_df, accepted_names, synonym_to_valid):
    """
    For each record, classify:
      - accepted_match: clean_name is an accepted WoRMS name
      - synonym_match:  clean_name is a synonym that maps to a valid name
      - no_worms_match: everything else
    """
    match_types = []
    valid_names = []

    for _, row in meta_df.iterrows():
        name = row.get("clean_name", "")
        if not name:
            match_types.append("no_worms_match")
            valid_names.append("")
            continue

        if name in accepted_names:
            match_types.append("accepted_match")
            valid_names.append(name)
        elif name in synonym_to_valid:
            match_types.append("synonym_match")
            valid_names.append(synonym_to_valid[name])
        else:
            match_types.append("no_worms_match")
            valid_names.append("")

    meta_df = meta_df.copy()
    meta_df["worms_match_type"] = match_types
    meta_df["worms_valid_name"] = valid_names
    return meta_df


def main():
    print("Loading WoRMS extant backbone...")
    accepted_df = load_accepted(ACCEPTED_EXTANT)
    syn_df = load_synonyms(SYNONYMS_EXTANT)

    print(f"Accepted extant species: {len(accepted_df):,}")
    print(f"Synonym rows (extant):   {len(syn_df):,}")

    print("Loading COI metadata...")
    meta_df = load_metadata(COI_METADATA_IN)
    print(f"COI metadata records:     {len(meta_df):,}")

    accepted_names, synonym_to_valid = build_name_maps(accepted_df, syn_df)

    print("Classifying records against WoRMS...")
    meta_matched = classify_records(meta_df, accepted_names, synonym_to_valid)

    # quick summary
    print("\nMatch type breakdown:")
    print(meta_matched["worms_match_type"].value_counts())

    meta_matched.to_csv(COI_METADATA_OUT, index=False)
    print(f"\nWrote matched metadata to: {COI_METADATA_OUT}")


if __name__ == "__main__":
    main()
