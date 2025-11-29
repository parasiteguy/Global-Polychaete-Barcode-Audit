

#!/usr/bin/env python3
"""
Merge WoRMS taxonomy into COI metadata and summarize by family & genus.

Inputs:
  - worms_accepted_species_extant.csv
  - polychaete_COI_metadata_matched_extant.csv

Outputs:
  - polychaete_COI_metadata_matched_extant_with_taxa.csv
  - COI_metadata_by_family_extant.csv
  - COI_metadata_by_genus_extant.csv
"""

from pathlib import Path
import pandas as pd

# ==== EDIT THESE PATHS IF NEEDED ====
ACCEPTED_EXTANT = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/worms_accepted_species_extant.csv")
META_IN        = Path("/Users/w00451725/Desktop/2026 Projects/Polychaete database project (Davinack and Daniels)/New FIles (without fossils)/polychaete_COI_metadata_matched_extant.csv")
META_WITH_TAXA = "polychaete_COI_metadata_matched_extant_with_taxa.csv"
FAM_OUT        = "COI_metadata_by_family_extant.csv"
GEN_OUT        = "COI_metadata_by_genus_extant.csv"
# ====================================


def load_accepted(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    required = ["ScientificName", "Family", "Genus"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Accepted species table missing columns: {missing}")

    # Restrict to accepted extant species if Status present
    if "Status" in df.columns:
        df = df[df["Status"] == "accepted"].copy()

    print(f"Loaded {len(df):,} extant accepted species from WoRMS.")
    return df[["ScientificName", "Family", "Genus"]]


def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    required = ["worms_match_type", "worms_valid_name",
                "specimen_voucher", "identified_by"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Metadata file missing required columns: {missing}")

    # Keep only records that matched WoRMS (accepted or synonym)
    df = df[df["worms_match_type"] != "no_worms_match"].copy()

    print(f"Loaded {len(df):,} valid COI records with WoRMS match.")
    return df


def non_empty(series: pd.Series) -> pd.Series:
    """Return boolean mask of non-empty, non-NA strings."""
    return series.notna() & (series.str.strip() != "")


def summarize_by_group(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    groups = df.groupby(group_col, dropna=True)

    rows = []
    for grp, sub in groups:
        n = len(sub)
        n_species = sub["worms_valid_name"].dropna().nunique()

        v_mask = non_empty(sub["specimen_voucher"])
        id_mask = non_empty(sub["identified_by"])

        n_voucher = v_mask.sum()
        n_identifier = id_mask.sum()

        voucher_pct = 100 * n_voucher / n if n > 0 else 0.0
        id_pct = 100 * n_identifier / n if n > 0 else 0.0

        rows.append(
            {
                group_col: grp,
                "n_records": n,
                "n_species": n_species,
                "n_with_voucher": n_voucher,
                "voucher_pct": voucher_pct,
                "n_with_identified_by": n_identifier,
                "identified_by_pct": id_pct,
            }
        )

    out = pd.DataFrame(rows)
    out = out.sort_values("n_records", ascending=False)
    return out


def main():
    # 1. Load inputs
    accepted_df = load_accepted(ACCEPTED_EXTANT)
    meta_df = load_metadata(META_IN)

    # 2. Merge Family/Genus from WoRMS onto metadata
    merged = meta_df.merge(
        accepted_df,
        how="left",
        left_on="worms_valid_name",
        right_on="ScientificName",
        suffixes=("", "_worms")
    )

    # Drop duplicate ScientificName column from WoRMS
    merged = merged.drop(columns=["ScientificName"], errors="ignore")

    # Save enriched metadata (with Family/Genus)
    merged.to_csv(META_WITH_TAXA, index=False)
    print(f"\nWrote merged metadata with taxonomy to: {META_WITH_TAXA}")

    # 3. Summarize by Family
    fam_summary = summarize_by_group(merged, "Family")
    fam_summary.to_csv(FAM_OUT, index=False)
    print(f"\nFamily-level metadata summary written to: {FAM_OUT}")
    print(fam_summary.head(10).to_string(index=False))

    # 4. Summarize by Genus
    gen_summary = summarize_by_group(merged, "Genus")
    gen_summary.to_csv(GEN_OUT, index=False)
    print(f"\nGenus-level metadata summary written to: {GEN_OUT}")
    print(gen_summary.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
