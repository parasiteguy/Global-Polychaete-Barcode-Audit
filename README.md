# Global Polychaete COI Barcode Audit

**The missing polychaete barcodes: a global assessment of COI coverage and metadata gaps in GenBank**

**Authors:** Andrew A. Davinack and Charlotte K. Daniels  
**Project repository:** https://github.com/parasiteguy/Global-Polychaete-Barcode-Audit

## Overview

This repository contains code and data associated with a global audit of mitochondrial cytochrome *c* oxidase subunit I (COI) records for Polychaeta in NCBI GenBank. The study compares GenBank taxonomic labels with an extant-species backbone derived from the World Register of Marine Species (WoRMS), measures the proportion of accepted species represented by at least one COI accession, examines uneven sequence representation, and assesses reporting of the GenBank `specimen_voucher` and `identified_by` fields.

The analyses concern **GenBank coverage of accepted extant species**. They do not estimate the proportion of undescribed polychaete diversity that has been sequenced, nor do they independently validate specimen identifications or nucleotide-sequence quality.

The GenBank retrieval used for the original study was performed in **November 2025**. WoRMS taxonomy and GenBank content may have changed since that snapshot.

## Summary of the final screened analysis

The revised analytical dataset discussed in the manuscript has the following record counts:

| Processing stage | Count |
| --- | ---: |
| Initially retrieved GenBank entries | 29,939 |
| Unique accessions after removing 476 repeated entries | 29,463 |
| Accessions after excluding 8 entries without identifiable COI annotations | 29,455 |
| Accessions assigned to species in the final accepted WoRMS backbone | 21,235 |
| Accepted extant polychaete species in the backbone | 11,918 |
| Species with at least one matched COI accession | 1,850 (15.5%) |

Among the 21,235 final species-assigned accessions, 10,257 (48.3%) contained a nonempty `specimen_voucher` field and 4,002 (18.8%) contained a nonempty `identified_by` field. Field presence measures **metadata reporting**, not confirmation that a voucher is accessible or a taxonomic determination is correct.

**Version warning:** Several currently hosted CSVs and scripts predate the final deduplication, COI-annotation screening, and accepted-species filtering. They are retained as historical/intermediate files and **must not be treated as the final screened dataset**. The final analysis needs its own archived accession list, processed metadata, code, and figure inputs before the repository can be described as fully end-to-end reproducible.

## File guide

The following files are in the repository's current top-level collection (or are the named files associated with the original analysis).

### Taxonomic backbone

| File | Purpose |
| --- | --- |
| [`build_worms_backbone_extant.py`](build_worms_backbone_extant.py) | Reads an original WoRMS Polychaeta CSV export, keeps species-rank entries marked as recent, separates accepted from non-accepted names, and writes accepted-species and non-accepted-name CSVs plus a SQLite database. |
| [`worms_accepted_species_extant.csv`](worms_accepted_species_extant.csv) | Fixed reference list of 11,918 accepted extant species used as the species-coverage denominator, with taxonomic attributes and AphiaIDs. |
| [`worms_synonyms_extant.csv`](worms_synonyms_extant.csv) | Table of non-accepted names and their available WoRMS valid-name information used during name reconciliation. Despite its filename, it also contains categories other than straightforward synonyms (for example, superseded combinations and names of uncertain status). Not every entry resolves to a species in the fixed accepted-species list. |

The backbone builder can also create `worms_polychaeta_backbone_extant.sqlite`, with `accepted_species` and `synonyms` tables. The SQLite database is a convenient snapshot, but the corresponding CSVs are sufficient for the existing CSV-based downstream scripts.

### GenBank retrieval and reconciliation

| File | Purpose |
| --- | --- |
| [`coi_sequence_mining.py`](coi_sequence_mining.py) | Updated, restartable Entrez retrieval implementation that searches accepted genera for multiple COI gene/title annotations and writes GenBank-format records with retrieval documentation. It preserves the original search expression by default but is **not** the identical script run in November 2025. |
| [`add_worms_matches_to_metadata_extant.py`](add_worms_matches_to_metadata_extant.py) | Matches cleaned GenBank species names against accepted WoRMS names and non-accepted-name mappings; assigns `accepted_match`, `synonym_match`, or `no_worms_match`. Its input `polychaete_COI_metadata_cleaned.csv` is not included in the file list above. |
| [`polychaete_COI_metadata_matched_extant.csv`](polychaete_COI_metadata_matched_extant.csv) | Original taxonomically reconciled metadata output, including accession, species-identification fields, and WoRMS matching results. This is a **pre-final-filtering** intermediate dataset, not the revised 21,235-accession dataset. |
| [`COI_accessions_valid_species.csv`](COI_accessions_valid_species.csv) | Historical one-column accession list associated with the preliminary accepted/synonym-matched records. It predates the final accession deduplication and COI screening; do not use it unmodified as the final accession inventory. |

### Coverage and metadata summaries

| File | Purpose |
| --- | --- |
| [`compute_coi_coverage_extant.py`](compute_coi_coverage_extant.py) | Produces per-species COI record counts by joining matched metadata to the accepted-species backbone. The currently hosted version works from the earlier dataset and should be updated to count screened, unique accessions. |
| [`summarize_COI_metadata_by_family_genus.py`](summarize_COI_metadata_by_family_genus.py) | Adds accepted family/genus information to matched COI records and summarizes accession/record counts and metadata reporting by family and genus. The existing version uses preliminary matches and does not implement all final filters. |
| [`voucherlinkedcomputation.py`](voucherlinkedcomputation.py) | Original standalone calculation of overall `specimen_voucher` and `identified_by` field reporting among preliminary WoRMS-matched records. Retained for methodological history; its original percentages should not replace the corrected manuscript percentages. |
| [`README.md`](README.md) | Project description, file provenance, methodology, and reproducibility status (this document). |

## Workflow and execution

The intended analysis has five stages:

1. Build an extant-only WoRMS taxonomic backbone from the original WoRMS CSV export.
2. Retrieve COI-like GenBank records for the accepted genera, then parse record-level metadata.
3. Reconcile GenBank taxonomic names against the accepted and non-accepted WoRMS names.
4. Remove repeated accessions and records lacking identifiable COI annotations; retain only records assigned to species in the fixed accepted-species list for species-level summaries.
5. Calculate coverage and metadata reporting, then generate the figures and supplementary tables.

### Software

The original analyses used Python 3.10 with pandas and Biopython; plotting used common scientific Python libraries. For basic backbone preparation and sequence retrieval:

```bash
python -m pip install pandas biopython
```

Additional numerical/plotting dependencies may be necessary for later analyses. Version-pinned dependencies and a fully tested one-command pipeline are **not yet provided** here.

### 1. Rebuild the WoRMS backbone (requires the original source export)

```bash
python build_worms_backbone_extant.py path/to/original_worms_polychaeta_export.csv
```

By default, this writes `worms_accepted_species_extant.csv`, `worms_synonyms_extant.csv`, and `worms_polychaeta_backbone_extant.sqlite`. Because the original WoRMS export is not among the files listed above, users can inspect the deposited fixed CSVs without rebuilding them.

### 2. Run the updated GenBank retrieval (creates a **new**, time-dependent dataset)

```bash
export NCBI_EMAIL="researcher@example.org"
python coi_sequence_mining.py \
  --accepted worms_accepted_species_extant.csv \
  --outdir coi_retrieval
```

The program writes `polychaete_COI_raw.gb` and associated retrieval logs/metadata within the requested output directory. NCBI responses reflect the database at the **new run date** and should not be assumed to match the November 2025 study snapshot. An NCBI API key, if used, must be supplied securely through the environment and never committed to GitHub.

### 3–5. Reconciliation, final screening, and summaries

The original scripts `add_worms_matches_to_metadata_extant.py`, `compute_coi_coverage_extant.py`, `summarize_COI_metadata_by_family_genus.py`, and `voucherlinkedcomputation.py` contain machine-specific input paths. Those paths must be updated before local execution. The first three were written for the earlier intermediate files and do **not**, by themselves, reproduce all the final screened figures and statistics above.

The original GenBank-to-metadata parsing script, the November 2025 GenBank snapshot, and a consolidated final filtering/plotting implementation must be archived or linked before claiming complete computational reproduction of the revised results. In particular, the final screening removed 476 repeated accession entries and eight records without identifiable COI annotations and enforced membership in the fixed accepted-species backbone. This was an annotation-based screening procedure, not independent molecular verification of every sequence.

## Interpretation and limitations

- **Coverage:** A species was considered represented if at least one retained COI-associated accession could be reconciled with the fixed accepted-species backbone. A complete mitochondrial genome with a COI annotation was counted as one accession; standardized barcode regions were not necessarily extracted.
- **Taxonomy:** Nomenclatural matching to WoRMS is not proof that the deposited specimen was correctly identified. Provisional `sp.`, `cf.`, and `aff.` designations can reflect legitimate taxonomic uncertainty.
- **Retrieval:** Searches based on accepted genera may miss records indexed only under historical genus combinations. Records outside GenBank (including BOLD-only records) are not included.
- **Metadata:** Nonempty `specimen_voucher` and `identified_by` fields indicate reporting in GenBank, not independently verified deposition or taxonomic expertise.
- **Reproducibility:** The accepted-species and synonym CSVs document the historical taxonomic snapshot; rerunning retrievals against live databases can alter counts.

## Data provenance and remaining archival components

Source databases: [WoRMS](https://www.marinespecies.org/) and [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/). The accepted-species and non-accepted-name tables originate from a WoRMS export supplied by VLIZ. The original full WoRMS export and the original `polychaete_COI_raw.gb` snapshot should be preserved with retrieval/export-date information and, where appropriate, deposited in a research data archive and linked here.

To complete the release accompanying the revised manuscript, archive the **original parsing script**, **final cleaned/species-assigned metadata**, **excluded-accession and reconciliation-exception tables**, **updated coverage and metadata summaries**, and **final figure-generation code**. The current files document substantial parts of the workflow but should not be represented as an entirely verified, one-command reproduction package.

## Citation and contact

Please cite the associated publication when bibliographic details are available, and cite WoRMS and NCBI GenBank appropriately for their underlying data. Until publication details are finalized, refer to the project as:

Davinack, A. A., and Daniels, C. K. *Major gaps in COI barcode coverage and taxonomic metadata across global polychaete diversity*. Manuscript in preparation/submission.

For questions or corrections, open an issue in this repository.
