#!/usr/bin/env python3
"""Retrieve polychaete COI GenBank records using a WoRMS genus backbone.

This is a reproducible replacement for the original genus-by-genus NCBI
E-utilities retrieval script. By default, the *original* COI search expression
and accepted-genus sampling strategy are retained.

Usage:
    export NCBI_EMAIL='researcher@example.org'
    python coi_sequence_mining.py --accepted worms_accepted_species.csv \
        --outdir coi_retrieval

See README_coi_sequence_mining.md for details and limitations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import logging
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError

import Bio
from Bio import Entrez, SeqIO

# Intentionally unchanged from the original script: changes to search terms
# can change the dataset and must be recorded as a new retrieval strategy.
COI_QUERY = (
    'COI[Gene] OR CO1[Gene] OR COXI[Gene] OR COX1[Gene] OR cox1[Gene] OR coxI[Gene] '
    'OR "cytochrome c oxidase subunit I"[Title] '
    'OR "cytochrome oxidase subunit I"[Title] '
    'OR "cytochrome c oxidase subunit 1"[Title] '
    'OR "cytochrome oxidase subunit 1"[Title] '
    'OR COI[Title] OR CO1[Title]'
)
OUTPUT_NAME = "polychaete_COI_raw.gb"
VERSION = "1.0.0"
LOG = logging.getLogger("coi_sequence_mining")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_json_atomic(path: Path, content: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(content, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(tmp, path)


def load_genera(path: Path) -> list[str]:
    """Read distinct genera from a CSV containing a `Genus` column."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        if not rows.fieldnames or "Genus" not in rows.fieldnames:
            raise ValueError(f"{path.name} must contain a column named 'Genus'")
        names = {(row.get("Genus") or "").strip() for row in rows}
    genera = sorted((g for g in names if g), key=lambda g: (g.casefold(), g))
    if not genera:
        raise ValueError(f"No genera found in {path.name}")
    return genera


def load_extra_genera(path: Path) -> list[str]:
    """Optional additional historical genera: Genus CSV or one name per line."""
    if path.suffix.lower() == ".csv":
        return load_genera(path)
    with path.open("r", encoding="utf-8-sig") as handle:
        names = [line.strip() for line in handle]
    names = [name for name in names if name and not name.startswith("#")]
    if not names:
        raise ValueError(f"No additional genera found in {path.name}")
    return sorted(set(names), key=lambda g: (g.casefold(), g))


def genus_stem(genus: str) -> str:
    """Safe deterministic filename, even if genus contains unusual characters."""
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", genus).strip("_")[:64] or "genus"
    fingerprint = hashlib.sha256(genus.encode("utf-8")).hexdigest()[:10]
    return f"{slug}_{fingerprint}"


def query_for_genus(genus: str) -> str:
    # Quote an organism name as in the original retrieval; do not change the
    # taxonomic search semantics without explicitly documenting a new run.
    if '"' in genus or "\n" in genus or "\r" in genus:
        raise ValueError(f"Invalid genus for NCBI query: {genus!r}")
    return f'"{genus}"[Organism] AND ({COI_QUERY})'


def describe_exception(exc: Exception) -> str:
    """Do not print raw NCBI HTTP URLs: they can contain an API key."""
    if isinstance(exc, HTTPError):
        return f"HTTP error {exc.code}"
    if isinstance(exc, URLError):
        return "NCBI connection error"
    if isinstance(exc, (RuntimeError, ValueError)):
        return str(exc)  # locally generated errors omit request URLs and API keys
    return type(exc).__name__


def is_permanent_http_error(exc: Exception) -> bool:
    return isinstance(exc, HTTPError) and 400 <= exc.code < 500 and exc.code != 429


def call_with_retry(operation, label: str, tries: int, delay: float):
    """Retry failed E-utilities requests, including truncated/invalid results."""
    for attempt in range(1, tries + 1):
        try:
            return operation()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as exc:
            kind = describe_exception(exc)
            if attempt >= tries or is_permanent_http_error(exc):
                raise RuntimeError(f"{label} failed ({kind}) after {attempt} attempt(s)") from None
            pause = min(60.0, delay * 2 ** (attempt - 1))
            LOG.warning("%s failed (%s), attempt %d/%d; retrying in %.1fs",
                        label, kind, attempt, tries, pause)
            time.sleep(pause)
    raise AssertionError("retry loop unreachable")


class RequestPacer:
    """Supplement Biopython's own automatic NCBI request throttling."""

    def __init__(self, interval: float):
        self.interval = interval
        self.last_start = 0.0

    def wait(self) -> None:
        now = time.monotonic()
        elapsed = now - self.last_start
        if elapsed < self.interval:
            time.sleep(self.interval - elapsed)
        self.last_start = time.monotonic()


def search_history(genus: str, pacer: RequestPacer, tries: int):
    """ESearch with history; no retmax limit or large ID lists to truncate."""
    query = query_for_genus(genus)

    def do_search():
        pacer.wait()
        with Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0) as handle:
            return Entrez.read(handle)

    response = call_with_retry(do_search, f"ESearch for {genus}", tries, 2.0)
    count = int(response["Count"])
    if count > 0 and (not response.get("WebEnv") or not response.get("QueryKey")):
        raise RuntimeError(f"ESearch for {genus} returned no NCBI history token")
    return {
        "query": query,
        "ncbi_query_translation": str(response.get("QueryTranslation", "")),
        "count": count,
        "webenv": str(response.get("WebEnv", "")),
        "query_key": str(response.get("QueryKey", "")),
    }


def fetch_batch(search: dict, start: int, size: int,
                pacer: RequestPacer, tries: int) -> str:
    """Fetch and parse-check one GenBank batch; return original GenBank text."""
    def do_fetch():
        pacer.wait()
        with Entrez.efetch(
            db="nucleotide", rettype="gb", retmode="text", retstart=start,
            retmax=size, webenv=search["webenv"], query_key=search["query_key"],
        ) as handle:
            raw = handle.read()
        # A successful HTTP status is not proof of a complete batch.
        # SeqIO.parse catches most truncated/malformed GenBank responses.
        record_count = sum(1 for _ in SeqIO.parse(io.StringIO(raw), "genbank"))
        if record_count != size:
            raise ValueError(f"expected {size} records, received {record_count}")
        return raw

    return call_with_retry(do_fetch, f"EFetch batch starting at {start}", tries, 2.0)


def retrieve_genus(genus: str, records_dir: Path, batch_size: int,
                   pacer: RequestPacer, tries: int) -> dict:
    """Commit a genus only after all its records have been retrieved."""
    stem = genus_stem(genus)
    path = records_dir / f"{stem}.gb"
    metadata_path = records_dir / f"{stem}.json"
    part = records_dir / f"{stem}.gb.part"
    part.unlink(missing_ok=True)   # incomplete previous attempt is never reused

    search = search_history(genus, pacer, tries)
    count = search["count"]
    LOG.info("%-25s %s matching nucleotide records", genus, count)
    retrieved = 0
    with part.open("w", encoding="utf-8", newline="\n") as out:
        for start in range(0, count, batch_size):
            size = min(batch_size, count - start)
            raw = fetch_batch(search, start, size, pacer, tries)
            out.write(raw)
            if raw and not raw.endswith("\n"):
                out.write("\n")
            retrieved += size
    if retrieved != count:
        raise RuntimeError(f"Unexpected count for {genus}: {retrieved}/{count}")

    os.replace(part, path)
    metadata = {
        "genus": genus,
        "query": search["query"],
        "ncbi_query_translation": search["ncbi_query_translation"],
        "ncbi_search_count": count,
        "downloaded_records": retrieved,
        "downloaded_at_utc": utc_now(),
        "filename": path.name,
        "sha256": sha256_file(path),
    }
    save_json_atomic(metadata_path, metadata)
    return metadata


def load_checkpoint(genus: str, records_dir: Path) -> dict | None:
    """Accept a checkpoint only if file and its checksum match the receipt."""
    stem = genus_stem(genus)
    metadata_path = records_dir / f"{stem}.json"
    path = records_dir / f"{stem}.gb"
    if not (metadata_path.is_file() and path.is_file()):
        return None
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if (metadata["genus"] == genus
                and metadata["filename"] == path.name
                and metadata["query"] == query_for_genus(genus)
                and metadata["downloaded_records"] == metadata["ncbi_search_count"]
                and metadata["sha256"] == sha256_file(path)):
            return metadata
    except (KeyError, ValueError, OSError, TypeError, json.JSONDecodeError):
        pass
    LOG.warning("Invalid/incomplete checkpoint for %s; retrieving afresh", genus)
    return None


def iter_genbank_blocks(path: Path):
    """Yield original GenBank text records without rewriting any annotations."""
    lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            lines.append(line)
            if line.strip() == "//":
                yield "".join(lines)
                lines = []
    if "".join(lines).strip():
        raise ValueError(f"Incomplete GenBank record in {path.name}")


def accession_from_block(block: str) -> str:
    match = re.search(r"^VERSION\s+(\S+)", block, flags=re.MULTILINE)
    if not match:
        # Fallback for older GenBank entries lacking VERSION.
        match = re.search(r"^ACCESSION\s+(\S+)", block, flags=re.MULTILINE)
    if not match or not block.startswith("LOCUS"):
        raise ValueError("Missing accession/version or LOCUS field in a GenBank record")
    return match.group(1)


def assemble_output(genera: list[str], records_dir: Path, outdir: Path) -> tuple[list[dict], int]:
    """Merge genus files, deduplicating identical accession.version identifiers."""
    summary = []
    seen: set[str] = set()
    target = outdir / OUTPUT_NAME
    part = outdir / (OUTPUT_NAME + ".part")
    with part.open("w", encoding="utf-8", newline="\n") as combined:
        for genus in genera:
            metadata = load_checkpoint(genus, records_dir)
            if metadata is None:
                raise RuntimeError(f"Cannot assemble dataset: {genus} lacks a valid checkpoint")
            local_count = added = 0
            for block in iter_genbank_blocks(records_dir / metadata["filename"]):
                local_count += 1
                accession = accession_from_block(block)
                if accession not in seen:
                    seen.add(accession)
                    combined.write(block)
                    added += 1
            if local_count != metadata["downloaded_records"]:
                raise RuntimeError(f"Count mismatch for {genus}: "
                                   f"found {local_count} records, expected {metadata['downloaded_records']}")
            summary.append({
                "genus": genus,
                "ncbi_search_count": metadata["ncbi_search_count"],
                "downloaded_records": local_count,
                "unique_records_added": added,
                "duplicate_accessions_skipped": local_count - added,
                "retrieved_at_utc": metadata["downloaded_at_utc"],
                "ncbi_query_translation": metadata["ncbi_query_translation"],
            })
    os.replace(part, target)
    return summary, len(seen)


def save_summary(path: Path, summary: list[dict]) -> None:
    fields = ("genus", "ncbi_search_count", "downloaded_records", "unique_records_added",
              "duplicate_accessions_skipped", "retrieved_at_utc", "ncbi_query_translation")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Retrieve COI GenBank records for WoRMS polychaete genera."
    )
    parser.add_argument("--accepted", type=Path, required=True,
                        help="WoRMS accepted-species CSV containing a Genus column")
    parser.add_argument("--outdir", type=Path, default=Path("coi_retrieval"),
                        help="Output directory (default: coi_retrieval)")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"),
                        help="NCBI contact email (or set NCBI_EMAIL)")
    parser.add_argument("--extra-genera-file", type=Path,
                        help="Optional historical genera: text lines or CSV with Genus column")
    parser.add_argument("--only-genus", action="append", default=[],
                        help="Development/subset run: restrict to named genera (repeatable)")
    parser.add_argument("--batch-size", type=int, default=200,
                        help="Records per EFetch call, 1–500 (default: 200)")
    parser.add_argument("--delay", type=float, default=0.4,
                        help="Minimum seconds between calls (default: 0.4)")
    parser.add_argument("--retries", type=int, default=4,
                        help="Maximum attempts per NCBI operation (default: 4)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume a matching earlier run, checking committed genus files")
    args = parser.parse_args(argv)
    if not args.email or "@" not in args.email:
        parser.error("provide --email or set the NCBI_EMAIL environment variable")
    if args.batch_size < 1 or args.batch_size > 500:
        parser.error("--batch-size must be between 1 and 500")
    if args.retries < 1:
        parser.error("--retries must be at least 1")
    if args.delay < 0:
        parser.error("--delay must be nonnegative")
    if args.delay < (0.11 if os.environ.get("NCBI_API_KEY") else 0.34):
        parser.error("--delay is too low for NCBI's standard rate limits")
    if not args.accepted.is_file():
        parser.error(f"input CSV does not exist: {args.accepted}")
    if args.extra_genera_file and not args.extra_genera_file.is_file():
        parser.error(f"extra genera file does not exist: {args.extra_genera_file}")
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    Entrez.email = args.email
    Entrez.api_key = os.environ.get("NCBI_API_KEY") or None
    Entrez.tool = "PolychaeteCOIAudit"
    # Use our explicit error handling rather than Biopython's additional
    # implicit retry loop. Biopython still performs automatic rate-limiting.
    Entrez.max_tries = 1
    accepted_genera = load_genera(args.accepted)
    extra_genera = load_extra_genera(args.extra_genera_file) if args.extra_genera_file else []
    genera = sorted(set(accepted_genera) | set(extra_genera), key=lambda g: (g.casefold(), g))
    if args.only_genus:
        unknown = sorted(set(args.only_genus) - set(genera))
        if unknown:
            raise ValueError(f"--only-genus names not found: {', '.join(unknown)}")
        genera = [g for g in genera if g in set(args.only_genus)]
    for genus in genera:
        query_for_genus(genus)  # validate before starting any requests

    settings = {
        "script_version": VERSION,
        "accepted_csv": args.accepted.name,
        "accepted_sha256": sha256_file(args.accepted),
        "extra_genera_file": args.extra_genera_file.name if args.extra_genera_file else None,
        "extra_genera_sha256": sha256_file(args.extra_genera_file) if args.extra_genera_file else None,
        "search_template": '"{genus}"[Organism] AND ({COI_QUERY})',
        "coi_query": COI_QUERY,
        "genera": genera,
        "subset_run": bool(args.only_genus),
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    settings_file = args.outdir / "run_settings.json"
    if settings_file.exists():
        if not args.resume:
            raise RuntimeError(f"{settings_file} already exists; use --resume or a new --outdir")
        prior = json.loads(settings_file.read_text(encoding="utf-8"))
        if {key: prior.get(key) for key in settings} != settings:
            raise RuntimeError("Cannot resume: input files or search settings changed; use a new --outdir")
    else:
        if args.resume:
            raise RuntimeError("--resume requested, but no run_settings.json exists")
        if any(args.outdir.iterdir()):
            raise RuntimeError("Output directory is not empty; choose a clean --outdir")
        save_json_atomic(settings_file, {**settings, "started_at_utc": utc_now()})

    # A previous completed manifest must not survive a subsequently failed rerun.
    (args.outdir / "run_manifest.json").unlink(missing_ok=True)
    records_dir = args.outdir / "genus_records"
    records_dir.mkdir(exist_ok=True)
    pacer = RequestPacer(args.delay)
    LOG.info("Genus queries: %d (%d accepted, %d extra); subset=%s",
             len(genera), len(accepted_genera), len(extra_genera), bool(args.only_genus))
    LOG.info("Run begun/continued; API key %s", "present" if Entrez.api_key else "not provided")

    for index, genus in enumerate(genera, 1):
        existing = load_checkpoint(genus, records_dir) if args.resume else None
        if existing is not None:
            LOG.info("[%d/%d] Reusing verified checkpoint: %s (%d records)",
                     index, len(genera), genus, existing["downloaded_records"])
            continue
        LOG.info("[%d/%d] Retrieving: %s", index, len(genera), genus)
        try:
            retrieve_genus(genus, records_dir, args.batch_size, pacer, args.retries)
        except Exception as exc:
            # The original script printed an error and silently continued;
            # this version never presents an incomplete dataset as finished.
            LOG.error("Retrieval stopped at genus %s: %s", genus, describe_exception(exc))
            LOG.error("Completed genera remain saved; run again with --resume")
            return 1

    summary, total_unique = assemble_output(genera, records_dir, args.outdir)
    save_summary(args.outdir / "retrieval_summary.csv", summary)
    merged = args.outdir / OUTPUT_NAME
    manifest = {
        "completed": True,
        "completed_at_utc": utc_now(),
        "started_at_utc": json.loads(settings_file.read_text(encoding="utf-8"))["started_at_utc"],
        "script_version": VERSION,
        "python_version": sys.version.split()[0],
        "biopython_version": Bio.__version__,
        "accepted_csv": args.accepted.name,
        "accepted_sha256": settings["accepted_sha256"],
        "extra_genera_file": settings["extra_genera_file"],
        "extra_genera_sha256": settings["extra_genera_sha256"],
        "search_template": settings["search_template"],
        "coi_query": COI_QUERY,
        "genera_queried": len(genera),
        "subset_run": settings["subset_run"],
        "ncbi_hits_total_before_deduplication": sum(row["ncbi_search_count"] for row in summary),
        "downloaded_records_before_deduplication": sum(row["downloaded_records"] for row in summary),
        "duplicate_accessions_removed": sum(row["duplicate_accessions_skipped"] for row in summary),
        "unique_records_written": total_unique,
        "combined_genbank_filename": merged.name,
        "combined_genbank_sha256": sha256_file(merged),
        "metadata_note": "The search is GenBank-only; genus-based and COI title/gene annotation-dependent. "
                         "Downloaded counts reflect NCBI at retrieval time, not the publication baseline.",
    }
    save_json_atomic(args.outdir / "run_manifest.json", manifest)
    LOG.info("Completed: %d unique GenBank records saved in %s", total_unique, merged)
    LOG.info("Saved retrieval_summary.csv and run_manifest.json in %s", args.outdir)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        # No raw URL or HTTP traceback: avoids accidentally exposing credentials.
        LOG.error("%s", describe_exception(exc) if isinstance(exc, (HTTPError, URLError)) else str(exc))
        sys.exit(1)
