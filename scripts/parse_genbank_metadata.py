"""Parse GenBank metadata JSONL and cross-check with genomic FASTA accessions."""

import csv
import datetime
import gzip
import json
import math
import re
import sys


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821


_DATE_RE = re.compile(r"^\d{4}(-\d{2}(-\d{2})?)?$")

# Single source of truth for metadata fields: (key, display_name, tsv_column)
METADATA_FIELDS = [
    ("date", "date", "date"),
    ("region", "geographic region", "geographic_region"),
    ("location", "geographic location", "geographic_location"),
    ("pango", "pango lineage", "pango"),
]


def parse_jsonl_metadata(jsonl_path):
    """Parse metadata from JSONL, returning dict keyed by accession."""
    print(f"Parsing JSONL metadata from {jsonl_path}...", flush=True)
    metadata = {}
    with gzip.open(jsonl_path, "rt") as f:
        for i, line in enumerate(f, 1):
            record = json.loads(line)
            accession = record["accession"]
            if accession in metadata:
                raise ValueError(f"Duplicate accession in JSONL: {accession}")

            date = record.get("isolate", {}).get("collectionDate")
            if date is not None and not _DATE_RE.match(date):
                raise ValueError(f"Unexpected date format for {accession}: {date!r}")

            geographic_region = record.get("location", {}).get("geographicRegion")
            geographic_location = record.get("location", {}).get("geographicLocation")
            pango = record.get("virus", {}).get("pangolinClassification")

            metadata[accession] = {
                "date": date,
                "region": geographic_region,
                "location": geographic_location,
                "pango": pango,
            }

            if i % 100_000 == 0:
                print(f"Parsed {i:,} JSONL records", flush=True)

    print(f"Finished parsing {i:,} JSONL records total", flush=True)
    return metadata


def parse_fasta_accessions(fasta_path):
    """Parse accessions from FASTA headers, returning ordered list."""
    print(f"Parsing FASTA accessions from {fasta_path}...", flush=True)
    accessions = []
    seen = set()
    with gzip.open(fasta_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                accession = line[1:].split()[0]
                if accession in seen:
                    raise ValueError(f"Duplicate accession in FASTA: {accession}")
                seen.add(accession)
                accessions.append(accession)
                if len(accessions) % 100_000 == 0:
                    print(f"Parsed {len(accessions):,} FASTA accessions", flush=True)

    print(f"Finished parsing {len(accessions):,} FASTA accessions total", flush=True)
    return accessions


def compute_summary(metadata, fasta_accessions):
    """Compute summary statistics."""
    lines = []

    jsonl_accessions = set(metadata)
    fasta_set = set(fasta_accessions)
    jsonl_only = jsonl_accessions - fasta_set
    fasta_only = fasta_set - jsonl_accessions

    # Collect field values and null counts
    values_by_key = {key: [] for key, _, _ in METADATA_FIELDS}
    null_counts = {key: 0 for key, _, _ in METADATA_FIELDS}

    for record in metadata.values():
        for key, _, _ in METADATA_FIELDS:
            if record[key] is None:
                null_counts[key] += 1
            else:
                values_by_key[key].append(record[key])

    total = len(fasta_accessions)
    lines.append(f"Total sequences (FASTA): {total:,}")
    lines.append(f"Total JSONL records: {len(metadata):,}")
    for key, display, _ in METADATA_FIELDS:
        lines.append(f"Unique {display}s: {len(set(values_by_key[key])):,}")

    dates = values_by_key["date"]
    if dates:
        sorted_dates = sorted(dates)
        lines.append(f"Date range: {sorted_dates[0]} to {sorted_dates[-1]}")

    lines.append("")

    # Top 10 for each field (excluding date)
    for key, display, _ in METADATA_FIELDS:
        if key == "date":
            continue
        values = values_by_key[key]
        counts = {}
        for v in values:
            counts[v] = counts.get(v, 0) + 1
        top = sorted(counts.items(), key=lambda x: -x[1])[:10]
        lines.append(f"Top 10 {display} values:")
        for val, count in top:
            lines.append(f"  {val}: {count:,}")
        lines.append(f"  (null: {null_counts[key]:,})")
        lines.append("")

    lines.append(
        "Null counts: "
        + ", ".join(f"{key}={null_counts[key]:,}" for key, _, _ in METADATA_FIELDS)
    )
    lines.append("")

    # Date percentiles (YYYY-MM-DD only)
    full_dates = sorted(datetime.date.fromisoformat(d) for d in dates if len(d) == 10)
    if full_dates:
        lines.append(
            f"Date percentiles (from {len(full_dates):,} fully-resolved dates):"
        )
        for p in [1, 10, 25, 50, 75, 90, 99]:
            idx = min(math.floor(p / 100 * len(full_dates)), len(full_dates) - 1)
            lines.append(f"  {p}th: {full_dates[idx]}")
        lines.append("")

    # Cross-check
    lines.append("Accession cross-check:")
    lines.append(f"  In both JSONL and FASTA: {len(jsonl_accessions & fasta_set):,}")
    lines.append(f"  JSONL-only: {len(jsonl_only):,}")
    lines.append(f"  FASTA-only: {len(fasta_only):,}")
    if jsonl_only:
        sample = sorted(jsonl_only)[:5]
        lines.append(f"  Sample JSONL-only: {', '.join(sample)}")
    if fasta_only:
        sample = sorted(fasta_only)[:5]
        lines.append(f"  Sample FASTA-only: {', '.join(sample)}")
    failed_cross_check = fasta_only or jsonl_only

    return "\n".join(lines), failed_cross_check


def main():
    metadata = parse_jsonl_metadata(snakemake.input["data_report"])  # noqa: F821
    fasta_accessions = parse_fasta_accessions(
        snakemake.input["genomic_fasta"]  # noqa: F821
    )

    # Write summary
    summary_text, failed_cross_check = compute_summary(metadata, fasta_accessions)
    print(f"\n{summary_text}", flush=True)
    if failed_cross_check:
        raise ValueError("Some accessions found in only JSONL or only FASTA")
    with open(snakemake.output["summary"], "w") as f:  # noqa: F821
        f.write(summary_text + "\n")

    # Write TSV in FASTA order
    columns = ["accession"] + [col for _, _, col in METADATA_FIELDS]
    print("\nWriting output TSV...", flush=True)
    with gzip.open(snakemake.output["metadata"], "wt") as f:  # noqa: F821
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for accession in fasta_accessions:
            if accession in metadata:
                record = metadata[accession]
                writer.writerow(
                    [accession] + [record[key] for key, _, _ in METADATA_FIELDS]
                )
            else:
                writer.writerow([accession] + [None] * len(METADATA_FIELDS))

    print("Done.", flush=True)


main()
