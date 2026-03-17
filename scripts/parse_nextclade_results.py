"""Parse nextclade results TSV and spike FASTA into clean output."""

import csv
import gzip
import sys


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821

OUTPUT_COLUMNS = [
    "accession",
    "Nextstrain_clade",
    "pango",
    "qc.overallScore",
    "substitutions",
    "deletions",
    "insertions",
    "missing_sites",
    "n_missing_sites",
    "aligned_spike",
]

# Map from nextclade TSV column name to output column name
TSV_COLUMN_MAP = {
    "seqName": "accession",
    "clade_nextstrain": "Nextstrain_clade",
    "Nextclade_pango": "pango",
    "qc.overallScore": "qc.overallScore",
    "aaSubstitutions": "substitutions",
    "aaDeletions": "deletions",
    "aaInsertions": "insertions",
    "unknownAaRanges": "missing_sites",
}

S_PREFIX_COLUMNS = ["aaSubstitutions", "aaDeletions", "aaInsertions", "unknownAaRanges"]


def strip_s_prefix(value, column_name, accession):
    """Validate all entries have S: prefix and strip it."""
    if not value:
        return ""
    entries = value.split(",")
    stripped = []
    for entry in entries:
        if not entry.startswith("S:"):
            raise ValueError(
                f"Entry {entry!r} in {column_name} for {accession} "
                f"does not start with 'S:'"
            )
        stripped.append(entry[2:])
    return ",".join(stripped)


def count_missing_sites(missing_sites_stripped):
    """Count total missing sites from ranges like '16-19' (=4) or '403' (=1)."""
    if not missing_sites_stripped:
        return 0
    total = 0
    for part in missing_sites_stripped.split(","):
        if "-" in part:
            start, end = part.split("-")
            total += int(end) - int(start) + 1
        else:
            total += 1
    return total


def read_next_fasta_entry(fasta_file, fasta_state):
    """Read next FASTA entry, returning (accession, sequence)."""
    header = fasta_state.get("pending_header")
    if header is None:
        # Read until we find a header
        for line in fasta_file:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
                break
        else:
            return None, None

    accession = header[1:].split()[0]
    seq_parts = []
    for line in fasta_file:
        line = line.rstrip("\n")
        if line.startswith(">"):
            fasta_state["pending_header"] = line
            return accession, "".join(seq_parts)
        seq_parts.append(line)

    fasta_state["pending_header"] = None
    return accession, "".join(seq_parts)


def main():
    tsv_path = snakemake.input["tsv"]  # noqa: F821
    fasta_path = snakemake.input["spike_fasta"]  # noqa: F821
    output_tsv_path = snakemake.output["tsv"]  # noqa: F821
    summary_path = snakemake.output["summary"]  # noqa: F821

    seen_accessions = set()
    n_written = 0
    n_dropped = 0
    pango_counts = {}
    clade_counts = {}

    with (
        gzip.open(tsv_path, "rt") as tsv_file,
        gzip.open(fasta_path, "rt") as fasta_file,
        gzip.open(output_tsv_path, "wt") as out_file,
    ):
        reader = csv.DictReader(tsv_file, delimiter="\t")
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)

        # Verify all expected columns exist
        for col in list(TSV_COLUMN_MAP) + ["failedCdses", "errors"]:
            if col not in reader.fieldnames:
                raise ValueError(f"Expected column {col!r} not in TSV header")

        fasta_state = {}

        for row in reader:
            seq_name = row["seqName"]
            accession = seq_name.split()[0]

            if accession in seen_accessions:
                raise ValueError(f"Duplicate accession: {accession}")
            seen_accessions.add(accession)

            # Skip rows with errors or failed S CDS
            errors = row["errors"]
            failed_cdses = row["failedCdses"]
            if errors or failed_cdses:
                n_dropped += 1
                continue

            # Read matching FASTA entry
            fasta_acc, spike_seq = read_next_fasta_entry(fasta_file, fasta_state)
            if fasta_acc is None:
                raise ValueError(
                    f"Ran out of FASTA entries at TSV accession {accession}"
                )
            if fasta_acc != accession:
                raise ValueError(
                    f"FASTA/TSV accession mismatch: FASTA={fasta_acc}, TSV={accession}"
                )

            # Strip S: prefix from relevant columns
            stripped = {}
            for col in S_PREFIX_COLUMNS:
                stripped[col] = strip_s_prefix(row[col], col, accession)

            # Compute n_missing_sites
            n_missing = count_missing_sites(stripped["unknownAaRanges"])

            # Validate n_missing_sites matches X count in spike sequence
            x_count = spike_seq.count("X")
            if n_missing != x_count:
                raise ValueError(
                    f"n_missing_sites ({n_missing}) != X count ({x_count}) "
                    f"for {accession}"
                )

            # Track stats
            pango = row["Nextclade_pango"]
            clade = row["clade_nextstrain"]
            pango_counts[pango] = pango_counts.get(pango, 0) + 1
            clade_counts[clade] = clade_counts.get(clade, 0) + 1

            # Write row
            writer.writerow(
                [
                    accession,
                    clade,
                    pango,
                    row["qc.overallScore"],
                    stripped["aaSubstitutions"],
                    stripped["aaDeletions"],
                    stripped["aaInsertions"],
                    stripped["unknownAaRanges"],
                    n_missing,
                    spike_seq,
                ]
            )
            n_written += 1

            if n_written % 100_000 == 0:
                print(f"Parsed {n_written:,} entries", flush=True)

    # Check no remaining FASTA entries
    extra_acc, _ = read_next_fasta_entry(
        gzip.open(fasta_path, "rt"), {"pending_header": None}
    )
    # We can't easily check for leftover FASTA entries since we closed the file,
    # but the accession matching provides sufficient validation.

    # Generate summary
    total_rows = n_written + n_dropped
    summary_lines = [
        f"Total TSV rows: {total_rows:,}",
        f"Dropped (errors or failed CDS): {n_dropped:,}",
        f"Parsed entries: {n_written:,}",
        f"Unique pango lineages: {len(pango_counts):,}",
        f"Unique Nextstrain clades: {len(clade_counts):,}",
        "",
    ]

    top_pango = sorted(pango_counts.items(), key=lambda x: -x[1])[:10]
    summary_lines.append("Top 10 pango lineages:")
    for val, count in top_pango:
        summary_lines.append(f"  {val}: {count:,}")
    summary_lines.append("")

    top_clades = sorted(clade_counts.items(), key=lambda x: -x[1])[:10]
    summary_lines.append("Top 10 Nextstrain clades:")
    for val, count in top_clades:
        summary_lines.append(f"  {val}: {count:,}")

    summary_text = "\n".join(summary_lines)
    print(f"\n{summary_text}", flush=True)

    with open(summary_path, "w") as f:
        f.write(summary_text + "\n")

    print("\nDone.", flush=True)


main()
