"""Merge genbank metadata, nextclade results, and UShER metadata into one file."""

import sys

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821

# ---------------------------------------------------------------------------
# Column selection and renaming for each input source
# ---------------------------------------------------------------------------

GENBANK_COLUMNS = {
    "accession": "accession",
    "date": "date_genbank",
    "geographic_region": "geographic_region",
    "pango": "pango_lineage_genbank",
}

NEXTCLADE_COLUMNS = {
    "accession": "accession",
    "Nextstrain_clade": "nextstrain_clade_nextclade",
    "pango": "pango_lineage_nextclade",
    "qc.overallScore": "nextclade_qc_score",
    "substitutions": "substitutions",
    "deletions": "deletions",
    "insertions": "insertions",
    "missing_sites": "missing_sites",
    "n_missing_sites": "n_missing_sites",
    "aligned_spike": "aligned_spike",
}

USHER_COLUMNS = {
    "genbank_accession": "accession",
    "date": "date_usher",
    "country": "country",
    "pango_lineage_usher": "pango_lineage_usher",
    "Nextstrain_clade_usher": "nextstrain_clade_usher",
}

# Desired column order after merge (date_month columns inserted in add_date_month)
COLUMN_ORDER = [
    "accession",
    # dates
    "date_genbank",
    "date_month_genbank",
    "date_usher",
    "date_month_usher",
    # geography
    "geographic_region",
    "country",
    # pango lineages
    "pango_lineage_genbank",
    "pango_lineage_nextclade",
    "pango_lineage_usher",
    # nextstrain clades
    "nextstrain_clade_nextclade",
    "nextstrain_clade_usher",
    # sequence quality
    "nextclade_qc_score",
    "n_missing_sites",
    # mutations and missing sites
    "substitutions",
    "deletions",
    "insertions",
    "missing_sites",
    # sequence last
    "aligned_spike",
]

# Columns to skip in per-column summaries
SKIP_SUMMARY_COLUMNS = {"substitutions", "deletions", "insertions", "missing_sites"}

# Columns that get only non-null vs null counts
ALIGNED_SPIKE_COLUMNS = {"aligned_spike"}

# Numeric / date columns that get percentile summaries
DATE_COLUMNS = {"date_genbank", "date_usher"}
DATE_MONTH_COLUMNS = {"date_month_genbank", "date_month_usher"}
NUMERIC_COLUMNS = {"nextclade_qc_score", "n_missing_sites"}

PERCENTILES = [0.01, 0.1, 1, 10, 50, 90, 99, 99.9, 99.99]

# Cross-source comparison pairs
COMPARISON_PAIRS = [
    ("pango_lineage_genbank", "pango_lineage_nextclade"),
    ("pango_lineage_genbank", "pango_lineage_usher"),
    ("pango_lineage_nextclade", "pango_lineage_usher"),
    ("nextstrain_clade_nextclade", "nextstrain_clade_usher"),
    ("date_genbank", "date_usher"),
    ("date_month_genbank", "date_month_usher"),
]

VALID_PANGO_PATTERN = r"[A-Z\d.]+"
PANGO_COLUMNS = {
    "pango_lineage_genbank": "genbank",
    "pango_lineage_nextclade": "nextclade",
    "pango_lineage_usher": "usher",
}


def nullify_invalid_pango(df, pango_columns, valid_pattern):
    """Set invalid Pango lineage names to None and report what was cleaned.

    Returns a list of summary lines describing what was nullified.
    """
    lines = []
    for col_name, source_name in pango_columns.items():
        col = df[col_name]
        matches = col.str.fullmatch(valid_pattern)
        invalid = col.notna() & ~matches.astype(bool)
        n_invalid = invalid.sum()
        if n_invalid:
            top_invalid = col[invalid].value_counts().head(10)
            msg = f"Nullifying {n_invalid:,} invalid {source_name} Pango lineage values:"
            print(msg, flush=True)
            lines.append(msg)
            for val, count in top_invalid.items():
                detail = f"  {val}: {count:,}"
                print(detail, flush=True)
                lines.append(detail)
            df.loc[invalid, col_name] = None
        else:
            msg = f"No invalid {source_name} Pango lineage values found."
            print(msg, flush=True)
            lines.append(msg)
    return lines


def read_and_rename(path, column_map, source_name):
    """Read a TSV, select and rename columns, validate unique accessions."""
    print(f"Reading {source_name} from {path}...", flush=True)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(df):,} rows, columns: {list(df.columns)}", flush=True)

    missing = set(column_map) - set(df.columns)
    if missing:
        raise ValueError(f"{source_name}: missing expected columns {missing}")

    df = df[list(column_map)].rename(columns=column_map)
    print(f"  Renamed columns: {list(df.columns)}", flush=True)

    # Drop rows with missing accession
    n_missing_acc = df["accession"].isna().sum()
    if n_missing_acc:
        print(f"  Dropping {n_missing_acc:,} rows with missing accession", flush=True)
        df = df.dropna(subset=["accession"])

    dups = df["accession"].duplicated()
    if dups.any():
        examples = df.loc[dups, "accession"].head(5).tolist()
        raise ValueError(f"{source_name}: duplicate accessions, e.g. {examples}")

    print(f"  {len(df):,} rows after column selection", flush=True)
    return df


def date_to_month(date_series):
    """Convert date strings to YYYY-MM, rounding to nearest month.

    For YYYY-MM-DD dates, day >= 16 rounds up to next month.
    For YYYY-MM dates (no day), keeps as-is.
    For YYYY-only dates, returns NaN.
    """

    def _to_month(val):
        if pd.isna(val) or len(val) < 7 or val[4] != "-":
            return None
        if len(val) == 7:
            # Already YYYY-MM
            return val
        # YYYY-MM-DD: round to nearest month
        year = int(val[:4])
        month = int(val[5:7])
        day = int(val[8:10])
        if day >= 16:
            month += 1
            if month > 12:
                month = 1
                year += 1
        return f"{year:04d}-{month:02d}"

    return date_series.apply(_to_month)


def pairwise_comparison(merged, col_a, col_b):
    """Compare two columns: both non-null equal, both non-null different, etc."""
    a_null = merged[col_a].isna()
    b_null = merged[col_b].isna()
    both_non_null = ~a_null & ~b_null
    both_null = a_null & b_null
    one_null = a_null ^ b_null
    equal = both_non_null & (merged[col_a] == merged[col_b])
    different = both_non_null & (merged[col_a] != merged[col_b])
    return {
        "both_non_null_equal": int(equal.sum()),
        "both_non_null_different": int(different.sum()),
        "one_null": int(one_null.sum()),
        "both_null": int(both_null.sum()),
    }


def generate_summary(genbank_df, nextclade_df, usher_df, merged, pango_nullify_lines):
    """Generate a detailed summary of the merge."""
    lines = []

    # --- Invalid Pango lineage cleanup ---
    lines.append("=" * 70)
    lines.append("INVALID PANGO LINEAGE CLEANUP")
    lines.append("=" * 70)
    lines.extend(pango_nullify_lines)
    lines.append("")

    # --- Input file descriptions ---
    lines.append("=" * 70)
    lines.append("INPUT FILE DESCRIPTIONS")
    lines.append("=" * 70)
    lines.append(f"Genbank metadata: {len(genbank_df):,} rows")
    lines.append(f"Nextclade results: {len(nextclade_df):,} rows")
    lines.append(f"UShER metadata: {len(usher_df):,} rows")
    lines.append(f"Merged output: {len(merged):,} rows")
    lines.append("")

    # --- Accession overlap ---
    lines.append("=" * 70)
    lines.append("ACCESSION OVERLAP")
    lines.append("=" * 70)
    sources = {
        "genbank": set(genbank_df["accession"]),
        "nextclade": set(nextclade_df["accession"]),
        "usher": set(usher_df["accession"]),
    }
    source_names = list(sources)
    for i in range(len(source_names)):
        for j in range(i + 1, len(source_names)):
            name_a = source_names[i]
            name_b = source_names[j]
            set_a = sources[name_a]
            set_b = sources[name_b]
            shared = set_a & set_b
            only_a = set_a - set_b
            only_b = set_b - set_a
            lines.append(f"\n{name_a} vs {name_b}:")
            lines.append(f"  Shared accessions: {len(shared):,}")
            lines.append(f"  Only in {name_a}: {len(only_a):,}")
            if only_a:
                examples = sorted(only_a)[:5]
                lines.append(f"    Examples: {', '.join(examples)}")
            lines.append(f"  Only in {name_b}: {len(only_b):,}")
            if only_b:
                examples = sorted(only_b)[:5]
                lines.append(f"    Examples: {', '.join(examples)}")
    lines.append("")

    # --- Column summaries ---
    lines.append("=" * 70)
    lines.append("COLUMN SUMMARIES")
    lines.append("=" * 70)
    for col in merged.columns:
        if col == "accession" or col in SKIP_SUMMARY_COLUMNS:
            continue

        lines.append(f"\n--- {col} ---")
        non_null = merged[col].notna().sum()
        null_count = merged[col].isna().sum()

        if col in ALIGNED_SPIKE_COLUMNS:
            lines.append(f"  Non-null: {non_null:,}")
            lines.append(f"  Null: {null_count:,}")
            continue

        if col in DATE_COLUMNS:
            lines.append(f"  Non-null: {non_null:,}, Null: {null_count:,}")
            valid_dates = pd.to_datetime(
                merged[col], format="mixed", errors="coerce"
            ).dropna()
            if len(valid_dates) > 0:
                valid_dates_sorted = valid_dates.sort_values()
                lines.append(
                    f"  Range: {valid_dates_sorted.iloc[0].date()} "
                    f"to {valid_dates_sorted.iloc[-1].date()}"
                )
                lines.append(
                    f"  Median: {valid_dates_sorted.iloc[len(valid_dates_sorted) // 2].date()}"
                )
                for p in PERCENTILES:
                    idx = min(
                        int(p / 100 * len(valid_dates_sorted)),
                        len(valid_dates_sorted) - 1,
                    )
                    lines.append(
                        f"  {p}th percentile: {valid_dates_sorted.iloc[idx].date()}"
                    )
            continue

        if col in DATE_MONTH_COLUMNS:
            # Categorical-style summary for YYYY-MM values
            lines.append(f"  Non-null: {non_null:,}, Null: {null_count:,}")
            if non_null > 0:
                sorted_vals = merged[col].dropna().sort_values()
                lines.append(
                    f"  Range: {sorted_vals.iloc[0]} to {sorted_vals.iloc[-1]}"
                )
                lines.append(f"  Median: {sorted_vals.iloc[len(sorted_vals) // 2]}")
                for p in PERCENTILES:
                    idx = min(int(p / 100 * len(sorted_vals)), len(sorted_vals) - 1)
                    lines.append(f"  {p}th percentile: {sorted_vals.iloc[idx]}")
            continue

        if col in NUMERIC_COLUMNS:
            lines.append(f"  Non-null: {non_null:,}, Null: {null_count:,}")
            numeric_vals = pd.to_numeric(merged[col], errors="coerce").dropna()
            if len(numeric_vals) > 0:
                lines.append(
                    f"  Range: {numeric_vals.min():.4g} to {numeric_vals.max():.4g}"
                )
                lines.append(f"  Median: {numeric_vals.median():.4g}")
                for p in PERCENTILES:
                    lines.append(
                        f"  {p}th percentile: {numeric_vals.quantile(p / 100):.4g}"
                    )
            continue

        # Categorical column: unique count, null count, top 5
        unique_count = merged[col].nunique()
        lines.append(f"  Unique values: {unique_count:,}")
        lines.append(f"  Null: {null_count:,}")
        top_values = merged[col].value_counts().head(5)
        lines.append("  Top 5 values:")
        for val, count in top_values.items():
            lines.append(f"    {val}: {count:,}")

    lines.append("")

    # --- Cross-source comparisons ---
    lines.append("=" * 70)
    lines.append("CROSS-SOURCE COMPARISONS")
    lines.append("=" * 70)
    for col_a, col_b in COMPARISON_PAIRS:
        lines.append(f"\n{col_a} vs {col_b}:")
        result = pairwise_comparison(merged, col_a, col_b)
        lines.append(f"  Both non-null and equal: {result['both_non_null_equal']:,}")
        lines.append(
            f"  Both non-null and different: {result['both_non_null_different']:,}"
        )
        lines.append(f"  One null: {result['one_null']:,}")
        lines.append(f"  Both null: {result['both_null']:,}")

    return "\n".join(lines)


def main():
    genbank_path = snakemake.input["genbank"]  # noqa: F821
    nextclade_path = snakemake.input["nextclade"]  # noqa: F821
    usher_path = snakemake.input["usher"]  # noqa: F821
    output_tsv_path = snakemake.output["tsv"]  # noqa: F821
    summary_path = snakemake.output["summary"]  # noqa: F821

    # Step 1: Read and rename
    genbank_df = read_and_rename(genbank_path, GENBANK_COLUMNS, "Genbank")
    nextclade_df = read_and_rename(nextclade_path, NEXTCLADE_COLUMNS, "Nextclade")
    usher_df = read_and_rename(usher_path, USHER_COLUMNS, "UShER")

    # Validate no column name collisions after renaming (except accession)
    all_renamed = (
        list(GENBANK_COLUMNS.values())
        + list(NEXTCLADE_COLUMNS.values())
        + list(USHER_COLUMNS.values())
    )
    non_accession = [c for c in all_renamed if c != "accession"]
    if len(non_accession) != len(set(non_accession)):
        raise ValueError(f"Column name collision after renaming: {non_accession}")

    # Step 2: Three-way outer merge
    print("\nMerging on accession (outer join)...", flush=True)
    merged = genbank_df.merge(nextclade_df, on="accession", how="outer")
    merged = merged.merge(usher_df, on="accession", how="outer")
    print(f"Merged: {len(merged):,} rows, {len(merged.columns)} columns", flush=True)

    # Validate no unexpected duplicate columns
    if merged.columns.duplicated().any():
        dups = merged.columns[merged.columns.duplicated()].tolist()
        raise ValueError(f"Unexpected duplicate columns after merge: {dups}")

    # Step 2b: Add date_month columns
    print("Adding date_month columns...", flush=True)
    merged["date_month_genbank"] = date_to_month(merged["date_genbank"])
    merged["date_month_usher"] = date_to_month(merged["date_usher"])

    # Step 2b2: Nullify invalid Pango lineage names
    pango_nullify_lines = nullify_invalid_pango(merged, PANGO_COLUMNS, VALID_PANGO_PATTERN)

    # Step 2c: Reorder columns
    merged = merged[COLUMN_ORDER]
    print(f"Column order: {list(merged.columns)}", flush=True)

    # Step 3: Write output
    print(f"\nWriting merged TSV to {output_tsv_path}...", flush=True)
    merged.to_csv(output_tsv_path, sep="\t", index=False, compression="gzip")
    print(f"Wrote {len(merged):,} rows", flush=True)

    # Step 4: Generate summary
    print("\nGenerating summary...", flush=True)
    summary_text = generate_summary(genbank_df, nextclade_df, usher_df, merged, pango_nullify_lines)
    print(f"\n{summary_text}", flush=True)

    with open(summary_path, "w") as f:
        f.write(summary_text + "\n")

    print("\nDone.", flush=True)


main()
