"""Aggregate CovSpectrum daily counts into monthly lineage frequencies.

Reads the raw CovSpectrum TSV (columns: date, pangoLineage, count), drops rows
with missing date or lineage, converts daily dates to months using the same
rounding convention as the Genbank pipeline (day >= 16 rounds up to next month),
and computes monthly lineage frequencies.
"""

import sys

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821


def date_to_month(date_str):
    """Convert a YYYY-MM-DD date string to YYYY-MM, rounding to nearest month.

    Day >= 16 rounds up to next month, matching the Genbank pipeline convention
    in merge_genbank_metadata_and_spikeseqs.py.
    """
    parts = date_str.split("-")
    if len(parts) != 3:
        raise ValueError(f"Expected YYYY-MM-DD date format, got {date_str!r}")
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    if day >= 16:
        month += 1
        if month > 12:
            month = 1
            year += 1
    return f"{year:04d}-{month:02d}"


def main():
    df = pd.read_csv(
        snakemake.input.tsv,  # noqa: F821
        sep="\t",
        dtype={"date": str, "pangoLineage": str, "count": int},
    )
    print(f"Read {len(df)} rows from CovSpectrum data")
    print(f"Columns: {list(df.columns)}")

    # Validate expected columns
    expected_cols = {"date", "pangoLineage", "count"}
    if set(df.columns) != expected_cols:
        raise ValueError(
            f"Expected columns {expected_cols}, got {set(df.columns)}"
        )

    # Drop rows with null/empty date or pangoLineage
    n_before = len(df)
    df = df[df["date"].notna() & (df["date"] != "")]
    n_dropped_date = n_before - len(df)
    print(f"Dropped {n_dropped_date} rows with missing date")

    n_before = len(df)
    df = df[df["pangoLineage"].notna() & (df["pangoLineage"] != "")]
    n_dropped_lineage = n_before - len(df)
    print(f"Dropped {n_dropped_lineage} rows with missing pangoLineage")

    # Convert dates to months
    df["month"] = df["date"].apply(date_to_month)
    df = df.rename(columns={"pangoLineage": "lineage"})

    # Group by (lineage, month) summing counts
    monthly = (
        df.groupby(["lineage", "month"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "lineage_counts"})
    )

    # Compute total counts per month
    month_totals = (
        monthly.groupby("month", as_index=False)["lineage_counts"]
        .sum()
        .rename(columns={"lineage_counts": "total_counts"})
    )
    monthly = monthly.merge(month_totals, on="month")

    # Compute monthly frequency
    monthly["monthly_frequency"] = monthly["lineage_counts"] / monthly["total_counts"]

    # Select and sort
    monthly = monthly[
        ["lineage", "month", "lineage_counts", "total_counts", "monthly_frequency"]
    ]
    monthly = monthly.sort_values(
        ["month", "monthly_frequency"],
        ascending=[True, False],
    ).reset_index(drop=True)

    monthly.to_csv(
        snakemake.output.tsv,  # noqa: F821
        sep="\t",
        index=False,
        float_format="%.3g",
    )
    print(f"Wrote {len(monthly)} rows to {snakemake.output.tsv}")  # noqa: F821

    # Write summary
    n_lineages = monthly["lineage"].nunique()
    n_months = monthly["month"].nunique()
    total_counts = monthly["lineage_counts"].sum()
    min_month = monthly["month"].min()
    max_month = monthly["month"].max()

    summary_lines = [
        "CovSpectrum frequency parsing summary",
        "======================================",
        f"Total lineages: {n_lineages}",
        f"Total months: {n_months}",
        f"Date range: {min_month} to {max_month}",
        f"Total sequence counts: {total_counts}",
        f"Rows dropped (missing date): {n_dropped_date}",
        f"Rows dropped (missing lineage): {n_dropped_lineage}",
    ]
    summary_text = "\n".join(summary_lines) + "\n"
    with open(snakemake.output.summary, "w") as f:  # noqa: F821
        f.write(summary_text)
    print(summary_text)


main()
