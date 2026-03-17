"""Merge Genbank and CovSpectrum monthly lineage frequencies.

Performs a full outer join on (lineage, month), fills missing counts with 0,
and produces parallel columns for each source.
"""

import sys

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821


def main():
    genbank = pd.read_csv(snakemake.input.genbank, sep="\t")  # noqa: F821
    covspectrum = pd.read_csv(snakemake.input.covspectrum, sep="\t")  # noqa: F821
    print(f"Genbank: {len(genbank)} rows, {genbank['lineage'].nunique()} lineages")
    print(f"CovSpectrum: {len(covspectrum)} rows, {covspectrum['lineage'].nunique()} lineages")

    # Build total_counts lookup per month for each source
    gb_month_totals = (
        genbank.groupby("month")["total_counts"].first().to_dict()
    )
    cs_month_totals = (
        covspectrum.groupby("month")["total_counts"].first().to_dict()
    )

    # Full outer join on (lineage, month)
    merged = genbank.merge(
        covspectrum,
        on=["lineage", "month"],
        how="outer",
        suffixes=("_genbank", "_covspectrum"),
    )
    print(f"Merged: {len(merged)} rows")

    # Fill missing lineage_counts with 0
    merged["lineage_counts_genbank"] = merged["lineage_counts_genbank"].fillna(0).astype(int)
    merged["lineage_counts_covspectrum"] = merged["lineage_counts_covspectrum"].fillna(0).astype(int)

    # Fill total_counts from source lookups (a month always exists in at least one source)
    merged["total_counts_genbank"] = merged["month"].map(gb_month_totals).fillna(0).astype(int)
    merged["total_counts_covspectrum"] = merged["month"].map(cs_month_totals).fillna(0).astype(int)

    # Compute monthly frequencies
    merged["monthly_frequency_genbank"] = merged.apply(
        lambda r: r["lineage_counts_genbank"] / r["total_counts_genbank"]
        if r["total_counts_genbank"] > 0
        else 0.0,
        axis=1,
    )
    merged["monthly_frequency_covspectrum"] = merged.apply(
        lambda r: r["lineage_counts_covspectrum"] / r["total_counts_covspectrum"]
        if r["total_counts_covspectrum"] > 0
        else 0.0,
        axis=1,
    )

    # Select output columns
    result = merged[
        [
            "lineage",
            "month",
            "lineage_counts_covspectrum",
            "lineage_counts_genbank",
            "total_counts_covspectrum",
            "total_counts_genbank",
            "monthly_frequency_covspectrum",
            "monthly_frequency_genbank",
        ]
    ]

    # Sort by month ascending, then max of the two frequencies descending
    result = result.assign(
        _max_freq=result[["monthly_frequency_covspectrum", "monthly_frequency_genbank"]].max(axis=1)
    )
    result = result.sort_values(
        ["month", "_max_freq"],
        ascending=[True, False],
    ).drop(columns="_max_freq").reset_index(drop=True)

    result.to_csv(
        snakemake.output.tsv,  # noqa: F821
        sep="\t",
        index=False,
        float_format="%.3g",
    )
    print(f"Wrote {len(result)} rows to {snakemake.output.tsv}")  # noqa: F821

    # Summary stats
    n_lineages = result["lineage"].nunique()
    n_months = result["month"].nunique()
    min_month = result["month"].min()
    max_month = result["month"].max()
    n_both = ((result["lineage_counts_genbank"] > 0) & (result["lineage_counts_covspectrum"] > 0)).sum()
    n_gb_only = ((result["lineage_counts_genbank"] > 0) & (result["lineage_counts_covspectrum"] == 0)).sum()
    n_cs_only = ((result["lineage_counts_genbank"] == 0) & (result["lineage_counts_covspectrum"] > 0)).sum()

    print("\nSummary:")
    print(f"  Total lineages: {n_lineages}")
    print(f"  Total months: {n_months}")
    print(f"  Date range: {min_month} to {max_month}")
    print(f"  Rows with data in both sources: {n_both}")
    print(f"  Rows with Genbank only: {n_gb_only}")
    print(f"  Rows with CovSpectrum only: {n_cs_only}")

    # Write summary (same content as log)
    sys.stdout.flush()
    with open(snakemake.log[0]) as log_f:  # noqa: F821
        summary_text = log_f.read()
    with open(snakemake.output.summary, "w") as f:  # noqa: F821
        f.write(summary_text)


main()
