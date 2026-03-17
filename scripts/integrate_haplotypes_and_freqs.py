"""Combine spike haplotypes with frequency summaries; produce collapsed version.

Reads the integrated spike haplotypes and integrated monthly frequencies,
computes per-lineage frequency summary stats, and produces two outputs:
- uncollapsed: one row per lineage with frequency summaries
- collapsed: lineages with identical complete_spike grouped together,
  frequencies re-aggregated, representative lineage chosen by max frequency
"""

import math
import sys

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821


def _median_month(group_sorted, counts_col):
    """Find the month containing the median sequence by count.

    Parameters
    ----------
    group_sorted : DataFrame
        Rows for one lineage, sorted by month ascending.
    counts_col : str or Series
        Column name in group_sorted, or a Series aligned with group_sorted.

    Returns
    -------
    str or None
        The month (YYYY-MM) containing the median sequence, or None if
        total counts is 0 (lineage absent from this source).
    """
    if isinstance(counts_col, str):
        counts = group_sorted[counts_col]
    else:
        counts = counts_col
    total = counts.sum()
    if total <= 0:
        return None
    cum = counts.cumsum()
    return group_sorted.loc[cum >= math.ceil(total / 2), "month"].iloc[0]


def compute_freq_summaries(freqs):
    """Compute per-lineage frequency summary statistics.

    Parameters
    ----------
    freqs : DataFrame
        Must have columns: lineage, month, lineage_counts_covspectrum,
        lineage_counts_genbank, monthly_frequency_covspectrum,
        monthly_frequency_genbank.

    Returns
    -------
    DataFrame with one row per lineage and columns:
        max_monthly_frequency_covspectrum, max_monthly_frequency_genbank,
        median_month_covspectrum, median_month_genbank,
        total_counts_covspectrum, total_counts_genbank,
        max_monthly_frequency, median_month, total_counts
    """
    summaries = []
    for lineage, group in freqs.groupby("lineage"):
        max_freq_cs = group["monthly_frequency_covspectrum"].max()
        max_freq_gb = group["monthly_frequency_genbank"].max()

        total_cs = group["lineage_counts_covspectrum"].sum()
        total_gb = group["lineage_counts_genbank"].sum()

        group_sorted = group.sort_values("month")

        median_month_cs = _median_month(group_sorted, "lineage_counts_covspectrum")
        median_month_gb = _median_month(group_sorted, "lineage_counts_genbank")

        # Combined median_month: use max(covspectrum, genbank) counts per month
        # to avoid double-counting (CovSpectrum includes GenBank/ENA/DDBJ data)
        combined_counts = group_sorted[
            ["lineage_counts_covspectrum", "lineage_counts_genbank"]
        ].max(axis=1)
        median_month_combined = _median_month(group_sorted, combined_counts)
        if median_month_combined is None:
            raise ValueError(
                f"Lineage {lineage} has zero counts in both sources; "
                f"expected > 0 after upstream filtering"
            )

        summaries.append(
            {
                "lineage": lineage,
                "max_monthly_frequency_covspectrum": max_freq_cs,
                "max_monthly_frequency_genbank": max_freq_gb,
                "median_month_covspectrum": median_month_cs,
                "median_month_genbank": median_month_gb,
                "total_counts_covspectrum": total_cs,
                "total_counts_genbank": total_gb,
                "max_monthly_frequency": max(max_freq_cs, max_freq_gb),
                "median_month": median_month_combined,
                "total_counts": max(total_cs, total_gb),
            }
        )

    return pd.DataFrame(summaries)


def join_haplotypes_and_summaries(haplotypes, freq_summaries, haplotype_cols):
    """Left join haplotypes with frequency summaries on lineage.

    Returns DataFrame with the final column order (without equivalent_lineages).
    """
    joined = haplotypes.merge(freq_summaries, on="lineage", how="left")

    freq_cols = [
        "max_monthly_frequency",
        "max_monthly_frequency_covspectrum",
        "max_monthly_frequency_genbank",
        "median_month",
        "median_month_covspectrum",
        "median_month_genbank",
        "total_counts",
        "total_counts_covspectrum",
        "total_counts_genbank",
    ]
    col_order = ["lineage", *freq_cols, *haplotype_cols]
    return joined[col_order]


def sort_output(df):
    """Sort by median_month desc, max_monthly_frequency desc, lineage asc."""
    return df.sort_values(
        ["median_month", "max_monthly_frequency", "lineage"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def main():
    haplotypes = pd.read_csv(snakemake.input.haplotypes, sep="\t")  # noqa: F821
    freqs = pd.read_csv(snakemake.input.freqs, sep="\t")  # noqa: F821
    print(f"Haplotypes: {len(haplotypes)} lineages")
    print(f"Frequencies: {len(freqs)} rows, {freqs['lineage'].nunique()} lineages")

    # --- Drop lineages with no frequency data ---
    freq_summaries = compute_freq_summaries(freqs)
    lineages_with_freqs = set(freq_summaries["lineage"])
    n_before = len(haplotypes)
    no_freq_lineages = haplotypes[~haplotypes["lineage"].isin(lineages_with_freqs)]
    haplotypes = haplotypes[haplotypes["lineage"].isin(lineages_with_freqs)].copy()
    n_dropped = n_before - len(haplotypes)
    print(f"\nDropped {n_dropped} lineages with no frequency counts in either source")
    if n_dropped > 0:
        print("  Concordance of dropped lineages:")
        for val, cnt in (
            no_freq_lineages["lanl_genbank_concordance"].value_counts().items()
        ):
            print(f"    {val}: {cnt}")

    haplotype_cols = [c for c in haplotypes.columns if c != "lineage"]

    # --- Uncollapsed ---
    uncollapsed = join_haplotypes_and_summaries(
        haplotypes, freq_summaries, haplotype_cols
    )
    uncollapsed = sort_output(uncollapsed)

    print(f"\nUncollapsed: {len(uncollapsed)} lineages")

    uncollapsed.to_csv(
        snakemake.output.uncollapsed,  # noqa: F821
        sep="\t",
        index=False,
        float_format="%.3g",
    )
    print(f"Wrote uncollapsed to {snakemake.output.uncollapsed}")  # noqa: F821

    # --- Collapsed ---
    # For each group, pick representative (highest max_monthly_frequency, break ties alphabetically)
    haplo_with_freq = haplotypes.merge(freq_summaries, on="lineage", how="left")
    haplo_with_freq["max_monthly_frequency"] = haplo_with_freq[
        "max_monthly_frequency"
    ].fillna(-1)

    n_unique_spikes = haplo_with_freq["complete_spike"].nunique()
    print(
        f"\nCollapsing: {len(haplotypes)} lineages -> {n_unique_spikes} unique spikes"
    )

    collapsed_rows = []
    for _, group in haplo_with_freq.groupby("complete_spike"):
        # Sort by max_monthly_frequency desc, then lineage asc to pick representative
        group_sorted = group.sort_values(
            ["max_monthly_frequency", "lineage"],
            ascending=[False, True],
        )
        representative = group_sorted.iloc[0]
        other_lineages = sorted(group_sorted.iloc[1:]["lineage"].tolist())

        # Re-aggregate frequencies for the group
        group_lineages = group["lineage"].tolist()
        group_freqs = freqs[freqs["lineage"].isin(group_lineages)]

        if len(group_freqs) > 0:
            # Sum lineage_counts per month, keep total_counts as-is (they're per-month totals)
            agg_freqs = group_freqs.groupby("month", as_index=False).agg(
                {
                    "lineage_counts_covspectrum": "sum",
                    "lineage_counts_genbank": "sum",
                    "total_counts_covspectrum": "first",
                    "total_counts_genbank": "first",
                }
            )
            agg_freqs["monthly_frequency_covspectrum"] = agg_freqs.apply(
                lambda r: (
                    r["lineage_counts_covspectrum"] / r["total_counts_covspectrum"]
                    if r["total_counts_covspectrum"] > 0
                    else 0.0
                ),
                axis=1,
            )
            agg_freqs["monthly_frequency_genbank"] = agg_freqs.apply(
                lambda r: (
                    r["lineage_counts_genbank"] / r["total_counts_genbank"]
                    if r["total_counts_genbank"] > 0
                    else 0.0
                ),
                axis=1,
            )
            agg_freqs["lineage"] = representative["lineage"]
            agg_freq_summary = compute_freq_summaries(agg_freqs)
        else:
            agg_freq_summary = pd.DataFrame()

        row_dict = {
            "lineage": representative["lineage"],
            "equivalent_lineages": ",".join(other_lineages) if other_lineages else "",
        }
        if len(agg_freq_summary) > 0:
            row_dict.update(agg_freq_summary.iloc[0].drop("lineage").to_dict())
        for col in haplotype_cols:
            row_dict[col] = representative[col]
        collapsed_rows.append(row_dict)

    collapsed_df = pd.DataFrame(collapsed_rows)

    # Arrange columns
    freq_cols = [
        "max_monthly_frequency",
        "max_monthly_frequency_covspectrum",
        "max_monthly_frequency_genbank",
        "median_month",
        "median_month_covspectrum",
        "median_month_genbank",
        "total_counts",
        "total_counts_covspectrum",
        "total_counts_genbank",
    ]
    collapsed_col_order = [
        "lineage",
        "equivalent_lineages",
        *freq_cols,
        *haplotype_cols,
    ]
    collapsed_df = collapsed_df.reindex(columns=collapsed_col_order)
    collapsed_df = sort_output(collapsed_df)

    n_with_equiv = (collapsed_df["equivalent_lineages"] != "").sum()
    print(f"Collapsed: {len(collapsed_df)} unique spikes")
    print(f"  With equivalent lineages: {n_with_equiv}")

    collapsed_df.to_csv(
        snakemake.output.collapsed,  # noqa: F821
        sep="\t",
        index=False,
        float_format="%.3g",
    )
    print(f"Wrote collapsed to {snakemake.output.collapsed}")  # noqa: F821

    # Write summary (same content as log)
    sys.stdout.flush()
    with open(snakemake.log[0]) as log_f:  # noqa: F821
        summary_text = log_f.read()
    with open(snakemake.output.summary, "w") as f:  # noqa: F821
        f.write(summary_text)


main()
