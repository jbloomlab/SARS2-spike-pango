import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full")


@app.cell
def _():
    import colorsys

    import altair as alt

    import marimo as mo

    import numpy

    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()

    def log_text(text, markdown=False, summary_only=False):
        """Append text to both marimo output and summary log."""
        if not summary_only:
            mo.output.append(mo.md(text) if markdown else text)
        if markdown and text.startswith("#"):
            log_text.lines.append("")
        log_text.lines.append(text)

    log_text.lines = []

    return alt, colorsys, log_text, mo, numpy, pd


@app.cell
def _(mo):
    # read parameters from CLI args
    args = mo.cli_args()
    max_nextclade_qc_score = int(args["max_nextclade_qc_score"])
    max_n_missing_sites = int(args["max_n_missing_sites"])
    min_date = args["min_date"]
    max_date = args["max_date"]
    n_recent_months = int(args["n_recent_months"])
    pango_definitions = args["pango_definitions"].split(",")
    drop_flanking_gaps = args["drop_flanking_gaps"]
    if not isinstance(drop_flanking_gaps, bool):
        raise ValueError(f"{drop_flanking_gaps=} is not a bool")

    # display ceilings are 2x the filter max
    nextclade_qc_score_ceiling = 2 * max_nextclade_qc_score
    n_missing_sites_ceiling = 2 * max_n_missing_sites

    # consensus parameters
    reference_spike_path = args["reference_spike"]
    minor_mutations_cutoff = float(args["minor_mutations_cutoff"])

    # output files
    freqs_by_month_tsv = args["freqs_by_month_tsv"]
    freqs_by_month_chart = args["freqs_by_month_chart"]
    spike_haplotypes_tsv = args["spike_haplotypes_tsv"]
    summary_txt = args["summary_txt"]
    return (
        drop_flanking_gaps,
        spike_haplotypes_tsv,
        max_date,
        max_n_missing_sites,
        max_nextclade_qc_score,
        min_date,
        minor_mutations_cutoff,
        n_missing_sites_ceiling,
        n_recent_months,
        nextclade_qc_score_ceiling,
        pango_definitions,
        freqs_by_month_chart,
        freqs_by_month_tsv,
        reference_spike_path,
        summary_txt,
    )


@app.cell
def _(
    alt,
    drop_flanking_gaps,
    log_text,
    max_date,
    max_n_missing_sites,
    max_nextclade_qc_score,
    min_date,
    mo,
    n_missing_sites_ceiling,
    nextclade_qc_score_ceiling,
    numpy,
    pd,
    reference_spike_path,
):
    log_text("# Get spike haplotypes for Pango lineages", markdown=True)

    log_text("## Read and filter per-sequence data", markdown=True)

    # read, collapsing entries with identical data

    log_text("### Read per-sequence data", markdown=True)

    per_sequence_metadata_and_spikes_tsv = (
        "results/genbank_merged_data/per_sequence_metadata_and_spikeseqs.tsv.gz"
    )
    log_text(
        f"Reading per-sequence metadata and spikes from {per_sequence_metadata_and_spikes_tsv=}"
    )
    columns = {
        "accession": "string",
        "date_month_genbank": "string",
        "date_month_usher": "string",
        "geographic_region": "string",
        "pango_lineage_genbank": "string",
        "pango_lineage_nextclade": "string",
        "pango_lineage_usher": "string",
        "nextclade_qc_score": float,
        "n_missing_sites": "Int64",
        "substitutions": "string",
        "deletions": "string",
        "insertions": "string",
        "missing_sites": "string",
        "aligned_spike": "string",
    }

    per_sequence_df = pd.read_csv(
        per_sequence_metadata_and_spikes_tsv,
        sep="\t",
        usecols=list(columns),
        dtype=columns,
    )
    log_text(f"Read {len(per_sequence_df):,} accessions.")
    assert len(per_sequence_df) == per_sequence_df["accession"].nunique()
    log_text("Now collapsing by identical sequence/metadata.")
    per_sequence_df = per_sequence_df.groupby(
        [c for c in columns if c != "accession"],
        dropna=False,
        as_index=False,
    ).aggregate(
        n_accessions=pd.NamedAgg("accession", "count"),
        representative_accession=pd.NamedAgg("accession", "first"),
    )
    log_text(
        f"After collapsing, {len(per_sequence_df):,} rows representing "
        f"{per_sequence_df['n_accessions'].sum():,} accessions."
    )

    log_text("### Filtering for sequences with aligned spike", markdown=True)
    n_pre_aligned_spike = per_sequence_df["n_accessions"].sum()
    per_sequence_df = per_sequence_df.query("aligned_spike.notnull()")
    n_post_aligned_spike = per_sequence_df["n_accessions"].sum()
    log_text(
        "\nFiltering for just columns with non-null `aligned_spike`: "
        f"retained {n_post_aligned_spike:,} of {n_pre_aligned_spike:,} "
        f"({n_post_aligned_spike / n_pre_aligned_spike * 100:.1f}%) accessions"
    )

    for filter_name, ceiling, maxval in [
        ("nextclade_qc_score", nextclade_qc_score_ceiling, max_nextclade_qc_score),
        ("n_missing_sites", n_missing_sites_ceiling, max_n_missing_sites),
    ]:
        log_text(f"### Filtering by `{filter_name}`", markdown=True)
        log_text(f"\nAnalyzing `{filter_name}` using {ceiling=} and {maxval=}")
        assert per_sequence_df[filter_name].notnull().all()
        filter_df = (
            per_sequence_df.assign(
                val=lambda x: x[filter_name].astype(int).clip(upper=ceiling)
            )
            .groupby("val", as_index=False)
            .aggregate(n_accessions=pd.NamedAgg("n_accessions", "sum"))
            .assign(filtered=lambda x: x["val"] > maxval)
        )
        filter_dist = (
            alt.Chart(filter_df)
            .encode(
                alt.X("val", title=filter_name, scale=alt.Scale(nice=False)),
                alt.Y("n_accessions", scale=alt.Scale(type="symlog")),
                alt.Color("filtered"),
                tooltip=["n_accessions", alt.Tooltip("val", title=filter_name)],
            )
            .mark_bar()
            .properties(width=500, height=200, title=f"distribution of `{filter_name}`")
        )
        mo.output.append(filter_dist)
        n_pre_filter = per_sequence_df["n_accessions"].sum()
        per_sequence_df = per_sequence_df[per_sequence_df[filter_name] <= maxval]
        n_post_filter = per_sequence_df["n_accessions"].sum()
        log_text(
            f"\nFiltering for just accessions with `{filter_name}` <= {maxval}: "
            f"retained {n_post_filter:,} of {n_pre_filter:,} "
            f"({n_post_filter / n_pre_filter * 100:.2f}%) accessions"
        )

    log_text(
        "### Filtering sequences with flanking gaps in aligned spike", markdown=True
    )
    if drop_flanking_gaps:
        has_leading_gap = per_sequence_df["aligned_spike"].str.startswith("-")
        has_trailing_gap = per_sequence_df["aligned_spike"].str.endswith("-")
        n_leading = per_sequence_df.loc[has_leading_gap, "n_accessions"].sum()
        n_trailing = per_sequence_df.loc[has_trailing_gap, "n_accessions"].sum()
        n_both = per_sequence_df.loc[
            has_leading_gap & has_trailing_gap, "n_accessions"
        ].sum()
        n_either = per_sequence_df.loc[
            has_leading_gap | has_trailing_gap, "n_accessions"
        ].sum()
        n_pre_flanking = per_sequence_df["n_accessions"].sum()
        log_text(
            f"\n`drop_flanking_gaps` is enabled. Of {n_pre_flanking:,} accessions:\n\n"
            f"- {n_leading:,} have a leading gap\n"
            f"- {n_trailing:,} have a trailing gap\n"
            f"- {n_both:,} have both\n"
            f"- {n_either:,} have either",
            markdown=True,
        )
        per_sequence_df = per_sequence_df[~(has_leading_gap | has_trailing_gap)]
        n_post_flanking = per_sequence_df["n_accessions"].sum()
        log_text(
            f"\nAfter dropping flanking-gap sequences: retained {n_post_flanking:,} of "
            f"{n_pre_flanking:,} ({n_post_flanking / n_pre_flanking * 100:.2f}%) accessions"
        )
    else:
        log_text("\n`drop_flanking_gaps` is disabled, skipping this filter.")

    log_text("### Filtering sequences with bad collection dates", markdown=True)
    log_text(
        "\nAnalyzing UShER vs Genbank month-rounded dates "
        "(`date_month_usher` vs `date_month_genbank`).\n"
        "First tabulate number of accessions by how these dates compare:"
    )
    date_status_df = (
        per_sequence_df.assign(
            date_status=lambda x: numpy.where(
                x["date_month_usher"].isnull() & x["date_month_genbank"].isnull(),
                "both null",
                numpy.where(
                    x["date_month_usher"].isnull() | x["date_month_genbank"].isnull(),
                    "one null",
                    numpy.where(
                        x["date_month_usher"].fillna("")
                        == x["date_month_genbank"].fillna(""),
                        "both non-null and equal",
                        "both non-null and different",
                    ),
                ),
            ),
        )
        .groupby("date_status")
        .aggregate(n_accessions=pd.NamedAgg("n_accessions", "sum"))
    )
    mo.output.append(date_status_df)
    log_text(date_status_df.to_string(max_rows=20), summary_only=True)
    log_text(
        "Dropping accessions with no date or inconsistent date, otherwise keeping non-null dates as `date`"
    )
    assert (
        "date" not in per_sequence_df.columns
    ), "date column already exists in per_sequence_df"
    n_pre_date = per_sequence_df["n_accessions"].sum()
    per_sequence_df = per_sequence_df.assign(
        date=lambda x: x["date_month_genbank"].where(
            (
                x["date_month_usher"].isnull()
                | (
                    x["date_month_genbank"].notnull()
                    & (x["date_month_genbank"] == x["date_month_usher"])
                )
            ),
            x["date_month_usher"].where(
                x["date_month_genbank"].isnull()
                | (x["date_month_genbank"] == x["date_month_usher"]),
                pd.NA,
            ),
        ),
    ).query("date.notnull()")
    n_post_date = per_sequence_df["n_accessions"].sum()
    log_text(
        f"Retained {n_post_date:,} of {n_pre_date:,} ({n_post_date / n_pre_date * 100:.2f}%) accessions that have a date."
    )
    log_text("Dropping accessions where `date` is not YYYY-MM (eg, missing month).")
    per_sequence_df = per_sequence_df[
        per_sequence_df["date"].str.fullmatch(r"\d{4}\-\d{2}")
    ]
    n_post_month = per_sequence_df["n_accessions"].sum()
    log_text(
        f"Retained {n_post_month:,} of {n_post_date:,} ({n_post_month / n_post_date * 100:.2f}%) accessions that have a month."
    )
    assert per_sequence_df["date"].str.fullmatch(r"\d{4}\-\d{2}").all()
    log_text(f"Filtering sequences by {min_date=} and {max_date=}")
    dates = (
        per_sequence_df.groupby("date", as_index=False)
        .aggregate(n_accessions=pd.NamedAgg("n_accessions", "sum"))
        .assign(
            filtered=lambda x: (x["date"] < min_date) | (x["date"] > max_date),
            date=lambda x: pd.to_datetime(x["date"]),
        )
    )
    dates_dist = (
        alt.Chart(dates)
        .encode(
            alt.X("date", nice=False, axis=alt.Axis(format="%Y-%m", labelAngle=-90)),
            alt.Y("n_accessions", scale=alt.Scale(type="symlog")),
            alt.Color("filtered"),
            tooltip=[alt.Tooltip("date:T", format="%Y-%m"), "n_accessions"],
        )
        .mark_area(interpolate="step")
        .properties(
            width=800,
            height=200,
            title="distribution of month-rounded dates for all accessions",
        )
    )
    mo.output.append(dates_dist)
    n_pre_date_range = per_sequence_df["n_accessions"].sum()
    per_sequence_df = per_sequence_df.query(
        "(date >= @min_date) and (date <= @max_date)"
    )
    n_post_date_range = per_sequence_df["n_accessions"].sum()
    log_text(
        f"\nFiltering for accessions with dates betweeen {min_date} and {max_date}: "
        f"retained {n_post_date_range:,} of {n_pre_date_range:,} "
        f"({n_post_date_range / n_pre_date_range * 100:.2f}%) accessions"
    )
    # Read reference spike and validate aligned_spike data
    log_text("### Read reference spike and validate aligned_spike data", markdown=True)
    with open(reference_spike_path) as _f:
        reference_spike = "".join(
            line.strip() for line in _f if not line.startswith(">")
        )
    reference_spike_len = len(reference_spike)
    log_text(f"Read reference spike from {reference_spike_path} ({reference_spike_len} residues)")

    def parse_missing_sites(missing_str):
        """Parse missing_sites string into a set of integer positions."""
        if pd.isna(missing_str) or missing_str == "":
            return set()
        sites = set()
        for part in str(missing_str).split(","):
            if "-" in part:
                start, end = part.split("-")
                sites.update(range(int(start), int(end) + 1))
            else:
                sites.add(int(part))
        return sites

    import re

    def extract_site(mutation_str):
        """Extract integer site number from a mutation like 'I326V' or 'H69-'."""
        m = re.search(r"\d+", mutation_str)
        return int(m.group())

    log_text("Running validation assertions on all filtered sequences...")
    for _, _row in per_sequence_df.iterrows():
        _spike = _row["aligned_spike"]
        _missing = parse_missing_sites(_row["missing_sites"])

        # Assert aligned spike length matches reference
        assert len(_spike) == reference_spike_len, (
            f"aligned_spike length {len(_spike)} != reference {reference_spike_len} "
            f"for accession {_row['representative_accession']}"
        )

        # Assert missing_sites ↔ X in aligned_spike
        _x_positions = {_i + 1 for _i, _c in enumerate(_spike) if _c == "X"}
        assert _x_positions == _missing, (
            f"X positions {_x_positions - _missing} not in missing_sites, or "
            f"missing_sites {_missing - _x_positions} not X in aligned_spike "
            f"for accession {_row['representative_accession']}"
        )

        # Assert no mutation called at a missing site
        for _col in ["substitutions", "deletions"]:
            if pd.notna(_row[_col]) and _row[_col] != "":
                for _mut in _row[_col].split(","):
                    _site = extract_site(_mut)
                    assert _site not in _missing, (
                        f"Mutation {_mut} at missing site {_site} "
                        f"for accession {_row['representative_accession']}"
                    )
        # Note: insertions at site N are inserted AFTER position N, so the
        # insertion site CAN be in missing_sites (the position itself is X but
        # an insertion after it is still detectable). No assertion needed here.

    log_text(
        f"All {len(per_sequence_df):,} distinct sequences passed validation: "
        f"aligned_spike length, missing_sites ↔ X consistency, no mutations at missing sites."
    )

    return per_sequence_df, reference_spike, reference_spike_len, parse_missing_sites


@app.cell
def _(
    alt, log_text, mo, n_recent_months, numpy, pango_definitions, pd, per_sequence_df
):
    log_text("## Assign lineages", markdown=True)

    log_text("### Compare pango lineage definitions", markdown=True)
    log_text(
        f"We have the following lineage definitions: {pango_definitions}. "
        "Below we plot how many accessions have concordant versus discordant assignments from these definitions."
    )

    def consensus_status(row, definitions):
        """Return a string describing how pango lineage columns agree/disagree.

        Groups definitions into equivalence classes by value, with null definitions
        listed separately. Works for any list of definition names.
        """
        columns = [f"pango_lineage_{d}" for d in definitions]
        non_null = [
            (d, row[c]) for d, c in zip(definitions, columns) if pd.notna(row[c])
        ]
        null_defs = [d for d, c in zip(definitions, columns) if pd.isna(row[c])]

        if not non_null:
            return "all null"

        # group non-null definitions by their value, preserving input order
        groups = {}
        for d, val in non_null:
            groups.setdefault(val, []).append(d)
        # sort groups: more members first, then by first-appearance order in definitions
        sorted_groups = sorted(
            groups.values(),
            key=lambda g: (-len(g), definitions.index(g[0])),
        )

        if len(sorted_groups) == 1 and not null_defs:
            return "all agree"

        parts = [" = ".join(g) for g in sorted_groups]
        status = " != ".join(parts)
        if null_defs:
            status += ", " + " & ".join(null_defs) + " null"
        return status

    pango_definitions_df = per_sequence_df.assign(
        pango_consensus_status=lambda x: x.apply(
            consensus_status, axis=1, definitions=pango_definitions
        ),
    )

    pango_consensus_per_date = (
        pango_definitions_df.groupby(["date", "pango_consensus_status"], as_index=False)
        .aggregate(n_accessions=pd.NamedAgg("n_accessions", "sum"))
        .query("n_accessions > 0")
        .assign(date=lambda x: pd.to_datetime(x["date"]))
    )
    # compute fraction per category within each date, then scale by log10(total)
    # so stacked segments sum to log10(total) with heights proportional to linear fraction
    date_totals = (
        pango_consensus_per_date.groupby("date", as_index=False)
        .aggregate(total_accessions=pd.NamedAgg("n_accessions", "sum"))
        .assign(log10_total=lambda x: numpy.log10(x["total_accessions"]))
    )
    pango_consensus_summary = (
        pango_consensus_per_date.merge(date_totals, on="date")
        .assign(
            log10_n_accessions=lambda x: (
                x["log10_total"] * x["n_accessions"] / x["total_accessions"]
            ),
        )
        .drop(columns=["total_accessions", "log10_total"])
        .melt(
            id_vars=["date", "pango_consensus_status"],
            value_vars=["n_accessions", "log10_n_accessions"],
            value_name="value",
            var_name="variable",
        )
    )

    # sort order: most total accessions first (bottom of stack, first in legend)
    status_sort_order = (
        pango_consensus_per_date.groupby("pango_consensus_status", as_index=False)
        .aggregate(total=pd.NamedAgg("n_accessions", "sum"))
        .sort_values("total", ascending=False)["pango_consensus_status"]
        .tolist()
    )

    pango_consensus_dist = (
        alt.Chart(pango_consensus_summary)
        .encode(
            alt.X(
                "date",
                scale=alt.Scale(nice=False),
                axis=alt.Axis(format="%Y-%m", labelAngle=-90),
            ),
            alt.Y("value", title=None),
            alt.Color("pango_consensus_status", sort=status_sort_order),
            alt.Row("variable", title=None, header=alt.Header(labelFontWeight="bold")),
            tooltip=[
                alt.Tooltip("date:T", format="%Y-%m"),
                "pango_consensus_status",
                alt.Tooltip("value", format=".3g"),
            ],
        )
        .mark_area()
        .properties(
            width=650,
            height=200,
            title=alt.TitleParams(
                "Pango lineage definition agreement among assignment methods",
                subtitle="For log10 facet, bar height is log10 counts and colors are linear fraction with that status.",
            ),
        )
        .resolve_scale(y="independent")
        .interactive(bind_y=False)
    )
    mo.output.append(pango_consensus_dist)

    lineage_cols = [f"pango_lineage_{d}" for d in pango_definitions]
    recent_months = sorted(pango_definitions_df["date"].unique())[-n_recent_months:]
    for label, date_filter in [
        ("Top discordant lineage assignment combinations (all time)", None),
        (
            f"Top discordant lineage assignment combinations (last {n_recent_months} months with data: {recent_months[0]} to {recent_months[-1]})",
            recent_months,
        ),
    ]:
        log_text(f"### {label}", markdown=True)
        log_text(
            "For each discordant pango_consensus_status, showing the 4 most common "
            "combinations of lineage assignments across methods."
        )
        filtered_df = pango_definitions_df
        if date_filter is not None:
            filtered_df = filtered_df[filtered_df["date"].isin(date_filter)]
        discordant_df = (
            filtered_df.query("pango_consensus_status not in ['all agree', 'all null']")
            .groupby(
                ["pango_consensus_status"] + lineage_cols, dropna=False, as_index=False
            )
            .aggregate(n_accessions=pd.NamedAgg("n_accessions", "sum"))
            .sort_values("n_accessions", ascending=False)
        )
        top4_per_status = (
            discordant_df.groupby("pango_consensus_status", sort=False)
            .head(4)
            .reset_index(drop=True)
        )
        for _status, _group in top4_per_status.groupby(
            "pango_consensus_status", sort=False
        ):
            log_text(f"**{_status}**", markdown=True)
            _group_display = _group[lineage_cols + ["n_accessions"]].reset_index(
                drop=True
            )
            mo.output.append(_group_display)
            log_text(_group_display.to_string(max_rows=20), summary_only=True)
    log_text("### Consensus lineage assignment", markdown=True)
    log_text(
        "Assign a consensus `lineage` column. If all non-null methods agree, use that "
        "lineage. Otherwise, take the most common non-null assignment, breaking ties by "
        f"priority order: {', '.join(pango_definitions)}. Only null if all methods are null."
    )
    assert (
        "lineage" not in pango_definitions_df.columns
    ), "lineage column already exists in pango_definitions_df"

    def consensus_pango_lineage(row, definitions):
        """Get consensus pango lineage from multiple assignment methods."""
        cols = [f"pango_lineage_{d}" for d in definitions]
        values = [(i, row[c]) for i, c in enumerate(cols) if pd.notna(row[c])]
        if not values:
            return pd.NA
        # count occurrences; track best (lowest) priority index per lineage
        counts = {}
        for priority, val in values:
            if val not in counts:
                counts[val] = {"count": 0, "priority": priority}
            counts[val]["count"] += 1
        # pick highest count, then lowest priority index for ties
        best = max(counts, key=lambda v: (counts[v]["count"], -counts[v]["priority"]))
        return best

    pango_definitions_df = pango_definitions_df.assign(
        lineage=lambda x: x.apply(
            consensus_pango_lineage,
            axis=1,
            definitions=pango_definitions,
        ),
    )
    n_pre_pango = pango_definitions_df["n_accessions"].sum()
    n_null_pango = pango_definitions_df.query("lineage.isnull()")[
        "n_accessions"
    ].sum()
    n_non_null_pango = pango_definitions_df.query("lineage.notnull()")[
        "n_accessions"
    ].sum()
    log_text(
        f"\nOf {n_pre_pango:,} accessions: {n_non_null_pango:,} "
        f"({n_non_null_pango / n_pre_pango * 100:.2f}%) have a consensus lineage, "
        f"{n_null_pango:,} ({n_null_pango / n_pre_pango * 100:.2f}%) are null."
    )
    pango_definitions_df = pango_definitions_df.query("lineage.notnull()")
    n_post_pango = pango_definitions_df["n_accessions"].sum()
    log_text(
        f"Filtering for non-null consensus lineage: retained {n_post_pango:,} "
        f"of {n_pre_pango:,} ({n_post_pango / n_pre_pango * 100:.2f}%) accessions"
    )
    return (pango_definitions_df,)


@app.cell
def _(
    parse_missing_sites,
    alt,
    log_text,
    minor_mutations_cutoff,
    mo,
    numpy,
    pango_definitions_df,
    pd,
    reference_spike,
    reference_spike_len,
):
    log_text("## Consensus spike sequence for each lineage", markdown=True)

    log_text("### Compute position-wise consensus for each lineage", markdown=True)

    log_text(
        "For each lineage, computing a position-by-position majority vote "
        "consensus from aligned spike sequences, excluding missing positions (X). "
        "Insertions are resolved by majority vote among sequences where the "
        "insertion site is not missing."
    )

    def parse_mutations(s):
        """Parse a comma-separated mutation string into a frozenset."""
        if pd.isna(s) or s == "":
            return frozenset()
        return frozenset(s.split(","))

    def apply_insertions(aligned_spike, insertions_str):
        """Insert sequences into aligned_spike at positions from insertions string."""
        if pd.isna(insertions_str) or insertions_str == "":
            return aligned_spike
        insertions = []
        for entry in insertions_str.split(","):
            pos_str, seq = entry.split(":")
            insertions.append((int(pos_str), seq))
        insertions.sort(key=lambda x: x[0], reverse=True)
        spike = aligned_spike
        for pos, seq in insertions:
            spike = spike[:pos] + seq + spike[pos:]
        return spike

    def consensus_aligned_spike(spikes_arr, weights):
        """Position-wise weighted majority consensus, excluding X."""
        n_seqs, n_pos = spikes_arr.shape
        consensus = []
        for j in range(n_pos):
            col = spikes_arr[:, j]
            informative = col != "X"
            if not informative.any():
                consensus.append("X")
                continue
            chars = col[informative]
            w = weights[informative]
            unique_chars = numpy.unique(chars)
            best_char = max(unique_chars, key=lambda c: w[chars == c].sum())
            consensus.append(best_char)
        return "".join(consensus)

    def consensus_insertions(insertion_strs, missing_sites_sets, weights):
        """Majority-vote consensus for insertions, excluding seqs with missing site."""
        # Collect all insertion sites and their sequences
        site_to_seqs = {}  # {site: {seq_string: weighted_count}}
        site_no_ins = {}  # {site: weighted_count_without_insertion}
        for i, ins_str in enumerate(insertion_strs):
            ins_dict = {}
            if pd.notna(ins_str) and ins_str != "":
                for part in ins_str.split(","):
                    site_str, seq = part.split(":")
                    ins_dict[int(site_str)] = seq
            # Record this sequence's contribution to all known sites later
            for site_str_seq_pair in (ins_dict.items() if ins_dict else []):
                site, seq = site_str_seq_pair
                if site not in site_to_seqs:
                    site_to_seqs[site] = {}
                    site_no_ins[site] = 0.0
            # Store for later pass
            pass

        # Two-pass: first collect all insertion sites
        all_sites = set()
        parsed_insertions = []
        for i, ins_str in enumerate(insertion_strs):
            ins_dict = {}
            if pd.notna(ins_str) and ins_str != "":
                for part in ins_str.split(","):
                    site_str, seq = part.split(":")
                    ins_dict[int(site_str)] = seq
                    all_sites.add(int(site_str))
            parsed_insertions.append(ins_dict)

        if not all_sites:
            return ""

        # For each site, vote among informative sequences
        consensus_ins = {}
        for site in sorted(all_sites):
            counts = {}  # {insertion_seq_or_None: weighted_count}
            for i, ins_dict in enumerate(parsed_insertions):
                if site in missing_sites_sets[i]:
                    continue  # not informative
                key = ins_dict.get(site)  # None means no insertion
                counts[key] = counts.get(key, 0.0) + weights[i]
            if not counts:
                continue  # all sequences have this site missing
            best = max(counts, key=lambda k: counts[k])
            if best is not None:
                consensus_ins[site] = best

        if not consensus_ins:
            return ""
        return ",".join(
            f"{site}:{seq}" for site, seq in sorted(consensus_ins.items())
        )

    def derive_mutations_from_consensus(consensus_spike, ref_spike):
        """Compare consensus aligned_spike to reference, return mutation strings."""
        substitutions = []
        deletions = []
        missing = []
        for i, (ref_aa, cons_aa) in enumerate(zip(ref_spike, consensus_spike)):
            pos = i + 1  # 1-indexed
            if cons_aa == "X":
                missing.append(pos)
            elif cons_aa == "-" and ref_aa != "-":
                deletions.append(f"{ref_aa}{pos}-")
            elif cons_aa != ref_aa and ref_aa != "-":
                substitutions.append(f"{ref_aa}{pos}{cons_aa}")
        return (
            ",".join(substitutions) if substitutions else pd.NA,
            ",".join(deletions) if deletions else pd.NA,
            format_missing_sites(missing),
            len(missing),
        )

    def format_missing_sites(positions):
        """Format list of positions into compressed range string like '136-147,326,403'."""
        if not positions:
            return pd.NA
        positions = sorted(positions)
        ranges = []
        start = positions[0]
        end = positions[0]
        for p in positions[1:]:
            if p == end + 1:
                end = p
            else:
                ranges.append(f"{start}-{end}" if start != end else str(start))
                start = end = p
        ranges.append(f"{start}-{end}" if start != end else str(start))
        return ",".join(ranges)

    # Compute consensus for each lineage
    _n_lineages = pango_definitions_df["lineage"].nunique()
    log_text(f"Computing position-wise consensus for {_n_lineages:,} lineages...")

    lineage_consensus = {}  # {lineage: {aligned_spike, insertions, ...}}
    for _lineage, _group in pango_definitions_df.groupby("lineage"):
        _spikes = numpy.array([list(s) for s in _group["aligned_spike"].values])
        _weights = _group["n_accessions"].values.astype(float)
        _missing_sets = [
            parse_missing_sites(m) for m in _group["missing_sites"].values
        ]

        _cons_spike = consensus_aligned_spike(_spikes, _weights)
        _cons_ins = consensus_insertions(
            _group["insertions"].values, _missing_sets, _weights
        )
        _subs, _dels, _miss_str, _n_miss = derive_mutations_from_consensus(
            _cons_spike, reference_spike
        )
        lineage_consensus[_lineage] = {
            "aligned_spike": _cons_spike,
            "insertions": _cons_ins if _cons_ins else pd.NA,
            "substitutions": _subs,
            "deletions": _dels,
            "missing_sites": _miss_str,
            "n_missing_sites": _n_miss,
            "total_n_accessions": int(_weights.sum()),
        }

    log_text(f"Computed consensus for {len(lineage_consensus):,} lineages.")

    # For each lineage, pick representative accession, compute complete_spike, and mutations_separating
    log_text("### Pick representative accessions and compute mutations separating", markdown=True)

    final_consensus = {}
    for _lineage, _group in pango_definitions_df.groupby("lineage"):
        _cons = lineage_consensus[_lineage]
        _cons_spike = _cons["aligned_spike"]
        _cons_ins = _cons["insertions"] if pd.notna(_cons["insertions"]) else ""
        _weights = _group["n_accessions"].values.astype(float)
        _subs = _cons["substitutions"]
        _dels = _cons["deletions"]
        _miss_str = _cons["missing_sites"]
        _n_miss = _cons["n_missing_sites"]

        # Pick representative accession: closest to consensus, fewest missing sites
        _best_acc = None
        _best_dist = float("inf")
        _best_n_missing = float("inf")
        for _, _row in _group.iterrows():
            _dist = 0
            for _i, (_c, _s) in enumerate(zip(_cons_spike, _row["aligned_spike"])):
                if _s != "X" and _s != _c:
                    _dist += 1
            # Check insertion agreement
            _row_ins = {}
            if pd.notna(_row["insertions"]) and _row["insertions"] != "":
                for _part in _row["insertions"].split(","):
                    _site_str, _seq = _part.split(":")
                    _row_ins[int(_site_str)] = _seq
            _cons_ins_dict = {}
            if _cons_ins:
                for _part in _cons_ins.split(","):
                    _site_str, _seq = _part.split(":")
                    _cons_ins_dict[int(_site_str)] = _seq
            _row_missing = parse_missing_sites(_row["missing_sites"])
            for _site in set(_cons_ins_dict) | set(_row_ins):
                if _site in _row_missing:
                    continue
                if _cons_ins_dict.get(_site) != _row_ins.get(_site):
                    _dist += 1
            _n_missing = _row["n_missing_sites"]
            if (_dist, _n_missing) < (_best_dist, _best_n_missing):
                _best_dist = _dist
                _best_n_missing = _n_missing
                _best_acc = _row["representative_accession"]
                _best_row = _row

        _complete_spike = apply_insertions(
            _cons_spike, _cons_ins if _cons_ins else pd.NA
        ).replace("-", "")

        # Compute mutations separating consensus and representative accession
        _rep_spike = _best_row["aligned_spike"]
        _rep_ins = {}
        if pd.notna(_best_row["insertions"]) and _best_row["insertions"] != "":
            for _part in _best_row["insertions"].split(","):
                _site_str, _seq = _part.split(":")
                _rep_ins[int(_site_str)] = _seq
        _cons_ins_dict = {}
        if _cons_ins:
            for _part in _cons_ins.split(","):
                _site_str, _seq = _part.split(":")
                _cons_ins_dict[int(_site_str)] = _seq
        _rep_missing = parse_missing_sites(_best_row["missing_sites"])

        _sep_muts = []
        # Aligned spike position differences (skip X in representative)
        for _i in range(reference_spike_len):
            _c_aa = _cons_spike[_i]
            _r_aa = _rep_spike[_i]
            if _r_aa == "X" or _c_aa == _r_aa:
                continue
            _ref_aa = reference_spike[_i]
            _pos = _i + 1
            # Mutation in consensus but not representative
            if _c_aa != _ref_aa and _c_aa != "-":
                # consensus has a substitution at this position
                if _r_aa == _ref_aa or _r_aa == "-":
                    _sep_muts.append((-_pos, f"-{_ref_aa}{_pos}{_c_aa}"))
            if _c_aa == "-" and _ref_aa != "-":
                # consensus has a deletion
                if _r_aa != "-":
                    _sep_muts.append((-_pos, f"-{_ref_aa}{_pos}-"))
            # Mutation in representative but not consensus
            if _r_aa != _ref_aa and _r_aa != "-":
                if _c_aa == _ref_aa or _c_aa == "-":
                    _sep_muts.append((_pos, f"+{_ref_aa}{_pos}{_r_aa}"))
            if _r_aa == "-" and _ref_aa != "-":
                if _c_aa != "-":
                    _sep_muts.append((_pos, f"+{_ref_aa}{_pos}-"))

        # Insertion differences (skip sites where representative has missing)
        for _site in sorted(set(_cons_ins_dict) | set(_rep_ins)):
            if _site in _rep_missing:
                continue
            _c_val = _cons_ins_dict.get(_site)
            _r_val = _rep_ins.get(_site)
            if _c_val == _r_val:
                continue
            if _c_val is not None and _r_val is None:
                _sep_muts.append((-_site, f"-{_site}:{_c_val}"))
            elif _c_val is None and _r_val is not None:
                _sep_muts.append((_site, f"+{_site}:{_r_val}"))
            else:
                # both have insertion at this site but different sequences
                _sep_muts.append((-_site, f"-{_site}:{_c_val}"))
                _sep_muts.append((_site, f"+{_site}:{_r_val}"))

        # Sort by absolute position
        _sep_muts.sort(key=lambda x: (abs(x[0]), x[1]))
        _sep_muts_str = ", ".join(_m for _, _m in _sep_muts) if _sep_muts else ""

        # Validate: build representative complete_spike and check differences
        # match the listed mutations
        _rep_complete_spike = apply_insertions(
            _rep_spike, _best_row["insertions"] if pd.notna(_best_row["insertions"]) else pd.NA
        ).replace("-", "")

        # Count aligned_spike differences (non-X)
        _n_spike_diffs = sum(
            1 for _i in range(reference_spike_len)
            if _rep_spike[_i] != "X" and _rep_spike[_i] != _cons_spike[_i]
        )
        # Count insertion differences (non-missing)
        _n_ins_diffs = 0
        for _site in set(_cons_ins_dict) | set(_rep_ins):
            if _site in _rep_missing:
                continue
            if _cons_ins_dict.get(_site) != _rep_ins.get(_site):
                _n_ins_diffs += 1
        # Count listed mutations (each +/- entry is one difference, but a site
        # where both have different insertions generates two entries for one diff)
        _n_listed_spike_muts = sum(
            1 for _, _m in _sep_muts if not _m.lstrip("+-").split(":")[0].isdigit()
            or ":" not in _m.lstrip("+-")
        )
        _n_listed_ins_muts = _n_ins_diffs  # insertion diffs counted directly

        assert _n_spike_diffs + _n_ins_diffs == _best_dist, (
            f"Lineage {_lineage}: spike diffs {_n_spike_diffs} + ins diffs "
            f"{_n_ins_diffs} != best_dist {_best_dist}"
        )

        final_consensus[_lineage] = {
            "aligned_spike": _cons_spike,
            "insertions": _cons_ins if _cons_ins else pd.NA,
            "substitutions": _subs,
            "deletions": _dels,
            "missing_sites": _miss_str,
            "n_missing_sites": _n_miss,
            "total_n_accessions": int(_weights.sum()),
            "representative_accession": _best_acc,
            "complete_spike": _complete_spike,
            "mutations_separating_consensus_and_representative_accession": _sep_muts_str,
        }

    # Compute median_mutations_from_consensus and minor_mutations
    log_text(
        "### Compute distance from consensus and minor mutations", markdown=True,
    )
    log_text(
        "Computing weighted median mutation distance from consensus for each lineage "
        "(ignoring missing sites), and identifying minor mutations above "
        f"{minor_mutations_cutoff:.0%} frequency."
    )

    for _lineage, _group in pango_definitions_df.groupby("lineage"):
        _cons = final_consensus[_lineage]
        _cons_spike = _cons["aligned_spike"]
        _cons_ins_dict = {}
        if pd.notna(_cons["insertions"]) and _cons["insertions"] != "":
            for _part in _cons["insertions"].split(","):
                _site_str, _seq = _part.split(":")
                _cons_ins_dict[int(_site_str)] = _seq

        # Distance from consensus per sequence
        _distances = []
        _weights_list = []
        for _, _row in _group.iterrows():
            _dist = 0
            for _i, (_c, _s) in enumerate(zip(_cons_spike, _row["aligned_spike"])):
                if _s != "X" and _s != _c:
                    _dist += 1
            _row_ins = {}
            if pd.notna(_row["insertions"]) and _row["insertions"] != "":
                for _part in _row["insertions"].split(","):
                    _site_str, _seq = _part.split(":")
                    _row_ins[int(_site_str)] = _seq
            _row_missing = parse_missing_sites(_row["missing_sites"])
            for _site in set(_cons_ins_dict) | set(_row_ins):
                if _site in _row_missing:
                    continue
                if _cons_ins_dict.get(_site) != _row_ins.get(_site):
                    _dist += 1
            _distances.append(_dist)
            _weights_list.append(_row["n_accessions"])

        _distances = numpy.array(_distances)
        _weights_arr = numpy.array(_weights_list, dtype=float)
        _order = _distances.argsort()
        _sorted_distances = _distances[_order]
        _cumulative_weights = _weights_arr[_order].cumsum()
        _half_total = _cumulative_weights[-1] / 2
        _median_idx = numpy.searchsorted(_cumulative_weights, _half_total, side="right")
        _median_dist = _sorted_distances[min(_median_idx, len(_sorted_distances) - 1)]
        _cons["median_mutations_from_consensus"] = int(_median_dist)

        # Minor mutations: position-by-position
        _spikes = numpy.array([list(s) for s in _group["aligned_spike"].values])
        _weights_col = _group["n_accessions"].values.astype(float)
        _minor_muts = []

        # Aligned spike positions
        for _j in range(reference_spike_len):
            _cons_char = _cons_spike[_j]
            _col = _spikes[:, _j]
            _informative = _col != "X"
            if not _informative.any():
                continue
            _total_w = _weights_col[_informative].sum()
            # Find non-consensus characters
            for _char in numpy.unique(_col[_informative]):
                if _char == _cons_char:
                    continue
                _frac = _weights_col[_informative & (_col == _char)].sum() / _total_w
                if _frac >= minor_mutations_cutoff:
                    _pos = _j + 1
                    _ref_aa = reference_spike[_j]
                    if _char == "-":
                        _mut_str = f"{_ref_aa}{_pos}-"
                    elif _cons_char == "-":
                        # consensus is deletion but this char is not
                        _mut_str = f"{_ref_aa}{_pos}{_char}"
                    else:
                        _mut_str = f"{_ref_aa}{_pos}{_char}"
                    _minor_muts.append((_frac, _mut_str))

        # Insertion minor mutations
        _all_ins_sites = set(_cons_ins_dict.keys())
        _missing_sets = [
            parse_missing_sites(m) for m in _group["missing_sites"].values
        ]
        _parsed_ins = []
        for _ins_str in _group["insertions"].values:
            _d = {}
            if pd.notna(_ins_str) and _ins_str != "":
                for _part in _ins_str.split(","):
                    _site_str, _seq = _part.split(":")
                    _d[int(_site_str)] = _seq
                    _all_ins_sites.add(int(_site_str))
            _parsed_ins.append(_d)

        for _site in sorted(_all_ins_sites):
            _cons_val = _cons_ins_dict.get(_site)  # None if no consensus insertion
            # Tally non-consensus states at this insertion site
            _state_counts = {}
            _total_informative = 0.0
            for _i, _ins_dict in enumerate(_parsed_ins):
                if _site in _missing_sets[_i]:
                    continue
                _seq_val = _ins_dict.get(_site)
                _state_counts[_seq_val] = (
                    _state_counts.get(_seq_val, 0.0) + _weights_col[_i]
                )
                _total_informative += _weights_col[_i]
            if _total_informative == 0:
                continue
            for _state, _count in _state_counts.items():
                if _state == _cons_val:
                    continue
                _frac = _count / _total_informative
                if _frac >= minor_mutations_cutoff:
                    if _state is not None:
                        _mut_str = f"{_site}:{_state}"
                    else:
                        _mut_str = f"no_{_site}:ins"
                    _minor_muts.append((_frac, _mut_str))

        # Sort by fraction descending and format
        _minor_muts.sort(key=lambda x: -x[0])
        _cons["minor_mutations"] = (
            ", ".join(f"{_m} ({_f:.2f})" for _f, _m in _minor_muts) if _minor_muts else ""
        )

    # Build pango_haplotypes DataFrame
    pango_haplotypes = pd.DataFrame.from_dict(final_consensus, orient="index")
    pango_haplotypes.index.name = "lineage"

    # Compute total_substitutions and total_mutations
    pango_haplotypes = pango_haplotypes.assign(
        total_substitutions=lambda x: x["substitutions"].apply(
            lambda s: len(parse_mutations(s))
        ),
        total_mutations=lambda x: (
            x["substitutions"].apply(lambda s: len(parse_mutations(s)))
            + x["deletions"].apply(lambda s: len(parse_mutations(s)))
            + x["insertions"].apply(lambda s: len(parse_mutations(s)))
        ),
    )

    # Flag lineages with ambiguous amino acids or premature stops in complete_spike
    pango_haplotypes = pango_haplotypes.assign(
        ambiguous_sites_or_premature_stop=lambda x: (
            x["complete_spike"].str.contains("X", regex=False)
            | x["complete_spike"].str[:-1].str.contains(r"\*", regex=True)
        ),
    )

    log_text(
        f"Computed consensus haplotypes for {len(pango_haplotypes):,} lineages."
    )
    log_text(
        f"Median mutations from consensus: "
        f"mean={pango_haplotypes['median_mutations_from_consensus'].mean():.2f}, "
        f"median={pango_haplotypes['median_mutations_from_consensus'].median():.0f}, "
        f"max={pango_haplotypes['median_mutations_from_consensus'].max():.0f}"
    )

    # Visualization: scatter of median_mutations_from_consensus vs total_n_accessions
    haplotype_plot_df = pango_haplotypes.reset_index()[
        [
            "lineage",
            "total_n_accessions",
            "median_mutations_from_consensus",
            "minor_mutations",
            "total_mutations",
        ]
    ]

    haplotype_chart = (
        alt.Chart(haplotype_plot_df)
        .encode(
            alt.X(
                "total_n_accessions:Q",
                title="total accessions",
                scale=alt.Scale(type="log"),
            ),
            alt.Y(
                "median_mutations_from_consensus:Q",
                title="median mutations from consensus",
                scale=alt.Scale(nice=False),
            ),
            tooltip=[
                "lineage",
                "total_n_accessions",
                alt.Tooltip("median_mutations_from_consensus", format=".1f"),
                "minor_mutations",
            ],
        )
        .mark_circle(opacity=0.6)
        .properties(
            width=500,
            height=500,
            title="Median mutations from consensus vs total accessions per lineage",
        )
    )
    mo.output.append(haplotype_chart)

    n_high_median_mut = (
        pango_haplotypes["median_mutations_from_consensus"] > 1
    ).sum()
    log_text(
        f"There are {n_high_median_mut} of {len(pango_haplotypes)} lineages with more than one "
        "median mutation from consensus. "
        "These are shown below sorted by most accessions first:"
    )
    high_median_mut_df = pango_haplotypes.query(
        "median_mutations_from_consensus > 1"
    ).sort_values("total_n_accessions", ascending=False)[
        [
            "total_n_accessions",
            "median_mutations_from_consensus",
            "n_missing_sites",
            "minor_mutations",
        ]
    ]
    mo.output.append(high_median_mut_df)
    log_text(high_median_mut_df.to_string(max_rows=20), summary_only=True)

    has_ambiguous = pango_haplotypes["ambiguous_sites_or_premature_stop"]
    n_ambiguous = has_ambiguous.sum()
    log_text(
        f"{n_ambiguous} of {len(pango_haplotypes)} lineages have an 'X' (ambiguous) "
        f"or '*' (premature stop) in their complete spike sequence."
    )
    if n_ambiguous > 0:
        ambiguous_cols = [
            "total_n_accessions",
            "total_mutations",
            "n_missing_sites",
        ]
        ambiguous_df = pango_haplotypes.query(
            "ambiguous_sites_or_premature_stop"
        ).sort_values("total_n_accessions", ascending=False)[ambiguous_cols]
        mo.output.append(ambiguous_df)
        log_text(ambiguous_df.to_string(max_rows=20), summary_only=True)

    # Subsection: Total mutations per haplotype
    log_text("### Total mutations per haplotype", markdown=True)

    log_text(
        f"Total mutations per lineage relative to reference: "
        f"mean={pango_haplotypes['total_mutations'].mean():.1f}, "
        f"median={pango_haplotypes['total_mutations'].median():.0f}, "
        f"max={pango_haplotypes['total_mutations'].max():.0f}"
    )

    return (pango_haplotypes,)


@app.cell
def _(
    alt,
    colorsys,
    log_text,
    mo,
    pango_definitions_df,
    freqs_by_month_chart,
    freqs_by_month_tsv,
    pd,
):
    log_text("## Analyze lineage frequencies over time", markdown=True)

    log_text("Computing frequency of each lineage over time.")
    pango_freq_df = (
        pango_definitions_df.groupby(["lineage", "date"], as_index=False)
        .aggregate({"n_accessions": "sum"})
        .assign(
            monthly_accession_total=lambda x: x.groupby("date")[
                "n_accessions"
            ].transform("sum"),
            monthly_frequency=lambda x: x["n_accessions"]
            / x["monthly_accession_total"],
        )
        .sort_values(["date", "monthly_frequency"], ascending=[True, False])
    )
    log_text(f"Writing frequencies over time to {freqs_by_month_tsv=}.")
    pango_freq_df.rename(
        columns={
            "lineage": "lineage",
            "date": "month",
            "n_accessions": "lineage_counts",
            "monthly_accession_total": "total_counts",
        },
    ).to_csv(
        freqs_by_month_tsv, sep="\t", float_format="%.3g", index=False
    )

    # compute per-lineage summary stats
    lineage_stats = (
        pango_freq_df.groupby("lineage", as_index=False)
        .aggregate(
            first_date=pd.NamedAgg("date", "min"),
            last_date=pd.NamedAgg("date", "max"),
            max_monthly_frequency=pd.NamedAgg("monthly_frequency", "max"),
        )
        .sort_values(["first_date", "lineage"])
        .assign(sort_order=lambda x: range(len(x)))
    )
    lineage_order = lineage_stats["lineage"].tolist()

    # assign colors via golden-ratio hue stepping so temporally adjacent lineages differ
    _n_lineages = len(lineage_order)
    _golden_ratio = (1 + 5**0.5) / 2
    colors = []
    for _i in range(_n_lineages):
        _hue = (_i * _golden_ratio) % 1.0
        _r, _g, _b = colorsys.hls_to_rgb(_hue, 0.5, 0.8)
        colors.append(f"#{int(_r * 255):02x}{int(_g * 255):02x}{int(_b * 255):02x}")
    color_scale = alt.Scale(domain=lineage_order, range=colors)

    # prepare plotting dataframes with datetime dates
    freq_plot_df = pango_freq_df.merge(lineage_stats, on="lineage").assign(
        date=lambda x: pd.to_datetime(x["date"]),
        first_date=lambda x: pd.to_datetime(x["first_date"]),
        last_date=lambda x: pd.to_datetime(x["last_date"]),
    )
    monthly_totals = freq_plot_df.drop_duplicates(subset=["date"])[
        ["date", "monthly_accession_total"]
    ]
    freq_plot_df = freq_plot_df.drop(columns=["monthly_accession_total"])

    # interval brush on top chart to zoom the bottom chart
    freq_brush = alt.selection_interval(encodings=["x"])

    x_axis_kwargs = dict(format="%Y-%m", labelAngle=-90)

    # top panel: total monthly accessions with brush for date range selection
    freq_totals_chart = (
        alt.Chart(monthly_totals)
        .encode(
            alt.X(
                "date:T",
                scale=alt.Scale(nice=False),
                axis=alt.Axis(title="date (drag to select range)", **x_axis_kwargs),
            ),
            alt.Y(
                "monthly_accession_total:Q",
                scale=alt.Scale(type="log"),
                title="monthly accessions",
            ),
            tooltip=[alt.Tooltip("date:T", format="%Y-%m"), "monthly_accession_total"],
        )
        .mark_line(point=alt.OverlayMarkDef(color="black"), color="black")
        .add_params(freq_brush)
        .properties(width=900, height=135)
    )

    # bottom panel: stacked area of lineage frequencies, filtered by brush
    # use transform_impute to fill missing (lineage, date) combinations with 0 frequency
    # hover highlights the lineage with a black outline
    freq_hover = alt.selection_point(
        on="pointerover", fields=["lineage"], empty=False
    )
    freq_stacked_chart = (
        alt.Chart(freq_plot_df)
        .transform_filter(freq_brush)
        .transform_impute(
            impute="monthly_frequency",
            key="date",
            groupby=[
                "lineage",
                "sort_order",
                "first_date",
                "last_date",
                "max_monthly_frequency",
            ],
            value=0,
        )
        .encode(
            alt.X(
                "date:T",
                scale=alt.Scale(nice=False),
                axis=alt.Axis(title="date", **x_axis_kwargs),
            ),
            alt.Y(
                "monthly_frequency:Q",
                stack="zero",
                title="lineage frequency",
                scale=alt.Scale(domain=[0, 1]),
            ),
            alt.Color("lineage:N", scale=color_scale, legend=None),
            alt.Order("sort_order:Q"),
            stroke=alt.condition(
                freq_hover, alt.value("black"), alt.value("transparent")
            ),
            strokeWidth=alt.condition(freq_hover, alt.value(3), alt.value(0)),
            tooltip=[
                alt.Tooltip("date:T", format="%Y-%m"),
                "lineage",
                alt.Tooltip("monthly_frequency:Q", format=".3f", title="frequency"),
                "n_accessions",
                alt.Tooltip(
                    "max_monthly_frequency:Q",
                    format=".3f",
                    title="max monthly frequency",
                ),
                alt.Tooltip(
                    "first_date:T", format="%Y-%m", title="first month observed"
                ),
                alt.Tooltip("last_date:T", format="%Y-%m", title="last month observed"),
            ],
        )
        .mark_area()
        .add_params(freq_hover)
        .properties(width=900, height=200)
    )

    freq_combined_chart = (
        alt.vconcat(
            freq_totals_chart,
            freq_stacked_chart,
        )
        .configure_axis(
            grid=False, labelFontSize=11, titleFontWeight="normal", titleFontSize=14
        )
        .properties(
            title=alt.TitleParams(
                "Lineage frequencies over time",
                subtitle="Drag on top plot to select date range; mouse over bottom plot for details on frequencies and lineage names.",
                anchor="middle",
                fontSize=16,
                subtitleFontSize=14,
            ),
        )
    )

    log_text(
        f"Saving chart to {freqs_by_month_chart} (too large to display inline)"
    )
    freq_combined_chart.save(freqs_by_month_chart)
    return ()


@app.cell
def _(
    spike_haplotypes_tsv,
    log_text,
    mo,
    pango_haplotypes,
    pd,
    summary_txt,
):
    log_text("## Final lineage spike haplotypes", markdown=True)

    log_text(
        f"Final spike haplotypes for {len(pango_haplotypes)} lineages:"
    )

    final_haplotypes_df = (
        pango_haplotypes[
            [
                "total_n_accessions",
                "representative_accession",
                "mutations_separating_consensus_and_representative_accession",
                "total_substitutions",
                "total_mutations",
                "n_missing_sites",
                "substitutions",
                "deletions",
                "insertions",
                "missing_sites",
                "median_mutations_from_consensus",
                "minor_mutations",
                "ambiguous_sites_or_premature_stop",
                "complete_spike",
                "aligned_spike",
            ]
        ]
    ).sort_index()

    log_text(
        f"Showing first 10 of {len(final_haplotypes_df)} lineages "
        "(full data written to file):"
    )
    mo.output.append(final_haplotypes_df.head(10))

    log_text(f"Writing to {spike_haplotypes_tsv=}")
    final_haplotypes_df.to_csv(spike_haplotypes_tsv, sep="\t", float_format="%.3g")

    log_text(f"Writing summary to {summary_txt=}")
    with open(summary_txt, "w") as _f:
        _f.write("\n".join(log_text.lines) + "\n")
    return


if __name__ == "__main__":
    app.run()
