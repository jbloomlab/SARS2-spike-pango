import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Get library lineages
    This notebook gets a subset of key or high frequency lineages that are designated as library lineages and written to a file.
    """)
    return


@app.cell
def _():
    import marimo as mo

    import altair as alt

    import pandas as pd

    mo.output.append(mo.md("# Choose library strains"))

    # input / output files
    lineage_tsv = "../results/lineages/spike_haplotypes_w_freqs_collapsed.tsv"
    library_tsv = "library_lineages.tsv"

    # parameters for choosing library lineages
    lineages_include = [ # list of lineages to manually include
        "B.1.351",  # Beta
    ]
    lineages_exclude = []  # list of lineages to manually exclude

    mo.output.append(mo.md("## Read all lineages"))
    lineages = pd.read_csv(lineage_tsv, sep="\t").assign(
        equivalent_lineages=lambda x: x["equivalent_lineages"].astype(str)
    ) 
    mo.output.append(f"Read {len(lineages)=} from {lineage_tsv=}")

    lineages = lineages.query("ambiguous_sites_or_premature_stop == 'no'")
    mo.output.append(f"Retained {len(lineages)=} after excluding ones w ambiguous sites, premature stops, or missing sites.")

    mo.output.append(mo.md("## Get lineages meeting queries"))
    queries = [
        "(max_monthly_frequency >= 0.1)",
        "(max_monthly_frequency >= 0.05) and (median_month >= '2023-01')",
        "(max_monthly_frequency >= 0.025) and (median_month >= '2025-10')",
    ]

    lineages["in_library"] = False
    for query in queries:
        mo.output.append(mo.md(f"### query: {query}"))
        matches = lineages.eval(query)
        n_match = matches.sum()
        newly_added = matches & ~lineages["in_library"]
        n_new = newly_added.sum()
        new_lineages_df = (
            lineages.loc[newly_added, ["lineage", "equivalent_lineages", "median_month", "max_monthly_frequency"]]
            .sort_values("median_month", ascending=False)
            .set_index("lineage")
        )
        lineages["in_library"] = lineages["in_library"] | matches
        mo.output.append(f"{n_match} of {len(lineages)} match query, with {n_new} new inclusions")    
        mo.output.append(new_lineages_df)


    # manually include lineages
    for lineage in lineages_include:
        if lineage not in lineages["lineage"].values:
            raise ValueError(f"lineages_include entry {lineage!r} not found in lineages")
    if lineages_include:
        include_mask = lineages["lineage"].isin(lineages_include)
        newly_added = include_mask & ~lineages["in_library"]
        n_new = newly_added.sum()
        new_lineages_df = (
            lineages.loc[newly_added, ["lineage", "equivalent_lineages", "median_month", "max_monthly_frequency"]]
            .sort_values("median_month", ascending=False)
            .set_index("lineage")
        )
        lineages["in_library"] = lineages["in_library"] | include_mask
        mo.output.append(mo.md(
            f"## Manual includes\n\n{len(lineages_include)} specified, {n_new} newly added."
        ))
        mo.output.append(new_lineages_df)

    # manually exclude lineages
    for lineage in lineages_exclude:
        if lineage not in lineages["lineage"].values:
            raise ValueError(f"lineages_exclude entry {lineage!r} not found in lineages")
    if lineages_exclude:
        exclude_mask = lineages["lineage"].isin(lineages_exclude)
        n_removed = (exclude_mask & lineages["in_library"]).sum()
        lineages.loc[exclude_mask, "in_library"] = False
        mo.output.append(mo.md(
            f"## Manual excludes\n\n{len(lineages_exclude)} specified, {n_removed} removed from library."
        ))

    mo.output.append(mo.md("## Analyze library and write to file"))

    library = lineages.query("in_library")
    mo.output.append(f"Library has {len(library)} lineages")

    discordant = library.query("lanl_genbank_concordance != 'equal'")
    if len(discordant):
        mo.output.append("Following library lineages are not equal in LANL and Genbank")
        mo.output.append(
            discordant[
                [
                    "lineage",
                    "equivalent_lineages",
                    "median_month",
                    "max_monthly_frequency",
                    "lanl_genbank_concordance",
                    "lanl_genbank_differences",
                ]
            ]
        )

    library = library[
        [
            'lineage',
            'equivalent_lineages',
            'max_monthly_frequency',
            'median_month',
            'total_substitutions',
            'total_mutations',
            'substitutions',
            'deletions',
            'insertions',
            'aligned_spike',
            'complete_spike',
        ]
    ]

    library_chart_df = library.drop(columns=["aligned_spike", "complete_spike"])
    library_chart = (
        alt.Chart(library_chart_df)
        .encode(
            alt.X("median_month:T"),
            alt.Y("max_monthly_frequency", scale=alt.Scale(type="log")),
            alt.Fill("total_mutations"),
            tooltip=[c for c in library_chart_df.columns if c != "median_month"],
        )
        .mark_circle(size=80, stroke="black")
        .properties(width=400, height=400, title="library lineages")
        .configure_axis(grid=False)
    )
    mo.output.append(library_chart)

    mo.output.append(f"Writing library strains to {library_tsv=}")
    (
        library
        .sort_values(["median_month", "max_monthly_frequency"], ascending=False)
        .to_csv(library_tsv, sep="\t", index=False, float_format="%.3g")
    )
    return (mo,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
