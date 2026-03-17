"""Merge Genbank and LANL spike haplotypes with concordance check.

Performs an outer join on lineage, classifies concordance between sources,
and produces a merged haplotype TSV preferring LANL values when sources differ.
"""

import sys

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821


def compute_mutation_diff(row_genbank, row_lanl, mut_type):
    """Compute mutations present in one source but not the other.

    Returns (genbank_only, lanl_only) as sorted lists of mutation strings.
    """
    gb_set = set(
        row_genbank[mut_type].split(",")
        if pd.notna(row_genbank[mut_type]) and row_genbank[mut_type]
        else []
    )
    lanl_set = set(
        row_lanl[mut_type].split(",")
        if pd.notna(row_lanl[mut_type]) and row_lanl[mut_type]
        else []
    )
    return sorted(gb_set - lanl_set), sorted(lanl_set - gb_set)


def main():
    genbank = pd.read_csv(snakemake.input.genbank, sep="\t")  # noqa: F821
    lanl = pd.read_csv(snakemake.input.lanl, sep="\t")  # noqa: F821
    print(f"Genbank: {len(genbank)} lineages")
    print(f"LANL: {len(lanl)} lineages")

    # Verify no duplicate lineages in either file
    if genbank["lineage"].duplicated().any():
        dups = genbank.loc[genbank["lineage"].duplicated(), "lineage"].tolist()
        raise ValueError(f"Duplicate lineages in Genbank: {dups}")
    if lanl["lineage"].duplicated().any():
        dups = lanl.loc[lanl["lineage"].duplicated(), "lineage"].tolist()
        raise ValueError(f"Duplicate lineages in LANL: {dups}")

    # Outer join on lineage
    merged = genbank.merge(
        lanl,
        on="lineage",
        how="outer",
        suffixes=("_genbank", "_lanl"),
        indicator=True,
    )
    print(f"Merged: {len(merged)} lineages")

    n_both = (merged["_merge"] == "both").sum()
    n_genbank_only = (merged["_merge"] == "left_only").sum()
    n_lanl_only = (merged["_merge"] == "right_only").sum()
    print(f"  In both sources: {n_both}")
    print(f"  Genbank only: {n_genbank_only}")
    print(f"  LANL only: {n_lanl_only}")

    # Classify concordance
    # For rows in both, compare complete_spike
    concordance = []
    for _, row in merged.iterrows():
        if row["_merge"] == "left_only":
            concordance.append("genbank_only")
        elif row["_merge"] == "right_only":
            concordance.append("lanl_only")
        elif row["complete_spike_genbank"] == row["complete_spike_lanl"]:
            concordance.append("equal")
        else:
            concordance.append("differ")
    merged["lanl_genbank_concordance"] = concordance

    n_equal = concordance.count("equal")
    n_differ = concordance.count("differ")
    print(f"  Equal: {n_equal}")
    print(f"  Differ: {n_differ}")

    # Note: concordance is based on complete_spike identity. Mutation strings
    # (substitutions, deletions, insertions) may differ even when complete_spike
    # is identical, because different mutation descriptions can produce the same
    # sequence (e.g., Genbank may call R158G while LANL calls E156G for
    # adjacent deletion+substitution events). Similarly, "differ" rows may have
    # identical mutation strings but different complete_spike due to missing
    # sites (X) in the Genbank consensus.

    # Shared columns to resolve between sources
    shared_cols = [
        "total_substitutions",
        "total_mutations",
        "substitutions",
        "deletions",
        "insertions",
        "aligned_spike",
        "complete_spike",
    ]

    # Compute lanl_genbank_differences for "differ" rows
    differences = []
    for _, row in merged.iterrows():
        if row["lanl_genbank_concordance"] != "differ":
            differences.append("")
            continue
        parts = []
        for mut_type in ["substitutions", "deletions", "insertions"]:
            gb_only, lanl_only = compute_mutation_diff(
                {mut_type: row[f"{mut_type}_genbank"]},
                {mut_type: row[f"{mut_type}_lanl"]},
                mut_type,
            )
            if gb_only:
                parts.append(f"genbank_only_{mut_type}: {','.join(gb_only)}")
            if lanl_only:
                parts.append(f"lanl_only_{mut_type}: {','.join(lanl_only)}")
        if not parts:
            parts.append("sequences_differ_at_ambiguous_sites")
        differences.append("; ".join(parts))
    merged["lanl_genbank_differences"] = differences

    # Resolve shared columns: use Genbank for equal/genbank_only, LANL for differ/lanl_only
    for col in shared_cols:
        gb_col = f"{col}_genbank"
        lanl_col = f"{col}_lanl"
        resolved = []
        for _, row in merged.iterrows():
            if row["lanl_genbank_concordance"] in ("equal", "genbank_only"):
                resolved.append(row[gb_col])
            else:
                resolved.append(row[lanl_col])
        merged[col] = resolved

    # Check for ambiguous sites (X) or premature stop (*) in the chosen complete_spike
    def _has_ambiguous_or_premature_stop(spike):
        if "X" in spike:
            return "yes"
        # Premature stop: * appearing before the last position
        if "*" in spike[:-1]:
            return "yes"
        return "no"

    merged["ambiguous_sites_or_premature_stop"] = merged["complete_spike"].apply(
        _has_ambiguous_or_premature_stop
    )

    # Add annotations from external TSVs
    add_annotations = dict(snakemake.params.add_annotations)  # noqa: F821
    annotation_cols = []
    for annot_name, annot_path in add_annotations.items():
        annot_df = pd.read_csv(annot_path, sep="\t")
        annot_lineages = set(annot_df["lineage"])
        missing = annot_lineages - set(merged["lineage"])
        if missing:
            raise ValueError(
                f"Annotation '{annot_name}' from {annot_path} has "
                f"{len(missing)} lineages not in the lineages being produced: "
                f"{sorted(missing)}"
            )
        merged[annot_name] = merged["lineage"].map(
            lambda x, al=annot_lineages: "yes" if x in al else "no"
        )
        n_yes = (merged[annot_name] == "yes").sum()
        n_no = (merged[annot_name] == "no").sum()
        print(f"\nAnnotation '{annot_name}' from {annot_path}:")
        print(f"  yes: {n_yes}, no: {n_no}")
        annotation_cols.append(annot_name)

    # Select output columns in order
    output_cols = [
        "lineage",
        "lanl_genbank_concordance",
        "total_substitutions",
        "total_mutations",
        "substitutions",
        "deletions",
        "insertions",
        "lanl_genbank_differences",
        "ambiguous_sites_or_premature_stop",
        *annotation_cols,
        "aligned_spike",
        "complete_spike",
    ]
    result = merged[output_cols].sort_values("lineage").reset_index(drop=True)

    result.to_csv(snakemake.output.tsv, sep="\t", index=False)  # noqa: F821
    print(f"\nWrote {len(result)} lineages to {snakemake.output.tsv}")  # noqa: F821

    # Print concordance summary
    print("\nConcordance summary:")
    for val in ["equal", "differ", "genbank_only", "lanl_only"]:
        n = (result["lanl_genbank_concordance"] == val).sum()
        print(f"  {val}: {n}")

    # Write summary (same content as log)
    sys.stdout.flush()
    with open(snakemake.log[0]) as log_f:  # noqa: F821
        summary_text = log_f.read()
    with open(snakemake.output.summary, "w") as f:  # noqa: F821
        f.write(summary_text)


main()
