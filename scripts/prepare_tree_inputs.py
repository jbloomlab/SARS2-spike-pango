"""Filter lineage haplotypes and prepare inputs for nextstrain-prot-titers-tree."""

import sys

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821

query_str = snakemake.params.query_str  # noqa: F821

print(f"Reading {snakemake.input.haplotypes_tsv}")  # noqa: F821
df = pd.read_csv(snakemake.input.haplotypes_tsv, sep="\t")  # noqa: F821
print(f"  Total lineages: {len(df)}")

if query_str is None or query_str == "all":
    print("  No filtering applied (query_str is null or 'all')")
else:
    n_before = len(df)
    df = df.query(query_str)
    print(f"  After query {query_str!r}: {len(df)} (dropped {n_before - len(df)})")

if len(df) == 0:
    raise ValueError(f"No lineages remain after filtering with query_str={query_str!r}")

print(f"\nWriting {len(df)} lineages to output files")

# Write alignment FASTA
with open(snakemake.output.alignment, "w") as f:  # noqa: F821
    for _, row in df.iterrows():
        f.write(f">{row['lineage']}\n{row['aligned_spike']}\n")
print(f"  Wrote alignment: {snakemake.output.alignment}")  # noqa: F821

# Write metadata TSV with strain/date columns expected by tree module
if "date" in df.columns:
    raise ValueError("'date' column already exists in haplotypes TSV")
if "strain" in df.columns:
    raise ValueError("'strain' column already exists in haplotypes TSV")
metadata = df.drop(columns=["aligned_spike", "complete_spike"]).copy()
metadata.insert(0, "strain", metadata.pop("lineage"))
metadata.insert(1, "date", metadata["median_month"])
metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)  # noqa: F821
print(f"  Wrote metadata: {snakemake.output.metadata}")  # noqa: F821
