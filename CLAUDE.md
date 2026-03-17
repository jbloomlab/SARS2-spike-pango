# CLAUDE.md

## Principles for this repository
- Fail fast: unexpected outcomes should give clear error messages, not silent passing. Eg, get keys from dicts explicitly, do not use `.get` with defaults unless there is a specific implementation reason.
- Keep `CLAUDE.md` up-to-date about repo with details for Claude Code, and `README.md` up-to-date with information for user.
- Use Python and Snakemake
- The repository should use a single `environment.yaml` file that should have versions pinned at least to the first minor version (eg, `python=3.13`).
- All input data should go in `./data`, and all generated results should go in `./results`.
- Each `snakemake` rule should have a log that goes in `results/logs/<rule_name_w_wildcards>.txt`, for these logs for shell rules do simple redirection of `&> {log}` when possible. In scripts, do a `sys.stderr = sys.stdout = open(snakemake.log[0], "w")`.
- Every snakemake rule must have a `conda: "environment.yaml"` directive so `snakemake --lint` passes cleanly.
- Pass `config` values to shell commands via `params`, not directly as `config[...]` in the shell block.
- Lint with `snakemake --lint` and `ruff` and format with `snakefmt` and `black`.

## Pipeline overview

### Genbank data on lineage spike haplotypes and lineage frequencies
1. **download_genbank_seqs** / **extract_genbank_seqs**: Download and extract SARS-CoV-2 sequences from Genbank.
2. **parse_genbank_metadata**: Parse JSONL metadata into clean TSV, cross-check with genomic FASTA.
3. **download_usher_metadata**: Download UShER metadata with Pango lineage annotations.
4. **get_nextclade_dataset** / **run_nextclade**: Classify sequences by Nextstrain clade and Pango lineage, extract aligned spike proteins.
5. **parse_nextclade_results**: Parse nextclade TSV + spike FASTA into clean output. Drops rows with `errors` or `failedCdses` (no spike translation available). Validates `S:` prefix on mutation columns and joins aligned spike sequences. Script: `scripts/parse_nextclade_results.py`.
6. **merge_genbank_metadata_and_spikeseqs**: Merge GenBank metadata, parsed nextclade results, and UShER metadata into a single per-sequence TSV with spike sequences in `results/genbank_merged_data/`. Script: `scripts/merge_genbank_metadata_and_spikeseqs.py`.
7. **get_genbank_lineage_haplotypes**: Marimo notebook that filters sequences, assigns consensus Pango lineages, and computes the consensus spike haplotype per lineage using position-by-position majority vote from aligned spike sequences (excluding missing positions encoded as X). Insertions are resolved by majority vote among sequences where the insertion site is not missing. Key output columns: `median_mutations_from_consensus` (weighted median distance from consensus ignoring missing sites), `minor_mutations` (non-consensus mutations above `minor_mutations_cutoff` with their frequencies). Config params include `reference_spike` (FASTA for calling mutations) and `minor_mutations_cutoff`. Notebook: `notebooks/get_genbank_lineage_haplotypes.py`. Outputs to `results/genbank_lineages/`: `spike_haplotypes.tsv` (sorted by lineage), `freqs_by_month.tsv` (columns: lineage, month, lineage_counts, total_counts, monthly_frequency).

### LANL data on lineage spike haplotypes
8. **download_lanl_haplotypes**: Download LANL consensus spike haplotype archive (`SPIKE.pcfx.Wuhan.tar.gz`) from `cov.lanl.gov`. Config key: `lanl_spike_haplotypes_url`. Output: `results/lanl_downloads/`.
9. **parse_lanl_haplotypes**: Extract `SPIKE.short.all.Wuhan.txt` from the tarball, parse `(consensus)` lines for each lineage regardless of HD value, validate mutations against Wuhan-Hu-1 reference spike, and reconstruct `aligned_spike` (reference-length with substitutions and `-` for deletions) and `complete_spike` (deletions removed, insertions inserted). Lineages with invalid names (e.g., "None") are skipped. Script: `scripts/parse_lanl_haplotypes.py`. Output: `results/lanl_lineages/spike_haplotypes.tsv` (columns: lineage, total_substitutions, total_mutations, substitutions, deletions, insertions, complete_spike, aligned_spike), `results/lanl_lineages/get_lineage_haplotypes_summary.txt`.

### CovSpectrum data on lineage frequencies
10. **download_covspectrum_counts**: Download aggregated lineage counts by date from CovSpectrum LAPIS API. Config key: `covspectrum_aggregated_url`. Output: `results/covspectrum_downloads/aggregated_counts.tsv`.
11. **parse_covspectrum_freqs**: Aggregate daily counts to monthly frequencies. Drops rows with missing date or lineage, rounds dates to months (day >= 16 rounds up, matching Genbank convention). Script: `scripts/parse_covspectrum_freqs.py`. Output: `results/covspectrum_lineages/freqs_by_month.tsv` (columns: lineage, month, lineage_counts, total_counts, monthly_frequency), `results/covspectrum_lineages/parse_covspectrum_freqs_summary.txt`.

### Integrated outputs
12. **integrate_spike_haplotypes**: Outer join Genbank and LANL `spike_haplotypes.tsv` on lineage, classify concordance (`equal`, `differ`, `genbank_only`, `lanl_only`), prefer LANL values when sources differ. Also adds boolean annotations from external TSVs configured under `add_annotations` in `config.yaml` (each maps an annotation name to a TSV file with a `lineage` column; lineages in the TSV get "yes", others "no"; errors if the TSV contains lineages not in the merged data). Script: `scripts/integrate_spike_haplotypes.py`. Output: `results/lineages/spike_haplotypes.tsv` (columns: lineage, lanl_genbank_concordance, total_substitutions, total_mutations, substitutions, deletions, insertions, lanl_genbank_differences, ambiguous_sites_or_premature_stop, annotation columns from `add_annotations`, aligned_spike, complete_spike), `results/lineages/integrate_spike_haplotypes_summary.txt`.

13. **integrate_freqs_by_month**: Full outer join of Genbank and CovSpectrum `freqs_by_month.tsv` on `(lineage, month)`. Missing counts filled with 0, `total_counts` looked up from source month totals. Script: `scripts/integrate_freqs_by_month.py`. Output: `results/lineages/freqs_by_month.tsv` (columns: lineage, month, lineage_counts_covspectrum, lineage_counts_genbank, total_counts_covspectrum, total_counts_genbank, monthly_frequency_covspectrum, monthly_frequency_genbank), `results/lineages/integrate_freqs_by_month_summary.txt`.

14. **integrate_haplotypes_and_freqs**: Combine integrated spike haplotypes with per-lineage frequency summaries. Drops lineages with no frequency counts in either source (LANL-only lineages absent from both sequence databases). Produces uncollapsed (one row per lineage) and collapsed (lineages with identical `complete_spike` grouped, frequencies re-aggregated) outputs. Representative lineage chosen by highest `max_monthly_frequency`. Script: `scripts/integrate_haplotypes_and_freqs.py`. Output: `results/lineages/spike_haplotypes_w_freqs_uncollapsed.tsv`, `results/lineages/spike_haplotypes_w_freqs_collapsed.tsv`, `results/lineages/integrate_haplotypes_and_freqs_summary.txt`.

### Trees
15. **prepare_tree_inputs** + **nextstrain-prot-titers-tree module**: Filter lineages from collapsed haplotypes using `query_str` from `config.yaml` `trees:` section (a pandas query string, or `null` for no filtering), write alignment FASTA + metadata TSV, then build Nextstrain Auspice JSON trees via the `nextstrain-prot-titers-tree` submodule. Uses `median_month` as the tree `date` column. Each tree's `nextstrain_prot_titers_tree_config` is passed to the submodule (with alignment, metadata, results_subdir, and auspice_json added by the pipeline). Script: `scripts/prepare_tree_inputs.py`. Output: `auspice/SARS2-spike-pango_{tree_name}.json`.
