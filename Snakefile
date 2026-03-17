import calendar
import datetime
import os


def current_round_year_month():
    """Get current year-month, rounded to nearest month boundary."""
    dt = datetime.datetime.today()
    days_in_month = calendar.monthrange(dt.year, dt.month)[1]
    if dt.day > days_in_month / 2:
        if dt.month == 12:
            return f"{dt.year+1:04d}-01"
        return f"{dt.year:04d}-{dt.month+1:02d}"
    return f"{dt.year:04d}-{dt.month:02d}"


configfile: "config.yaml"


rule all:
    input:
        "results/genbank_merged_data/per_sequence_metadata_and_spikeseqs.tsv.gz",
        "results/genbank_merged_data/per_sequence_metadata_and_spikeseqs_summary.txt",
        "results/genbank_lineages/get_lineage_haplotypes.html",
        "results/genbank_lineages/freqs_by_month.tsv",
        "results/genbank_lineages/freqs_by_month.html",
        "results/genbank_lineages/spike_haplotypes.tsv",
        "results/genbank_lineages/get_lineage_haplotypes_summary.txt",
        "results/lanl_lineages/spike_haplotypes.tsv",
        "results/lanl_lineages/get_lineage_haplotypes_summary.txt",
        "results/covspectrum_lineages/freqs_by_month.tsv",
        "results/covspectrum_lineages/parse_covspectrum_freqs_summary.txt",
        "results/lineages/spike_haplotypes.tsv",
        "results/lineages/integrate_spike_haplotypes_summary.txt",
        "results/lineages/freqs_by_month.tsv",
        "results/lineages/integrate_freqs_by_month_summary.txt",
        "results/lineages/spike_haplotypes_w_freqs_uncollapsed.tsv",
        "results/lineages/spike_haplotypes_w_freqs_collapsed.tsv",
        "results/lineages/integrate_haplotypes_and_freqs_summary.txt",
        expand(
            config["trees_prefix"] + "_{tree_name}.json",
            tree_name=config["trees"],
        ),


rule download_usher_metadata:
    """Download UShER SARS-CoV-2 public metadata with Pango lineage annotations."""
    output:
        metadata="results/usher_metadata/public-latest.metadata.tsv.gz",
        timestamp="results/usher_metadata/timestamp.txt",
    params:
        url=config["usher_metadata"],
    conda:
        "environment.yaml"
    log:
        "results/logs/download_usher_metadata.txt",
    shell:
        """
        curl -fL -o {output.metadata} {params.url} &> {log}
        date -Iseconds > {output.timestamp}
        """


rule download_genbank_seqs:
    """Download all SARS-CoV-2 genome, CDS, and protein sequences from Genbank."""
    output:
        zipfile="results/genbank_seqs/download.zip",
        timestamp="results/genbank_seqs/timestamp.txt",
    params:
        taxon_id=config["taxon_id"],
    conda:
        "environment.yaml"
    log:
        "results/logs/download_genbank_seqs.txt",
    shell:
        """
        datasets download virus genome taxon {params.taxon_id} \
            --include genome,cds,protein \
            --filename {output.zipfile} \
            &> {log}
        date -Iseconds > {output.timestamp}
        """


rule extract_genbank_seqs:
    """Extract and gzip a file from the downloaded Genbank zip."""
    input:
        zipfile="results/genbank_seqs/download.zip",
    output:
        seqfile="results/genbank_seqs/{filename}.gz",
    wildcard_constraints:
        filename=r".+\.(fna|faa|jsonl)",
    params:
        member=lambda wc: f"ncbi_dataset/data/{wc.filename}",
    threads: 4
    conda:
        "environment.yaml"
    log:
        "results/logs/extract_genbank_seqs_{filename}.txt",
    shell:
        "unzip -p {input.zipfile} {params.member} 2> {log} | pigz -p {threads} > {output.seqfile}"


rule parse_genbank_metadata:
    """Parse GenBank metadata JSONL and cross-check with genomic FASTA accessions."""
    input:
        data_report="results/genbank_seqs/data_report.jsonl.gz",
        genomic_fasta="results/genbank_seqs/genomic.fna.gz",
    output:
        metadata="results/genbank_seqs/parsed_metadata.tsv.gz",
        summary="results/genbank_seqs/parsed_metadata_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/parse_genbank_metadata.txt",
    script:
        "scripts/parse_genbank_metadata.py"


rule get_nextclade_dataset:
    """Download the nextclade reference dataset for SARS-CoV-2 classification."""
    output:
        timestamp="results/nextclade_dataset/timestamp.txt",
    params:
        dataset_name=config["nextclade_dataset"],
        dataset_dir=lambda _, output: os.path.dirname(output.timestamp),
    conda:
        "environment.yaml"
    log:
        "results/logs/get_nextclade_dataset.txt",
    shell:
        """
        nextclade dataset get \
            --name {params.dataset_name} \
            --output-dir {params.dataset_dir} \
            &> {log}
        date -Iseconds > {output.timestamp}
        """


rule run_nextclade:
    """Classify sequences by Nextstrain clade and Pango lineage using nextclade."""
    input:
        sequences="results/genbank_seqs/genomic.fna.gz",
        dataset_timestamp="results/nextclade_dataset/timestamp.txt",
    output:
        tsv="results/nextclade/nextclade.tsv.gz",
        spike="results/nextclade/nextclade.cds_translation.S.fasta.gz",
    params:
        dataset_dir=lambda _, input: os.path.dirname(input.dataset_timestamp),
        columns=",".join(
            [
                "seqName",
                "clade",
                "dynamic",
                "coverage",
                "alignmentScore",
                "qc.overallScore",
                "qc.overallStatus",
                "totalMissing",
                "errors",
                "insertions",
                "aaInsertions",
                "failedCdses",
                "unknownAaRanges",
                "aaSubstitutions",
                "aaDeletions",
                "aaInsertions",
                "frameShifts",
            ]
        ),
        translations=lambda _: "results/nextclade/nextclade.cds_translation.{cds}.fasta.gz",
    threads: 24
    conda:
        "environment.yaml"
    log:
        "results/logs/run_nextclade.txt",
    shell:
        """
        nextclade run \
            --input-dataset {params.dataset_dir} \
            --output-tsv {output.tsv} \
            --output-columns-selection {params.columns} \
            --cds-selection S \
            --output-translations {params.translations} \
            --jobs {threads} \
            {input.sequences} \
            &> {log}
        """


rule merge_genbank_metadata_and_spikeseqs:
    """Merge genbank metadata, nextclade results, and UShER metadata."""
    input:
        genbank="results/genbank_seqs/parsed_metadata.tsv.gz",
        nextclade="results/nextclade/parsed_nextclade_results.tsv.gz",
        usher="results/usher_metadata/public-latest.metadata.tsv.gz",
    output:
        tsv="results/genbank_merged_data/per_sequence_metadata_and_spikeseqs.tsv.gz",
        summary="results/genbank_merged_data/per_sequence_metadata_and_spikeseqs_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/merge_genbank_metadata_and_spikeseqs.txt",
    script:
        "scripts/merge_genbank_metadata_and_spikeseqs.py"


rule parse_nextclade_results:
    """Parse nextclade results TSV and spike FASTA into clean output."""
    input:
        tsv="results/nextclade/nextclade.tsv.gz",
        spike_fasta="results/nextclade/nextclade.cds_translation.S.fasta.gz",
    output:
        tsv="results/nextclade/parsed_nextclade_results.tsv.gz",
        summary="results/nextclade/parsed_nextclade_results_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/parse_nextclade_results.txt",
    script:
        "scripts/parse_nextclade_results.py"


rule get_genbank_lineage_haplotypes:
    """Compute consensus spike haplotypes and frequencies per Pango lineage from Genbank data."""
    input:
        notebook="notebooks/get_genbank_lineage_haplotypes.py",
        data="results/genbank_merged_data/per_sequence_metadata_and_spikeseqs.tsv.gz",
        reference_spike=config["lineage_haplotypes"]["reference_spike"],
    output:
        html="results/genbank_lineages/get_lineage_haplotypes.html",
        freqs_tsv="results/genbank_lineages/freqs_by_month.tsv",
        freqs_chart="results/genbank_lineages/freqs_by_month.html",
        haplotypes_tsv="results/genbank_lineages/spike_haplotypes.tsv",
        summary="results/genbank_lineages/get_lineage_haplotypes_summary.txt",
    params:
        drop_flanking_gaps=config["lineage_haplotypes"]["drop_flanking_gaps"],
        max_nextclade_qc_score=config["lineage_haplotypes"]["max_nextclade_qc_score"],
        max_n_missing_sites=config["lineage_haplotypes"]["max_n_missing_sites"],
        min_date=config["lineage_haplotypes"]["min_date"],
        max_date=current_round_year_month(),
        n_recent_months=config["lineage_haplotypes"]["n_recent_months"],
        pango_definitions=",".join(config["lineage_haplotypes"]["pango_definitions"]),
        minor_mutations_cutoff=config["lineage_haplotypes"]["minor_mutations_cutoff"],
    conda:
        "environment.yaml"
    log:
        "results/logs/get_genbank_lineage_haplotypes.txt",
    shell:
        """
        marimo export html --no-include-code \
            -o {output.html} \
            {input.notebook} \
            -- \
            -drop_flanking_gaps {params.drop_flanking_gaps} \
            -max_nextclade_qc_score {params.max_nextclade_qc_score} \
            -max_n_missing_sites {params.max_n_missing_sites} \
            -min_date {params.min_date} \
            -max_date {params.max_date} \
            -n_recent_months {params.n_recent_months} \
            -pango_definitions {params.pango_definitions} \
            -reference_spike {input.reference_spike} \
            -minor_mutations_cutoff {params.minor_mutations_cutoff} \
            -freqs_by_month_tsv {output.freqs_tsv} \
            -freqs_by_month_chart {output.freqs_chart} \
            -spike_haplotypes_tsv {output.haplotypes_tsv} \
            -summary_txt {output.summary} \
            &> {log}
        """


rule download_lanl_haplotypes:
    """Download LANL consensus spike haplotypes for Pango lineages."""
    output:
        tarball="results/lanl_downloads/lanl_spike_haplotypes.tar.gz",
        timestamp="results/lanl_downloads/timestamp.txt",
    params:
        url=config["lanl_spike_haplotypes_url"],
    conda:
        "environment.yaml"
    log:
        "results/logs/download_lanl_haplotypes.txt",
    shell:
        """
        curl -fL -o {output.tarball} {params.url} &> {log}
        date -Iseconds > {output.timestamp}
        """


rule parse_lanl_haplotypes:
    """Parse LANL consensus spike haplotypes into standardized TSV."""
    input:
        tarball="results/lanl_downloads/lanl_spike_haplotypes.tar.gz",
        reference_spike=config["lineage_haplotypes"]["reference_spike"],
    output:
        tsv="results/lanl_lineages/spike_haplotypes.tsv",
        summary="results/lanl_lineages/get_lineage_haplotypes_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/parse_lanl_haplotypes.txt",
    script:
        "scripts/parse_lanl_haplotypes.py"


rule integrate_spike_haplotypes:
    """Merge Genbank and LANL spike haplotypes with concordance check."""
    input:
        config["add_annotations"].values(),
        genbank="results/genbank_lineages/spike_haplotypes.tsv",
        lanl="results/lanl_lineages/spike_haplotypes.tsv",
    output:
        tsv="results/lineages/spike_haplotypes.tsv",
        summary="results/lineages/integrate_spike_haplotypes_summary.txt",
    params:
        add_annotations=config["add_annotations"],
    conda:
        "environment.yaml"
    log:
        "results/logs/integrate_spike_haplotypes.txt",
    script:
        "scripts/integrate_spike_haplotypes.py"


rule integrate_freqs_by_month:
    """Merge Genbank and CovSpectrum monthly lineage frequencies."""
    input:
        genbank="results/genbank_lineages/freqs_by_month.tsv",
        covspectrum="results/covspectrum_lineages/freqs_by_month.tsv",
    output:
        tsv="results/lineages/freqs_by_month.tsv",
        summary="results/lineages/integrate_freqs_by_month_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/integrate_freqs_by_month.txt",
    script:
        "scripts/integrate_freqs_by_month.py"


rule integrate_haplotypes_and_freqs:
    """Combine spike haplotypes with frequency summaries; produce collapsed version."""
    input:
        haplotypes="results/lineages/spike_haplotypes.tsv",
        freqs="results/lineages/freqs_by_month.tsv",
    output:
        uncollapsed="results/lineages/spike_haplotypes_w_freqs_uncollapsed.tsv",
        collapsed="results/lineages/spike_haplotypes_w_freqs_collapsed.tsv",
        summary="results/lineages/integrate_haplotypes_and_freqs_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/integrate_haplotypes_and_freqs.txt",
    script:
        "scripts/integrate_haplotypes_and_freqs.py"


rule download_covspectrum_counts:
    """Download aggregated lineage counts by date from CovSpectrum."""
    output:
        tsv="results/covspectrum_downloads/aggregated_counts.tsv",
        timestamp="results/covspectrum_downloads/timestamp.txt",
    params:
        url=config["covspectrum_aggregated_url"],
    conda:
        "environment.yaml"
    log:
        "results/logs/download_covspectrum_counts.txt",
    shell:
        """
        curl -fL -o {output.tsv} '{params.url}' &> {log}
        date -Iseconds > {output.timestamp}
        """


rule parse_covspectrum_freqs:
    """Aggregate CovSpectrum daily counts into monthly lineage frequencies."""
    input:
        tsv="results/covspectrum_downloads/aggregated_counts.tsv",
    output:
        tsv="results/covspectrum_lineages/freqs_by_month.tsv",
        summary="results/covspectrum_lineages/parse_covspectrum_freqs_summary.txt",
    conda:
        "environment.yaml"
    log:
        "results/logs/parse_covspectrum_freqs.txt",
    script:
        "scripts/parse_covspectrum_freqs.py"


rule prepare_tree_inputs:
    """Filter lineages and prepare alignment, metadata, and site map for tree building."""
    input:
        haplotypes_tsv="results/lineages/spike_haplotypes_w_freqs_collapsed.tsv",
    output:
        alignment="results/trees/{tree_name}/alignment.fasta",
        metadata="results/trees/{tree_name}/metadata.tsv",
    params:
        query_str=lambda wc: config["trees"][wc.tree_name]["query_str"],
    conda:
        "environment.yaml"
    log:
        "results/logs/prepare_tree_inputs_{tree_name}.txt",
    script:
        "scripts/prepare_tree_inputs.py"


for _tree_name, _tree_config in config["trees"].items():
    _module_name = "tree_" + _tree_name.replace("-", "_")

    module:
        name: _module_name
        snakefile:
            "nextstrain-prot-titers-tree/Snakefile"
        config:
            _tree_config["nextstrain_prot_titers_tree_config"] | {
                "alignment": f"results/trees/{_tree_name}/alignment.fasta",
                "metadata": f"results/trees/{_tree_name}/metadata.tsv",
                "results_subdir": f"results/trees/{_tree_name}/tree_results",
                "auspice_json": f"{config['trees_prefix']}_{_tree_name}.json",
            }

    use rule * from _module_name as _module_name*
