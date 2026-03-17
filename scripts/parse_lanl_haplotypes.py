"""Parse LANL consensus spike haplotypes into standardized TSV.

Extracts consensus (HD=0) spike haplotypes from the LANL common forms file,
validates mutations against the Wuhan-Hu-1 reference, reconstructs full spike
sequences, and writes a TSV matching the Genbank haplotype output schema.
"""

import re
import sys
import tarfile

import pandas as pd
from Bio import SeqIO


sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821


def parse_consensus_lines(tar_path):
    """Extract consensus lines from SPIKE.short.all.Wuhan.txt in the tarball.

    Returns list of (lineage, mutation_string) tuples. The mutation_string is
    the content inside [...] for HD=0 consensus lines.
    """
    with tarfile.open(tar_path, "r:gz") as tf:
        member = tf.getmember("SPIKE.short.all.Wuhan.txt")
        f = tf.extractfile(member)
        text = f.read().decode("utf-8")

    results = []
    for line in text.splitlines():
        if not line.endswith("(consensus)"):
            continue
        # Format: lineage  lineage_count  form_count  pct  HD [mutations] (consensus)
        match = re.match(
            r"^(\S+)\s+\d+\s+\d+\s+[\d.]+%\s+\d+\s+\[([^\]]*)\]\s+\(consensus\)$",
            line.strip(),
        )
        if match is None:
            raise ValueError(f"Failed to parse consensus line: {line!r}")
        lineage = match.group(1)
        mutation_str = match.group(2)
        results.append((lineage, mutation_str))

    if not results:
        raise ValueError("No consensus lines found in LANL file")

    return results


def parse_mutations(mutation_str):
    """Parse a LANL mutation string into substitutions, deletions, and insertions.

    Returns (substitutions, deletions, insertions) where each is a list of strings.
    Substitutions: e.g. ["T19I", "N501Y"]
    Deletions: e.g. ["H69-", "V70-"]
    Insertions: e.g. ["214:EPE"] (position:residues format matching Genbank pipeline)
    """
    substitutions = []
    deletions = []
    insertions = []

    if not mutation_str.strip():
        return substitutions, deletions, insertions

    # Validate that the mutation string contains only expected characters
    allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789,+-")
    unexpected = set(mutation_str) - allowed
    if unexpected:
        raise ValueError(
            f"Unexpected characters in mutation string: {unexpected!r} "
            f"in {mutation_str!r}"
        )

    for mut in mutation_str.split(","):
        mut = mut.strip()
        if mut.startswith("+"):
            # Insertion: +214EPE -> 214:EPE
            ins_match = re.match(r"^\+(\d+)([A-Z]+)$", mut)
            if ins_match is None:
                raise ValueError(f"Failed to parse insertion: {mut!r}")
            pos = ins_match.group(1)
            residues = ins_match.group(2)
            insertions.append(f"{pos}:{residues}")
        elif mut.endswith("-"):
            # Deletion: H69- -> H69-
            del_match = re.match(r"^([A-Z])(\d+)-$", mut)
            if del_match is None:
                raise ValueError(f"Failed to parse deletion: {mut!r}")
            deletions.append(mut)
        else:
            # Substitution: T19I
            sub_match = re.match(r"^([A-Z])(\d+)([A-Z])$", mut)
            if sub_match is None:
                raise ValueError(f"Failed to parse substitution: {mut!r}")
            substitutions.append(mut)

    return substitutions, deletions, insertions


def validate_and_apply_mutations(ref_seq, substitutions, deletions, insertions, lineage):
    """Validate mutations against reference and reconstruct spike sequences.

    Returns (aligned_spike, complete_spike).
    - aligned_spike: reference length, with substitutions and '-' for deletions
    - complete_spike: deletions removed, insertions inserted
    """
    seq = list(ref_seq)
    ref_len = len(ref_seq)

    # Apply substitutions
    for mut in substitutions:
        ref_aa = mut[0]
        pos = int(mut[1:-1])
        mut_aa = mut[-1]
        if pos < 1 or pos > ref_len:
            raise ValueError(
                f"Lineage {lineage}: substitution {mut} position {pos} "
                f"out of range (1-{ref_len})"
            )
        if ref_seq[pos - 1] != ref_aa:
            raise ValueError(
                f"Lineage {lineage}: substitution {mut} expects reference "
                f"'{ref_aa}' at position {pos}, but found '{ref_seq[pos - 1]}'"
            )
        seq[pos - 1] = mut_aa

    # Apply deletions
    for mut in deletions:
        ref_aa = mut[0]
        pos = int(mut[1:-1])
        if pos < 1 or pos > ref_len:
            raise ValueError(
                f"Lineage {lineage}: deletion {mut} position {pos} "
                f"out of range (1-{ref_len})"
            )
        if ref_seq[pos - 1] != ref_aa:
            raise ValueError(
                f"Lineage {lineage}: deletion {mut} expects reference "
                f"'{ref_aa}' at position {pos}, but found '{ref_seq[pos - 1]}'"
            )
        seq[pos - 1] = "-"

    aligned_spike = "".join(seq)

    # Build complete_spike: remove deletions, then insert insertions
    # First, build a list of (position, residues) for insertions, sorted by position
    # descending so we can insert from right to left without shifting indices
    complete_seq = aligned_spike.replace("-", "")

    # For insertions, we need to find the correct position in the complete sequence.
    # The insertion position refers to the reference coordinate. We need to map
    # reference positions to complete_seq positions (accounting for deleted positions).
    if insertions:
        # Build mapping from reference 1-based position to complete_seq 0-based index.
        # Position p in reference maps to the index in complete_seq of the character
        # that was at position p (if not deleted), or the index where it would be
        # (for insertion after that position).
        ref_to_complete = {}
        complete_idx = 0
        for ref_idx in range(ref_len):
            ref_to_complete[ref_idx + 1] = complete_idx
            if seq[ref_idx] != "-":
                complete_idx += 1

        # Parse and sort insertions by position descending (insert right-to-left)
        ins_list = []
        for ins in insertions:
            pos_str, residues = ins.split(":")
            pos = int(pos_str)
            if pos < 1 or pos > ref_len:
                raise ValueError(
                    f"Lineage {lineage}: insertion {ins} position {pos} "
                    f"out of range (1-{ref_len})"
                )
            ins_list.append((pos, residues))

        ins_list.sort(key=lambda x: x[0], reverse=True)

        complete_list = list(complete_seq)
        for pos, residues in ins_list:
            # Insert after the position in the complete sequence
            insert_idx = ref_to_complete[pos]
            if seq[pos - 1] != "-":
                insert_idx += 1
            complete_list[insert_idx:insert_idx] = list(residues)

        complete_spike = "".join(complete_list)
    else:
        complete_spike = complete_seq

    return aligned_spike, complete_spike


def main():
    ref_record = next(
        SeqIO.parse(snakemake.input.reference_spike, "fasta")  # noqa: F821
    )
    ref_seq = str(ref_record.seq)

    all_consensus_lines = parse_consensus_lines(snakemake.input.tarball)  # noqa: F821
    print(f"Parsed {len(all_consensus_lines)} consensus lines from LANL file")

    # Skip lineages with no name or name "None"
    invalid_names = {"None", ""}
    skipped_lineages = [
        (lin, muts) for lin, muts in all_consensus_lines if lin in invalid_names
    ]
    consensus_lines = [
        (lin, muts) for lin, muts in all_consensus_lines if lin not in invalid_names
    ]
    for lin, muts in skipped_lineages:
        print(
            f"Skipping lineage with invalid name {lin!r}: mutations=[{muts}]"
        )

    # Check for duplicate lineages
    lineage_names = [lin for lin, _ in consensus_lines]
    duplicates = [lin for lin in lineage_names if lineage_names.count(lin) > 1]
    if duplicates:
        raise ValueError(f"Duplicate lineages in LANL file: {set(duplicates)}")

    rows = []
    n_with_insertions = 0
    n_with_deletions = 0
    n_subs_only = 0
    n_no_mutations = 0

    for lineage, mutation_str in consensus_lines:
        substitutions, deletions, insertions = parse_mutations(mutation_str)

        aligned_spike, complete_spike = validate_and_apply_mutations(
            ref_seq, substitutions, deletions, insertions, lineage
        )

        # Check for ambiguous amino acids
        ambiguous = set(aligned_spike.replace("-", "")) - set("ACDEFGHIKLMNPQRSTVWY*")
        if ambiguous:
            raise ValueError(
                f"Lineage {lineage}: ambiguous amino acids in consensus: "
                f"{ambiguous}"
            )

        total_substitutions = len(substitutions)
        total_deletions = len(deletions)
        total_insertion_residues = sum(
            len(ins.split(":")[1]) for ins in insertions
        )
        total_mutations = total_substitutions + total_deletions + total_insertion_residues

        if insertions:
            n_with_insertions += 1
        if deletions:
            n_with_deletions += 1
        if total_mutations == 0:
            n_no_mutations += 1
        elif not deletions and not insertions:
            n_subs_only += 1

        rows.append(
            {
                "lineage": lineage,
                "total_substitutions": total_substitutions,
                "total_mutations": total_mutations,
                "substitutions": ",".join(substitutions) if substitutions else "",
                "deletions": ",".join(deletions) if deletions else "",
                "insertions": ",".join(insertions) if insertions else "",
                "complete_spike": complete_spike,
                "aligned_spike": aligned_spike,
            }
        )

    df = pd.DataFrame(rows).sort_values("lineage").reset_index(drop=True)

    # Validate aligned_spike length is always reference length
    bad_len = df[df["aligned_spike"].str.len() != len(ref_seq)]
    if len(bad_len) > 0:
        raise ValueError(
            f"aligned_spike length mismatch for lineages: "
            f"{bad_len['lineage'].tolist()}"
        )

    df.to_csv(snakemake.output.tsv, sep="\t", index=False)  # noqa: F821
    print(f"Wrote {len(df)} lineages to {snakemake.output.tsv}")  # noqa: F821

    # Write summary
    summary_lines = [
        "LANL spike haplotype parsing summary",
        "=====================================",
        f"Total consensus lines in LANL file: {len(all_consensus_lines)}",
        f"Skipped lineages (invalid name): {len(skipped_lineages)}",
        f"Total lineages parsed: {len(df)}",
        f"Lineages with no mutations (identical to reference): {n_no_mutations}",
        f"Lineages with substitutions only: {n_subs_only}",
        f"Lineages with deletions: {n_with_deletions}",
        f"Lineages with insertions: {n_with_insertions}",
    ]
    if skipped_lineages:
        summary_lines.append("\nSkipped lineages:")
        for lin, muts in skipped_lineages:
            summary_lines.append(f"  name={lin!r}: mutations=[{muts}]")

    summary_text = "\n".join(summary_lines) + "\n"
    with open(snakemake.output.summary, "w") as f:  # noqa: F821
        f.write(summary_text)
    print(summary_text)


main()
