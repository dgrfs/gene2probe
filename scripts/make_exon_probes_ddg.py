#!/usr/bin/env python3
"""Standalone CLI wrapper for the exon probe design notebook."""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Design exon probes for a gene using the same settings as "
            "notebooks/002_make_exon_probes_XIST.ipynb."
        )
    )
    parser.add_argument(
        "gene_ID",
        help="Target gene symbol/name to design probes for.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gene_ID = args.gene_ID

    import pandas as pd
    import pybedtools

    from gene2probe import (  # noqa: F401
        check_for_required_nts,
        detect_offtargets,
        filter_by_GC_content,
        generate_kmers,
        get_region_of_interest,
        get_sequence_stats,
        read_gtf,
        remove_overlaps,
        run_blast,
        write_fasta,
    )

    repo_root = Path(__file__).resolve().parents[1]

    # Match the notebook defaults exactly.
    mode = "transcript"
    out_dir = repo_root / "sample_run" / f"probeDesign_{gene_ID}_{mode}"
    out_dir.mkdir(parents=True, exist_ok=True)

    gtf = repo_root / "hg38_resources" / "hg38.ncbiRefSeq.gtf"
    fasta = repo_root / "hg38_resources" / "hg38.fa"
    snp_db = repo_root / "hg38_resources" / "hg38_snp151Common.bed"
    repeats = repo_root / "hg38_resources" / "hg38_rmsk.bed"
    gaps = repo_root / "hg38_resources" / "hg38_rmsk.bed"
    blast_db = repo_root / "hg38_resources" / "001_blastdb" / "hg38_ncbiRefSeq_transcripts_db"

    print("current working directory:", os.getcwd())

    blast_exec_path = Path.home() / ".miniforge3" / "envs" / "gene2probe_env" / "bin"
    if not blast_exec_path.is_dir():
        warnings.warn(
            f"BLAST executable directory not found: {blast_exec_path}",
            RuntimeWarning,
        )

    print("blast executable path:", blast_exec_path)

    probe_length = 50
    split_nt = 25
    min_GC = 0.44
    max_GC = 0.72
    required_nts = {24: "T"}
    probe_offset = 1000
    n_desired_probes = 3
    min_mismatches = 5

    LHS_pref = "CCTTGGCACCCGAGAATTCCA"
    LHS_suff = ""
    RHS_pref = "/5Phos/"
    RHS_suff = "CCCATATAAGAAA"

    gene_anno = read_gtf(str(gtf))
    roi_bed = get_region_of_interest(
        gene_anno,
        gene_ID,
        gene_id_type="gene_name",
        feature=mode,
    )

    kmers = generate_kmers(roi_bed, k=probe_length)
    kmers.to_csv(out_dir / "kmers_all.csv")

    kmers = remove_overlaps(kmers, str(repeats))
    kmers = remove_overlaps(kmers, str(gaps))
    kmers = remove_overlaps(kmers, str(snp_db))

    kmers_bed = pybedtools.BedTool.from_dataframe(kmers)
    kmers_seq = kmers_bed.sequence(fi=str(fasta), s=True)
    kmers_seq_stats = get_sequence_stats(kmers_seq.seqfn, probe_length, split_nt)
    kmers = pd.merge(kmers, kmers_seq_stats, left_index=True, right_index=True)

    if required_nts is not None:
        kmers["has_required_nts"] = check_for_required_nts(kmers, required_nts)
        kmers = kmers[kmers["has_required_nts"] == True].reset_index(drop=True)  # noqa: E712

    kmers.to_csv(out_dir / "kmers_candidates_unfiltered.csv")
    kmers = filter_by_GC_content(kmers, min_GC, max_GC)
    kmers.to_csv(out_dir / "kmers_candidates_filtered.csv")

    write_fasta(
        kmers["name"],
        kmers["transcript_seq"],
        str(out_dir / "kmers_candidates_filtered_transcript_seqs.fa"),
    )

    if split_nt is not None:
        kmers["transcript_seq_LHS"] = [seq[0:split_nt] for seq in kmers["transcript_seq"]]
        kmers["transcript_seq_RHS"] = [seq[split_nt:probe_length] for seq in kmers["transcript_seq"]]

        write_fasta(
            kmers["name"],
            kmers["transcript_seq_LHS"],
            str(out_dir / "kmers_candidates_filtered_transcript_seqs_LHS.fa"),
        )
        write_fasta(
            kmers["name"],
            kmers["transcript_seq_RHS"],
            str(out_dir / "kmers_candidates_filtered_transcript_seqs_RHS.fa"),
        )

    blast_res = {}
    blast_res["full"] = run_blast(
        fasta=str(out_dir / "kmers_candidates_filtered_transcript_seqs.fa"),
        blastdb=str(blast_db),
        path2blastn=str(blast_exec_path / "blastn"),
        outfile=str(out_dir / "kmers_candidates_filtered_blast_output.txt"),
    )
    if split_nt is not None:
        blast_res["LHS"] = run_blast(
            fasta=str(out_dir / "kmers_candidates_filtered_transcript_seqs_LHS.fa"),
            blastdb=str(blast_db),
            path2blastn=str(blast_exec_path / "blastn"),
            outfile=str(out_dir / "kmers_candidates_filtered_blast_output_LHS.txt"),
        )
        blast_res["RHS"] = run_blast(
            fasta=str(out_dir / "kmers_candidates_filtered_transcript_seqs_RHS.fa"),
            blastdb=str(blast_db),
            path2blastn=str(blast_exec_path / "blastn"),
            outfile=str(out_dir / "kmers_candidates_filtered_blast_output_RHS.txt"),
        )

    offtargets = []
    for key in blast_res:
        offtargets += detect_offtargets(blast_res[key], gene_ID, min_mismatches=min_mismatches)
    offtargets = list(set(offtargets))
    kmers = kmers[kmers["name"].isin(offtargets) == False].reset_index(drop=True)  # noqa: E712

    kmers = kmers.sort_values("longest_homopolymer", ascending=True).reset_index(drop=True)

    df = kmers.copy()
    selected_probes_list = []
    while len(selected_probes_list) < n_desired_probes and not df.empty:
        selected_probe = df.iloc[[0]].copy()
        selected_probes_list.append(selected_probe)

        start = int(selected_probe["start"].iloc[0])
        end = int(selected_probe["end"].iloc[0])

        df = (
            df.loc[
                (df["end"] < start - probe_offset)
                | (df["start"] > end + probe_offset)
            ]
            .reset_index(drop=True)
            .copy()
        )

    if selected_probes_list:
        selected_probes_df = pd.concat(selected_probes_list, axis=0).reset_index(drop=True)
    else:
        selected_probes_df = kmers.head(0).copy()

    if split_nt is not None and not selected_probes_df.empty:
        selected_probes_df["probe_seq_LHS"] = [
            LHS_pref + seq[0:split_nt] + LHS_suff for seq in selected_probes_df["probe_seq"]
        ]
        selected_probes_df["probe_seq_RHS"] = [
            RHS_pref + seq[split_nt:probe_length] + RHS_suff for seq in selected_probes_df["probe_seq"]
        ]

    selected_probes_df["gene_ID"] = gene_ID
    selected_probes_df.to_csv(out_dir / "kmers_selected_probes.csv")

    print(selected_probes_df)
    print(f"Saved outputs to: {out_dir}")


if __name__ == "__main__":
    main()
