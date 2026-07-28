# Notebook 002 Multi-Gene Update

This note records the changes made to `notebooks/002_make_exon_probes_XIST.ipynb` so the notebook can process multiple Gene IDs in one run.

## What Changed

- Replaced the single-gene workflow with a `gene_ids` list.
- Added a reusable `process_gene(gene_id)` function that runs the full probe-design pipeline for one gene.
- Kept the existing parameter block for probe design, filtering, BLAST, and adapter settings.
- Wrote outputs into gene-specific folders:
  - `../sample_run/probeDesign_{gene}_{mode}/`
- Added combined exports at the root output directory:
  - `probeDesign_{mode}_summary.csv`
  - `probeDesign_{mode}_selected_probes_all_genes.csv`

## Per-Gene Outputs

For each gene in `gene_ids`, the notebook now generates:

- `kmers_all.csv`
- `kmers_candidates_unfiltered.csv`
- `kmers_candidates_filtered.csv`
- `kmers_candidates_filtered_transcript_seqs.fa`
- `kmers_candidates_filtered_transcript_seqs_LHS.fa` and `kmers_candidates_filtered_transcript_seqs_RHS.fa` when probes are split
- `kmers_candidates_filtered_blast_output.txt`
- `kmers_candidates_filtered_blast_output_LHS.txt` and `kmers_candidates_filtered_blast_output_RHS.txt` when probes are split
- `kmers_selected_probes.csv`

## Behavior Notes

- The notebook still uses the existing `gene2probe` helper functions.
- The pipeline now loops over every gene in `gene_ids` and stores each result independently.
- The combined summary dataframe records basic counts for each gene, including how many k-mers survive each filtering stage.

