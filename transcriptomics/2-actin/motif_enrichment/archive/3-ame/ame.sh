#!/usr/bin/env bash

foreground="/Users/ronaldcutler/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop-unique/batch-2-actin/motif_enrichment/1-extract-seq/log2fc_greaterthan_5_genes_3utr.fasta"

background="/Users/ronaldcutler/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop-unique/batch-2-actin/motif_enrichment/1-extract-seq/log2fc_lessthan_0_genes_3utr.fasta"

cisbp="/Users/ronaldcutler/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/annotations/motif_databases/CISBP-RNA/Mus_musculus.dna_encoded.meme"

ame \
  --control "$background" \
  --method fisher \
  --scoring totalhits \
  --evalue-report-threshold 0.1 \
  --oc . \
  "$foreground" "$cisbp"