#!/usr/bin/env bash

foreground="/Users/ronaldcutler/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop-unique/batch-2-arc-rep1-removed/motif_enrichment/1-extract-seq/log2fc_greaterthan_5_padj_lessthan_0p05.fasta"
background="/Users/ronaldcutler/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop-unique/batch-2-arc-rep1-removed/motif_enrichment/1-extract-seq/log2fc_lessthan_0.fasta"

xstreme \
  --dna2rna \
  --p "$foreground" \
  --n "$background" \
  --minw 6 \
  --maxw 15 \
  --align right \
  --evt 0.05 \
  --oc xstreme_out