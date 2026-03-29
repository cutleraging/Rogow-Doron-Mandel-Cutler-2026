python filter_oRNAment.py \
/Users/ronaldcutler/EinsteinMed\ Dropbox/Ronald\ Cutler/Vijg-lab/Collaborations/Jackson\ Rogow/transcriptomics/annotations/oRNAment/Mus_musculus_cDNA_oRNAment.csv \
--no-header \
--gene-id 22602 \
--region "3;3" \
--min-score 0.75 \
--min-unpaired 0.5 \
-o oRNAment_arc.csv

python concat_rbp_pwms.py \
  --filtered oRNAment_arc.csv \
  --mapping /Users/ronaldcutler/EinsteinMed\ Dropbox/Ronald\ Cutler/Vijg-lab/Collaborations/Jackson\ Rogow/transcriptomics/annotations/oRNAment/RBP_id_encoding.csv \
  --pwm-dir /Users/ronaldcutler/EinsteinMed\ Dropbox/Ronald\ Cutler/Vijg-lab/Collaborations/Jackson\ Rogow/transcriptomics/annotations/oRNAment/PWMs \
  --out oRNAment_arc_rbp.meme

  grep -E '^MOTIF[[:space:]]' oRNAment_arc_rbp.meme | awk '{print $2}' > oRNAment_arc_rbp.txt