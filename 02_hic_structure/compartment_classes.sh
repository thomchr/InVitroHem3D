#!/bin/bash
# Assign consensus OCRs and genes to domain-change classes, and split each class into
# loop-engaged ("_looped") and remaining elements.
#      and one BED per TAD / sub-TAD remodeling class in compareTAD/class_beds/.
# Output is read by compartment_tables.R.
set -euo pipefail

DATA="${EHT_DATA_DIR:-data}"
OUT="$DATA/hiC/class_peaks"
CONSENSUS="$DATA/snATACseq/macs2_out/consensus_regions.bed"
LOOPS="$DATA/hiC/loops"          # HE_D9_plus7 / HPC_Pro_D9_plus8 .hicLoops.resMerge.FDR_adjusted_withP.anno.gene2OCR.txt
mkdir -p "$OUT"

# OCRs engaged in promoter loops, adherent and HPC
for s in HE_D9_plus7 HPC_Pro_D9_plus8 ; do
	cut -f1 "$LOOPS/$s.hicLoops.resMerge.FDR_adjusted_withP.anno.gene2OCR.txt" | sort -u | sed "s/:/\t/g" > "$OUT/$s.loop.peaks.bed"
done
cat "$OUT"/*.loop.peaks.bed | cut -f1-3 | sort -k1,1 -k2,2n -u > "$OUT/loop.peaks.bed"
tail -q -n +2 "$LOOPS"/*.hicLoops.resMerge.FDR_adjusted_withP.anno.gene2OCR.txt > "$OUT/gene2OCR.all.txt"

# one class per BED: compartments, then TAD / sub-TAD classes
classes=("$DATA/hiC/AtoB.bed" "$DATA/hiC/BtoA.bed")
if [ -d "$DATA/hiC/compareTAD/class_beds" ]; then classes+=("$DATA"/hiC/compareTAD/class_beds/*.bed); fi

for bed in "${classes[@]}" ; do
	c=$(basename "$bed" .bed)
	intersectBed -u -a "$CONSENSUS" -b "$bed" | cut -f1-3 > "$OUT/$c.all.peaks.bed"
	intersectBed -u -a "$OUT/$c.all.peaks.bed" -b "$OUT/loop.peaks.bed" > "$OUT/${c}_looped.peaks.bed"
	intersectBed -v -a "$OUT/$c.all.peaks.bed" -b "$OUT/loop.peaks.bed" > "$OUT/$c.peaks.bed"
	awk -F'\t' 'BEGIN{OFS="\t"} {split($1, a, ":"); print a[1], a[2], a[3], $0}' "$OUT/gene2OCR.all.txt" \
		| intersectBed -u -a stdin -b "$bed" | cut -f4- > "$OUT/$c.gene2OCR.txt"
done

# unchanged compartments: consensus peaks outside both switching classes
cat "$DATA/hiC/AtoB.bed" "$DATA/hiC/BtoA.bed" | sort -k1,1 -k2,2n > "$OUT/switched.bed"
intersectBed -v -a "$CONSENSUS" -b "$OUT/switched.bed" | cut -f1-3 > "$OUT/Unchange.peaks.bed"
