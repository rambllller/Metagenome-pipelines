#!/usr/bin/env bash
set -euo pipefail

############################################
# Edit these paths before running
############################################
RAW_DIR="/path/to/data/raw_fastq"
HOST_BWA_INDEX_PREFIX="/path/to/ref/host/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
OUTDIR="/path/to/results/00_preprocess"
THREADS=8

############################################
# fastp parameters matching the described QC
############################################
MIN_LEN=50
QUAL_PHRED=20

mkdir -p "${OUTDIR}/fastp" "${OUTDIR}/host_removed" "${OUTDIR}/logs"

command -v fastp >/dev/null 2>&1 || { echo "[ERROR] fastp not found"; exit 1; }
command -v bwa >/dev/null 2>&1 || { echo "[ERROR] bwa not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found"; exit 1; }

shopt -s nullglob
R1_FILES=("${RAW_DIR}"/*_R1.fastq.gz "${RAW_DIR}"/*_R1.fq.gz "${RAW_DIR}"/*_1.fastq.gz "${RAW_DIR}"/*_1.fq.gz)
[ ${#R1_FILES[@]} -gt 0 ] || { echo "[ERROR] No R1 files found in ${RAW_DIR}"; exit 1; }

for R1 in "${R1_FILES[@]}"; do
  base=$(basename "$R1")
  sample="${base%%_R1.fastq.gz}"
  sample="${sample%%_R1.fq.gz}"
  sample="${sample%%_1.fastq.gz}"
  sample="${sample%%_1.fq.gz}"

  if [[ -f "${RAW_DIR}/${sample}_R2.fastq.gz" ]]; then
    R2="${RAW_DIR}/${sample}_R2.fastq.gz"
  elif [[ -f "${RAW_DIR}/${sample}_R2.fq.gz" ]]; then
    R2="${RAW_DIR}/${sample}_R2.fq.gz"
  elif [[ -f "${RAW_DIR}/${sample}_2.fastq.gz" ]]; then
    R2="${RAW_DIR}/${sample}_2.fastq.gz"
  elif [[ -f "${RAW_DIR}/${sample}_2.fq.gz" ]]; then
    R2="${RAW_DIR}/${sample}_2.fq.gz"
  else
    echo "[WARN] Cannot find R2 for ${sample}, skip"
    continue
  fi

  clean_r1="${OUTDIR}/fastp/${sample}.clean.R1.fastq.gz"
  clean_r2="${OUTDIR}/fastp/${sample}.clean.R2.fastq.gz"
  json="${OUTDIR}/fastp/${sample}.fastp.json"
  html="${OUTDIR}/fastp/${sample}.fastp.html"
  log="${OUTDIR}/logs/${sample}.preprocess.log"

  echo "[INFO] ${sample}: fastp"
  fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$clean_r1" \
    -O "$clean_r2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred "$QUAL_PHRED" \
    --length_required "$MIN_LEN" \
    --thread "$THREADS" \
    --json "$json" \
    --html "$html" > "$log" 2>&1

  echo "[INFO] ${sample}: BWA host removal"
  bwa mem -t "$THREADS" "$HOST_BWA_INDEX_PREFIX" "$clean_r1" "$clean_r2" 2>> "$log" \
    | samtools view -@ "$THREADS" -b -f 12 -F 256 - 2>> "$log" \
    | samtools sort -@ "$THREADS" -n -o "${OUTDIR}/host_removed/${sample}.name_sorted.bam" - 2>> "$log"

  samtools fastq -@ "$THREADS" \
    -1 "${OUTDIR}/host_removed/${sample}.microbial.R1.fastq.gz" \
    -2 "${OUTDIR}/host_removed/${sample}.microbial.R2.fastq.gz" \
    -0 /dev/null -s /dev/null -n \
    "${OUTDIR}/host_removed/${sample}.name_sorted.bam" >> "$log" 2>&1

  rm -f "${OUTDIR}/host_removed/${sample}.name_sorted.bam"
  echo "[INFO] ${sample}: done"
done

echo "[INFO] All samples finished."
