#!/usr/bin/env bash
set -euo pipefail

############################################
# Edit these paths before running
############################################
INPUT_DIR="/path/to/results/00_preprocess/host_removed"
OUTDIR="/path/to/results/01_gene_catalog"
NR_DB="/path/to/nr_tax.dmnd"
THREADS=8

############################################
# Method-matching parameters
############################################
MEGAHIT_MIN_CONTIG_LEN=300
ORF_MIN_NT_LEN=100
CDHIT_IDENTITY=0.90
CDHIT_COVERAGE=0.90
DIAMOND_EVALUE="1e-5"
SOAP_MAX_MISMATCHES=7   # ~95% identity for 150 bp PE reads; adjust if read length differs

mkdir -p \
  "${OUTDIR}/prodigal" \
  "${OUTDIR}/catalog" \
  "${OUTDIR}/annotation" \
  "${OUTDIR}/soap" \
  "${OUTDIR}/per_sample_gene_abundance" \
  "${OUTDIR}/logs"

for cmd in megahit prodigal cd-hit-est cd-hit diamond soap samtools python3; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "[ERROR] $cmd not found"; exit 1; }
done
command -v 2bwt-builder >/dev/null 2>&1 || { echo "[ERROR] 2bwt-builder (SOAP2 index builder) not found"; exit 1; }

shopt -s nullglob
R1S=("${INPUT_DIR}"/*.microbial.R1.fastq.gz)
R2S=("${INPUT_DIR}"/*.microbial.R2.fastq.gz)
[ ${#R1S[@]} -gt 0 ] || { echo "[ERROR] No microbial R1 files found in ${INPUT_DIR}"; exit 1; }
[ ${#R1S[@]} -eq ${#R2S[@]} ] || { echo "[ERROR] R1/R2 count mismatch"; exit 1; }

R1_CSV=$(IFS=,; echo "${R1S[*]}")
R2_CSV=$(IFS=,; echo "${R2S[*]}")

############################################
# 1. Co-assembly
############################################
megahit \
  -1 "$R1_CSV" \
  -2 "$R2_CSV" \
  --presets meta-sensitive \
  --min-contig-len "$MEGAHIT_MIN_CONTIG_LEN" \
  -t "$THREADS" \
  -o "${OUTDIR}/assembly" > "${OUTDIR}/logs/megahit.log" 2>&1

############################################
# 2. ORF prediction
############################################
prodigal \
  -i "${OUTDIR}/assembly/final.contigs.fa" \
  -d "${OUTDIR}/prodigal/all_orfs.fna" \
  -a "${OUTDIR}/prodigal/all_orfs.faa" \
  -o "${OUTDIR}/prodigal/all_orfs.gff" \
  -f gff \
  -p meta 2> "${OUTDIR}/logs/prodigal.log"

python3 "$(dirname "$0")/helper_filter_prodigal_orfs.py" \
  --nucleotide "${OUTDIR}/prodigal/all_orfs.fna" \
  --protein "${OUTDIR}/prodigal/all_orfs.faa" \
  --min-nt-len "$ORF_MIN_NT_LEN" \
  --out-nucleotide "${OUTDIR}/prodigal/orfs_min${ORF_MIN_NT_LEN}.fna" \
  --out-protein "${OUTDIR}/prodigal/orfs_min${ORF_MIN_NT_LEN}.faa"

############################################
# 3. Non-redundant catalog
############################################
cd-hit \
  -i "${OUTDIR}/prodigal/orfs_min${ORF_MIN_NT_LEN}.faa" \
  -o "${OUTDIR}/catalog/nr_catalog.faa" \
  -c "$CDHIT_IDENTITY" \
  -aS "$CDHIT_COVERAGE" \
  -T "$THREADS" -M 0 > "${OUTDIR}/logs/cdhit.log" 2>&1

python3 "$(dirname "$0")/helper_extract_fasta_by_ids.py" \
  --id-fasta "${OUTDIR}/catalog/nr_catalog.faa" \
  --source-fasta "${OUTDIR}/prodigal/orfs_min${ORF_MIN_NT_LEN}.fna" \
  --output-fasta "${OUTDIR}/catalog/nr_catalog.fna"

############################################
# 4. NR annotation
############################################
diamond blastp \
  --db "$NR_DB" \
  --query "${OUTDIR}/catalog/nr_catalog.faa" \
  --out "${OUTDIR}/annotation/diamond_nr.tsv" \
  --evalue "$DIAMOND_EVALUE" \
  --max-target-seqs 1 \
  --threads "$THREADS" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore staxids sscinames > "${OUTDIR}/logs/diamond.log" 2>&1

############################################
# 5. Build SOAP2 index for NR catalog
############################################
2bwt-builder "${OUTDIR}/catalog/nr_catalog.fna" > "${OUTDIR}/logs/soap_index.log" 2>&1

############################################
# 6. Per-sample mapping and gene abundance
############################################
python3 "$(dirname "$0")/helper_write_fasta_lengths.py" \
  --fasta "${OUTDIR}/catalog/nr_catalog.fna" \
  --output "${OUTDIR}/catalog/nr_catalog.lengths.tsv"

for R1 in "${INPUT_DIR}"/*.microbial.R1.fastq.gz; do
  sample=$(basename "$R1")
  sample="${sample%.microbial.R1.fastq.gz}"
  R2="${INPUT_DIR}/${sample}.microbial.R2.fastq.gz"
  [ -f "$R2" ] || { echo "[WARN] Missing R2 for ${sample}, skip"; continue; }

  soap_pe="${OUTDIR}/soap/${sample}.pe.soap"
  soap_se="${OUTDIR}/soap/${sample}.se.soap"
  log="${OUTDIR}/logs/${sample}.soap.log"

  soap \
    -a "$R1" \
    -b "$R2" \
    -D "${OUTDIR}/catalog/nr_catalog.fna.index" \
    -o "$soap_pe" \
    -2 "$soap_se" \
    -m 0 -x 1000 -s 40 -l 32 -v "$SOAP_MAX_MISMATCHES" -r 1 > "$log" 2>&1

  python3 "$(dirname "$0")/helper_parse_soap_to_gene_abundance.py" \
    --soap-pe "$soap_pe" \
    --soap-se "$soap_se" \
    --gene-lengths "${OUTDIR}/catalog/nr_catalog.lengths.tsv" \
    --sample "$sample" \
    --output "${OUTDIR}/per_sample_gene_abundance/${sample}.gene_abundance.tsv"
done

############################################
# 7. Merge per-sample gene abundance tables
############################################
python3 "$(dirname "$0")/02_build_gene_abundance_matrix.py" \
  --input-dir "${OUTDIR}/per_sample_gene_abundance" \
  --output "${OUTDIR}/gene_abundance_matrix.tsv"

echo "[INFO] Gene catalog, NR annotation, and gene abundance matrix finished."
