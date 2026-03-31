# metagenome_nr_taxa_pipeline

A GitHub-ready shotgun metagenome pipeline matching the computational part of the described method:

`fastp -> BWA host removal -> MEGAHIT co-assembly -> Prodigal ORF prediction -> CD-HIT non-redundant gene catalog -> SOAP2 mapping -> DIAMOND NR annotation -> taxon abundance matrix -> LEfSe -> MaAsLin3`

## Run

```bash
bash scripts/00_preprocess_fastp_bwa.sh
bash scripts/01_gene_catalog_nr_annotation.sh
python scripts/03_sum_nr_taxa_abundance.py \
  --gene-matrix results/01_gene_catalog/gene_abundance_matrix.tsv \
  --diamond results/01_gene_catalog/annotation/diamond_nr.tsv \
  --output results/02_taxa_abundance/taxa_abundance.tsv

Rscript scripts/04_lefse_taxa.R \
  --input results/02_taxa_abundance/taxa_abundance.tsv \
  --metadata metadata/groups.tsv \
  --outdir results/03_lefse \
  --group-col Group

Rscript scripts/05_maaslin3_taxa.R \
  --input results/02_taxa_abundance/taxa_abundance.tsv \
  --metadata metadata/groups.tsv \
  --outdir results/04_maaslin3 \
  --fixed-effects Group \
  --reference Group,NO
```

## Notes

- Edit the path variables at the top of `00_preprocess_fastp_bwa.sh` and `01_gene_catalog_nr_annotation.sh` before running.
- Input files are expected as paired-end FASTQ files named `Sample_R1.fastq.gz` and `Sample_R2.fastq.gz`.
- `metadata/groups.tsv` must contain at least `SampleID` and the grouping column used in R.
