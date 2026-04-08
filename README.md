## Note

This repository contains the code, workflow, and documentation for the analyses presented in the paper:

**Bellwether hypothesis: dominant bacteria steer gut Microbial Nitrogen Metabolism via Additive Effects of Non-differential Species**

# metagenome_pipeline

A GitHub-ready shotgun metagenome workflow for taxonomic abundance analysis from paired-end metagenomic FASTQ files.

The computational workflow follows:

`fastp -> BWA host removal -> MEGAHIT co-assembly -> Prodigal ORF prediction -> ORF filtering -> CD-HIT non-redundant gene catalog -> SOAP2 mapping -> gene abundance matrix -> DIAMOND taxonomic annotation -> taxon abundance matrix -> LEfSe -> MaAsLin3`

This repository is organised in a stepwise structure so that each stage can be run, checked, and resumed independently.

---

## 1. Repository structure

```text
metagenome_nr_taxa_pipeline/
├── metadata/
│   └── groups.tsv.example
├── scripts/
│   ├── 00_preprocess_fastp_bwa.sh
│   ├── 01_gene_catalog_nr_annotation.sh
│   ├── 02_build_gene_abundance_matrix.py
│   ├── 03_sum_nr_taxa_abundance.py
│   ├── 04_lefse_taxa.R
│   ├── 05_maaslin3_taxa.R
│   ├── helper_filter_prodigal_orfs.py
│   ├── helper_extract_fasta_by_ids.py
│   ├── helper_write_fasta_lengths.py
│   └── helper_parse_soap_to_gene_abundance.py
└── README.md
```

### Main scripts

- `00_preprocess_fastp_bwa.sh`  
  Quality control with `fastp`, followed by host-read removal with `BWA`.

- `01_gene_catalog_nr_annotation.sh`  
  Co-assembly with `MEGAHIT`, ORF prediction with `Prodigal`, ORF filtering, non-redundant gene catalog construction with `CD-HIT`, per-sample mapping with `SOAP2`, gene abundance estimation, and DIAMOND taxonomic annotation.

- `02_build_gene_abundance_matrix.py`  
  Merges per-sample gene abundance tables into a single `gene_abundance_matrix.tsv`. This step is called automatically inside `01_gene_catalog_nr_annotation.sh`.

- `03_sum_nr_taxa_abundance.py`  
  Aggregates gene abundances to taxon-level abundance using the DIAMOND annotation table.

- `04_lefse_taxa.R`  
  Performs LEfSe on the taxon abundance matrix.

- `05_maaslin3_taxa.R`  
  Performs MaAsLin3 on the taxon abundance matrix.

### Helper scripts

The `helper_*.py` scripts are internal utilities used by the main shell pipeline and do not usually need to be run manually.

---

## 2. Software requirements

The pipeline was organised for Linux / WSL and assumes the following tools are available in `PATH`:

- `fastp`
- `bwa`
- `samtools`
- `megahit`
- `prodigal`
- `cd-hit`
- `cd-hit-est`
- `diamond`
- `soap`
- `2bwt-builder`
- `python3`
- `Rscript` (for LEfSe and MaAsLin3)

### Important note for SOAP2

`soap` and `2bwt-builder` may not come from the same package source on modern Linux systems.
A practical setup is:

- install `soap` from the system package manager if available;
- use a tested `2bwt-builder` binary that can build indices successfully on your system.

If `soap` works but `2bwt-builder` does not, `01_gene_catalog_nr_annotation.sh` will fail at the per-sample mapping stage.

---

## 3. Input requirements

### 3.1 Raw reads

Input files must be paired-end metagenomic FASTQ files named as:

```text
Sample_R1.fastq.gz
Sample_R2.fastq.gz
```

or equivalently `.fq.gz`.

If your raw files use another naming convention (for example `CFMT1R1.fq.gz` / `CFMT1R2.fq.gz`), rename or copy them into a compatible format before running the pipeline.

### 3.2 Metadata

`metadata/groups.tsv` must contain at least:

```text
SampleID    Group
```

where `SampleID` matches the sample prefix used in the abundance matrix.

---

## 4. Reference databases

This workflow requires two external reference resources:

### 4.1 Host reference genome for read removal

Used in `00_preprocess_fastp_bwa.sh`.

For samples, download reference genome FASTA and build a BWA index before running the pipeline.
A convenient choice is a single-file genome FASTA such as an Ensembl DNA FASTA file.

Example layout:

```text
ref/host/
├── Species....fa
├── Species....fa.amb
├── Species....fa.ann
├── Species....fa.bwt
├── Species....fa.pac
└── Species....fa.sa
```

Build the index with:

```bash
bwa index /path/to/host.fa
```

### 4.2 DIAMOND protein database for taxonomic annotation

Used in `01_gene_catalog_nr_annotation.sh`.

#### Recommended for the full taxonomic workflow

Use a DIAMOND database that includes taxonomy information, because `03_sum_nr_taxa_abundance.py` expects DIAMOND output fields including:

- `staxids`
- `sscinames`

Therefore, the DIAMOND database must be built with taxonomy files.

A typical strict setup is:

1. prepare a protein FASTA reference;
2. prepare taxonomy mapping files (`prot.accession2taxid`, `names.dmp`, `nodes.dmp`);
3. build the database with:

```bash
diamond makedb \
  --in reference_proteins.fa \
  --db reference_taxonomy.dmnd \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnodes nodes.dmp \
  --taxonnames names.dmp
```

## 5. Directory layout expected during analysis

A practical project layout is:

```text
project/
├── data/
│   └── raw_fastq/
├── metadata/
│   └── groups.tsv
├── ref/
│   └── host/
├── db/
├── results/
└── scripts/
```

The shell scripts create output subdirectories automatically under `results/`.

---

## 6. Pipeline steps and outputs

### Step 00. Quality control and host read removal

Run:

```bash
bash scripts/00_preprocess_fastp_bwa.sh
```

Main outputs:

```text
results/00_preprocess/
├── fastp/
│   ├── Sample.clean.R1.fastq.gz
│   ├── Sample.clean.R2.fastq.gz
│   ├── Sample.fastp.html
│   └── Sample.fastp.json
├── host_removed/
│   ├── Sample.microbial.R1.fastq.gz
│   └── Sample.microbial.R2.fastq.gz
└── logs/
```

Description:

- trims and filters paired-end reads with `fastp`;
- removes host-derived reads by mapping against the host genome with `BWA`;
- writes the retained microbial reads to `host_removed/`.

### Step 01. Assembly, gene catalog construction, mapping, and annotation

Run:

```bash
bash scripts/01_gene_catalog_nr_annotation.sh
```

Main outputs:

```text
results/01_gene_catalog/
├── assembly/
│   └── final.contigs.fa
├── prodigal/
│   ├── all_orfs.fna
│   └── all_orfs.faa
├── catalog/
│   ├── nr_catalog.fna
│   └── nr_catalog.faa
├── soap/
├── per_sample_gene_abundance/
│   └── Sample.gene_abundance.tsv
├── annotation/
│   └── diamond_nr.tsv
├── gene_abundance_matrix.tsv
└── logs/
```

Description:

- performs co-assembly of all host-filtered reads with `MEGAHIT`;
- predicts ORFs on contigs with `Prodigal`;
- filters ORFs and keeps coding sequences of the required minimum length;
- constructs a non-redundant gene catalog with `CD-HIT`;
- builds SOAP2 indices for the catalog and maps each sample back to the catalog;
- generates per-sample gene abundance tables and then merges them into `gene_abundance_matrix.tsv`;
- annotates the non-redundant protein catalog with DIAMOND.

### Step 02. Build the merged gene abundance matrix

This step is executed automatically inside Step 01. If needed, it can also be run manually:

```bash
python scripts/02_build_gene_abundance_matrix.py \
  --input-dir results/01_gene_catalog/per_sample_gene_abundance \
  --output results/01_gene_catalog/gene_abundance_matrix.tsv
```

### Step 03. Summarise gene abundances to taxon abundance

Run:

```bash
python scripts/03_sum_nr_taxa_abundance.py \
  --gene-matrix results/01_gene_catalog/gene_abundance_matrix.tsv \
  --diamond results/01_gene_catalog/annotation/diamond_nr.tsv \
  --output results/02_taxa_abundance/taxa_abundance.tsv
```

Main output:

```text
results/02_taxa_abundance/
└── taxa_abundance.tsv
```

Description:

- reads the merged gene abundance matrix;
- reads the DIAMOND top-hit annotation table;
- uses the scientific taxon name (`sscinames`) to sum gene abundances into a taxon-level abundance matrix.

**Important:** this step requires a taxonomy-enabled DIAMOND database. If `diamond_nr.tsv` was generated from a database without taxonomy fields, this step will fail.

### Step 04. LEfSe analysis

Run:

```bash
Rscript scripts/04_lefse_taxa.R \
  --input results/02_taxa_abundance/taxa_abundance.tsv \
  --metadata metadata/groups.tsv \
  --outdir results/03_lefse \
  --group-col Group
```

### Step 05. MaAsLin3 analysis

Run:

```bash
Rscript scripts/05_maaslin3_taxa.R \
  --input results/02_taxa_abundance/taxa_abundance.tsv \
  --metadata metadata/groups.tsv \
  --outdir results/04_maaslin3 \
  --fixed-effects Group \
  --reference Group,NO
```

---

## 7. Minimal run order

For a full taxonomic analysis:

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

## 8. Parameters to edit before running

Before running, edit the path variables at the top of:

- `scripts/00_preprocess_fastp_bwa.sh`
- `scripts/01_gene_catalog_nr_annotation.sh`

These usually include:

- raw FASTQ directory
- host genome BWA index prefix
- output directory
- DIAMOND database path
- thread number
