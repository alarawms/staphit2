# alarawms/staphit2

[![GitHub Actions CI Status](https://github.com/alarawms/staphit2/actions/workflows/nf-test.yml/badge.svg)](https://github.com/alarawms/staphit2/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/alarawms/staphit2/actions/workflows/linting.yml/badge.svg)](https://github.com/alarawms/staphit2/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**alarawms/staphit2** is an nf-core-compatible bioinformatics pipeline for comprehensive MRSA (*Staphylococcus aureus*) genomic surveillance. It takes raw Illumina paired-end reads and produces assembly-based typing, antimicrobial resistance profiling, structured virulence characterization, plasmid reconstruction, phylogenetic analysis, outbreak cluster detection, and automated reporting.

Built on a species-agnostic core architecture, the pipeline is designed to support additional bacterial pathogens in the future via species descriptor files.

```
Reads → TrimGalore → Rasusa (subsampling) → FastQC
                │
                ├→ SKESA / SPAdes → QUAST → CheckM2 → QC Gate
                │        │
                │        ├→ MLST + spaTyper + SCCmec + agr
                │        ├→ AMRFinderPlus (genes + point mutations)
                │        ├→ ABRicate (resfinder + vfdb + plasmidfinder)
                │        ├→ MOB-suite (plasmid reconstruction)
                │        └→ Prokka → Panaroo → IQ-TREE → SNP-dists
                │
                └→ KMA (read-based AMR)
                          │
                          └→ Aggregator → Summary → Outbreak Clustering
                                                  → Run Report
                                                  → Visualization
                                                  → MultiQC
```

### Analysis steps

| Step | Tool | Description |
|------|------|-------------|
| Read QC | TrimGalore + FastQC | Adapter trimming and quality assessment |
| Read subsampling | Rasusa | Normalize coverage to target depth (default: 100x) |
| Assembly | SKESA or SPAdes | *De novo* genome assembly |
| Assembly QC | QUAST + CheckM2 | Assembly metrics + completeness/contamination check |
| MLST | mlst | 7-locus multi-locus sequence typing |
| spa typing | spaTyper | *spa* gene repeat typing |
| SCCmec typing | SCCmec Typer | Staphylococcal cassette chromosome *mec* classification |
| agr typing | agr Typer | Accessory gene regulator group assignment |
| Species ID | Mash | Genomic distance-based species confirmation |
| AMR genes | AMRFinderPlus + ABRicate + KMA | Assembly- and read-based resistance detection |
| Point mutations | AMRFinderPlus POINT | gyrA, parC, rpoB, fusA mutations with predicted phenotype |
| Virulence profiling | staphit-virulence | PVL, TSST-1, IEC type, enterotoxins, operon completeness |
| Plasmid reconstruction | MOB-suite | Full plasmids with mobility and gene-to-replicon assignment |
| Pangenome | Panaroo | Core genome alignment |
| Phylogenetics | IQ-TREE or FastTree | Maximum-likelihood tree |
| Bayesian phylodynamics | BEAST2 (optional) | Time-calibrated Bayesian phylogeny with strict clock and coalescent prior |
| SNP distances | SNP-dists | Pairwise SNP distance matrix |
| Outbreak clustering | staphit-cluster | Tiered clustering (≤5/15/40 SNPs) with epi annotation |
| Reporting | staphit-report | Automated Markdown run report |
| Visualization | staphit-visualize | Publication-ready figures |
| Dashboard | staphit-dashboard | Interactive Plotly Dash surveillance dashboard (optional) |

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

### Input

Prepare a samplesheet CSV with paired-end FASTQ paths:

```csv
sample,fastq_1,fastq_2
MRSA_001,/data/MRSA_001_R1.fastq.gz,/data/MRSA_001_R2.fastq.gz
MRSA_002,/data/MRSA_002_R1.fastq.gz,/data/MRSA_002_R2.fastq.gz
```

### Running the pipeline

```bash
# Basic run
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results

# Full surveillance run with metadata and antibiogram
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --metadata sample_metadata.csv \
    --antibiogram antibiogram.csv

# Resume after adding samples
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    -resume
```

### Fetching public genomes directly via the pipeline (--cc / --st)

Instead of `--input`, pass `--cc` or `--st` and the pipeline fetches public isolates automatically: `pubmlst_fetch.py` queries PubMLST + NCBI SRA + ENA for matching accessions, then `fasterq-dump` downloads and assembles the samplesheet in-flight.

```bash
# All CC97 isolates (expands ST97, ST1153, ST1154, …)
nextflow run alarawms/staphit2 \
    -profile docker \
    --cc 97 \
    --outdir results/cc97

# Specific ST with geographic filter
nextflow run alarawms/staphit2 \
    -profile docker \
    --st 97 \
    --fetch_country "Saudi Arabia" \
    --outdir results/st97_sa

# Cap downloads for a pilot run
nextflow run alarawms/staphit2 \
    -profile docker \
    --cc 97 \
    --fetch_max 20 \
    --outdir results/cc97_pilot
```

> **Note:** Run from outside the staphit2 directory, or set a dedicated `--outdir`, to avoid Nextflow config conflicts.  
> SRA downloads require `fasterq-dump` (sra-tools container pulled automatically).

### Metadata tools

Convert lab data (Vitek PDFs, clinical spreadsheets) into pipeline-ready CSVs:

```bash
# Vitek 2 PDFs → antibiogram with full MIC data
python bin/staphit-metadata convert --from-vitek-pdf vitek_pdfs/ -o antibiogram.csv

# External clinical XLSX → enrich metadata
python bin/staphit-metadata convert --from-external-xlsx clinical.xlsx \
    --metadata sample_metadata.csv -o supplement.csv
```

### Fetching public genomes (pubmlst_fetch.py + nf-core/fetchngs)

`bin/pubmlst_fetch.py` retrieves public *S. aureus* WGS records by ST, clonal complex, country, continent, host, or year from three sources in parallel: PubMLST BIGSdb, NCBI Entrez SRA, and ENA Portal.

**Step 1 — fetch metadata and accession list**

```bash
# All CC97 isolates (expands to member STs: 97, 1153, 1154, …)
python bin/pubmlst_fetch.py \
    --cc 97 --source all --sra-only \
    --out cc97_metadata.tsv \
    --accessions cc97_accessions.txt

# Specific ST only
python bin/pubmlst_fetch.py \
    --st 97 --source all --sra-only \
    --out st97_metadata.tsv \
    --accessions st97_accessions.txt

# Filter by geography or host
python bin/pubmlst_fetch.py \
    --cc 97 --continent europe --host "Homo sapiens" --sra-only \
    --out cc97_eu_human.tsv \
    --accessions cc97_eu_human_accessions.txt
```

Key flags:

| Flag | Description |
|------|-------------|
| `--cc CC` | Clonal complex — expands to all member STs via PubMLST |
| `--st ST[,ST2,…]` | One or more exact sequence types |
| `--source all\|pubmlst\|ncbi\|ena` | Data sources to query (default: all) |
| `--sra-only` | Discard isolates without a public run accession |
| `--accessions FILE` | One SRA/ENA run ID per line — input for nf-core/fetchngs |
| `--out FILE` | Full metadata TSV with provenance fields |

**Step 2 — download FASTQs with nf-core/fetchngs**

[nf-core/fetchngs](https://nf-co.re/fetchngs) handles ENA and NCBI downloads reliably (tries Aspera, then FTP, then HTTPS):

```bash
nextflow run nf-core/fetchngs \
    --input cc97_accessions.txt \
    --outdir cc97_fetchngs \
    --nf_core_pipeline rnaseq \
    -profile singularity
```

This produces `cc97_fetchngs/samplesheet/samplesheet.csv` in the `sample,fastq_1,fastq_2` format expected by staphit2.

**Step 3 — run staphit2**

```bash
nextflow run alarawms/staphit2 \
    -profile singularity \
    --input cc97_fetchngs/samplesheet/samplesheet.csv \
    --outdir results/cc97 \
    -resume
```

### Interactive dashboard

```bash
pip install dash plotly pandas dash-bootstrap-components
python bin/staphit-dashboard results/ --host 0.0.0.0 --port 8050
```

Two modes: **Investigation** (linked tree + map + table) and **Analysis** (resistance, trends, virulence, clusters, plasmids, QC tabs).

### Annotated phylogenetic tree (local)

`bin/plot_tree.R` reads pipeline outputs directly and produces a PDF and SVG with six annotation strips (SCCmec, agr group, PVL, TSST, mecA, outbreak cluster) and tip points coloured by spa type.

#### System requirements (one-time, Fedora/RHEL)

```bash
sudo dnf install -y \
    cairo-devel fontconfig-devel freetype-devel \
    harfbuzz-devel fribidi-devel libpng-devel \
    libcurl-devel libxml2-devel libuv-devel pandoc
```

On Ubuntu/Debian:

```bash
sudo apt-get install -y \
    libcairo2-dev libfontconfig1-dev libfreetype-dev \
    libharfbuzz-dev libfribidi-dev libpng-dev \
    libcurl4-openssl-dev libxml2-dev libuv1-dev pandoc
```

R packages are installed automatically to `~/.R/library` on first run.

#### Usage

```bash
Rscript bin/plot_tree.R results/cc97 cc97
# → cc97_tree.pdf  (14 × ~26 in, scales with isolate count)
# → cc97_tree.svg  (same dimensions, text as SVG elements — editable in Inkscape)
```

Failures are non-fatal: the script exits 0 and writes a `PLOT_TREE_FAILED` file so a wrapping Nextflow process is never blocked.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Input/Output** | | |
| `--input` | | Samplesheet CSV (`--cc`/`--st` fetch mode bypasses this) |
| `--outdir` | `results` | Output directory |
| `--metadata` | | Sample metadata TSV (PHA4GE schema; provides dates for BEAST2) |
| `--antibiogram` | | Antibiogram CSV (NCBI long format) |
| **Public data fetch (--cc / --st mode)** | | |
| `--cc` | | Clonal complex to fetch (e.g. `97`). Expands to all member STs via PubMLST |
| `--st` | | Sequence type(s), comma-separated (e.g. `97` or `97,1153`) |
| `--fetch_source` | `all` | Metadata sources: `all`, `pubmlst`, `ncbi`, `ena` |
| `--fetch_country` | | Filter by country (e.g. `Saudi Arabia`) |
| `--fetch_continent` | | Filter by continent: `europe`, `asia`, `africa`, `americas` |
| `--fetch_host` | | Filter by host species (e.g. `Homo sapiens`) |
| `--fetch_year_from` | | Earliest isolation year |
| `--fetch_year_to` | | Latest isolation year |
| `--fetch_max` | | Cap number of downloads (useful for pilot runs) |
| **Assembly** | | |
| `--genome_size` | `2800000` | Expected genome size for read subsampling |
| `--target_depth` | `100` | Target coverage for Rasusa |
| **QC** | | |
| `--min_completeness` | `90` | CheckM2 minimum completeness (%) |
| `--max_contamination` | `5` | CheckM2 maximum contamination (%) |
| `--skip_qc_gate` | `false` | Skip CheckM2 quality filtering |
| **Typing** | | |
| `--sccmec_viz` | `false` | Generate SCCmec SVG/HTML element maps |
| **Phylogenetics** | | |
| `--phylo_method` | `panaroo` | `panaroo` (pangenome) or `snippy` (reference-based) |
| `--tree_builder` | `iqtree` | `iqtree` or `fasttree` |
| `--reference` | | Reference genome for snippy (GenBank format) |
| `--panaroo_clean` | `moderate` | Panaroo clean mode: `strict`, `moderate`, `sensitive` |
| `--panaroo_threshold` | `0.95` | Core genome threshold |
| `--snippy_mincov` | `10` | Snippy minimum read depth |
| `--snippy_minqual` | `100` | Snippy minimum mapping quality |
| **BEAST2 Bayesian phylodynamics** | | |
| `--use_beast` | `false` | Enable BEAST2 time-calibrated Bayesian phylogeny (requires `--phylo_method snippy` or `both`) |
| `--beast_chain_length` | `10000000` | MCMC chain length (steps) |
| `--beast_log_every` | `1000` | Log frequency (sample every N steps) |
| **Clustering** | | |
| `--cluster_snp_tiers` | `5,15,40` | SNP thresholds: direct, outbreak, related |
| `--cluster_cgmlst_tiers` | `10,24,50` | cgMLST thresholds |

## BEAST2 Bayesian phylodynamics

When `--use_beast` is passed, the pipeline runs a time-calibrated Bayesian phylogenetic analysis on the core SNP alignment produced by snippy-core. This is an optional step that runs after Phase 5 (Phylogeny) and requires `--phylo_method snippy` or `--phylo_method both`.

**What it does:**

- Reads `collection_date` from the metadata TSV (`--metadata`) and converts dates to decimal years for tip-date calibration
- If no metadata is provided, or the `collection_date` column is absent, BEAST2 runs without date calibration (unrooted clock model)
- Model: HKY + Gamma(4) site model, strict molecular clock, coalescent constant-size tree prior
- Produces a Maximum Clade Credibility (MCC) tree annotated with posterior node ages and HPD intervals

**How to enable:**

```bash
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --metadata metadata.tsv \
    --use_beast \
    --phylo_method snippy
```

**Optional tuning:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--beast_chain_length` | `10000000` | Total MCMC steps (increase for better convergence; check ESS > 200 in Tracer) |
| `--beast_log_every` | `1000` | Logging interval (lower = larger files but finer posterior sampling) |

**Outputs** (in `results/beast2/`):

| File | Description |
|------|-------------|
| `beast_mcc.tree` | MCC tree with posterior node age annotations (open in FigTree or R `treeio`) |
| `beast_run.log` | MCMC trace for convergence diagnostics (open in Tracer) |
| `beast_run.trees` | Full posterior tree distribution (10 % burn-in applied by TreeAnnotator) |

> **Note:** BEAST2 is computationally intensive. For large collections (>200 samples) consider reducing `--beast_chain_length` for a pilot run and checking ESS values in Tracer before a production run.

## Output

```
results/
├── trimgalore/         # Trimmed reads
├── rasusa/             # Subsampled reads
├── fastqc/             # Read quality reports
├── skesa/              # Assembled genomes
├── quast/              # Assembly QC
├── checkm2/            # Genome completeness/contamination
├── qc_gate/            # QC pass/fail summary
├── mlst/               # MLST types
├── spatyper/           # spa types
├── sccmec/             # SCCmec types (+SVG if --sccmec_viz)
├── agr_typing/         # agr groups
├── mash/               # Species confirmation
├── amrfinderplus/      # AMR genes + point mutations
├── abricate/           # resfinder + vfdb + plasmidfinder
├── kma/                # Read-based AMR
├── mob_recon/          # Plasmid reconstruction
├── plasmids/           # Plasmid summary
├── prokka/             # Genome annotations
├── panaroo/            # Pangenome analysis
├── iqtree/             # Phylogenetic tree
├── beast2/             # Bayesian phylogeny (only if --use_beast)
│   ├── beast_mcc.tree  # Maximum clade credibility annotated tree
│   ├── beast_run.log   # BEAST2 MCMC trace (ESS diagnostics in Tracer)
│   └── beast_run.trees # Posterior tree distribution (input for TreeAnnotator)
├── snpdists/           # SNP distance matrix
├── aggregated/         # Per-sample JSON reports + summary
├── clusters/           # Outbreak clusters + transmission pairs
├── report/             # Markdown run report
├── figures/            # Publication figures
├── multiqc/            # MultiQC report
└── pipeline_info/      # Execution reports
```

## Pipeline architecture

The pipeline separates **generic core** processes (QC, assembly, AMR, phylogeny, clustering) from **species-specific** analysis (typing, virulence, mutations). This is defined in:

```
subworkflows/local/core/          # Generic — works for any pathogen
subworkflows/local/species/       # Species-specific adapters
species/s_aureus.yml              # S. aureus descriptor
```

Adding a new pathogen requires only a species descriptor YAML + adapter subworkflows — the core pipeline stays unchanged.

## Tools and scripts

| Script | Purpose |
|--------|---------|
| `bin/staphit-metadata` | Metadata conversion (Vitek PDF/CSV, external XLSX) |
| `bin/pubmlst_fetch.py` | ST/CC-based public genome retrieval from PubMLST, NCBI SRA, and ENA |
| `bin/beast2_prep.py` | Generates BEAST2 2.7 XML from core SNP alignment + metadata dates |
| `bin/plot_tree.R` | Annotated phylogenetic tree (PDF + SVG) from pipeline outputs |
| `bin/staphit-aggregate` | Per-sample report aggregation |
| `bin/staphit-virulence` | Structured virulence profiling (PVL, TSST, IEC, operons) |
| `bin/staphit-mutations` | Point mutation extraction and phenotype prediction |
| `bin/staphit-cluster` | Tiered outbreak clustering with epi annotation |
| `bin/staphit-plasmids` | Plasmid profile parsing from MOB-suite |
| `bin/staphit-qc` | CheckM2 quality gate assessment |
| `bin/staphit-report` | Automated run report generator |
| `bin/staphit-visualize` | Publication figure generation |
| `bin/staphit-dashboard` | Interactive Plotly Dash surveillance dashboard |
| `bin/staphit-watch` | FASTQ directory monitor for continuous surveillance |

## Submitting to nf-core

This pipeline follows nf-core conventions and is on the path to community submission. Current status:

### Requirements for nf-core listing

| Requirement | Status |
|-------------|--------|
| nf-core template | ✅ Created with nf-core create |
| nextflow_schema.json | ⚠️ Needs custom params added |
| nf-test for all processes | ⚠️ Pending |
| CI/CD (GitHub Actions) | ✅ Template CI exists |
| `-profile test` with bundled data | ⚠️ Needs test dataset |
| Documentation (usage.md, output.md) | ⚠️ Pending |
| nf-core lint passes | ⚠️ Pending fixes |
| Standard nf-core modules where available | ✅ 11 nf-core modules used |
| Conda environment per module | ⚠️ Pending |
| Code review by nf-core community | Not started |

### Submission process

1. **Join nf-core Slack** — [https://nf-co.re/join](https://nf-co.re/join), introduce the pipeline in `#new-pipelines`
2. **Request a pipeline repo** — nf-core creates `nf-core/staphit2` under their org
3. **Pass nf-core lint** — `nf-core pipelines lint` must pass with no errors
4. **Add test profile** — bundled test data that runs in <10 minutes on CI
5. **Write docs** — `docs/usage.md` and `docs/output.md` following nf-core format
6. **Code review** — community review of the pipeline code
7. **First release** — tagged release with CHANGELOG, listed on nf-co.re/pipelines

### What makes this pipeline unique for nf-core

- **Surveillance-focused**: outbreak clustering, transmission pair detection, watch mode
- **Lab data integration**: Vitek 2 PDF parser, clinical XLSX converter, PHA4GE metadata
- **Structured virulence**: PVL/TSST/IEC typing beyond raw gene lists
- **Genotype-phenotype concordance**: compare predicted resistance against observed MICs
- **Interactive dashboard**: Plotly Dash with linked tree + map + table views
- **Multi-species scaffold**: species descriptor pattern ready for E. coli, Klebsiella, TB

## Credits

alarawms/staphit2 was developed by alarawms at KAIMRC/KAUST.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
