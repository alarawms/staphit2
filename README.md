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

### Metadata tools

Convert lab data (Vitek PDFs, clinical spreadsheets) into pipeline-ready CSVs:

```bash
# Vitek 2 PDFs → antibiogram with full MIC data
python bin/staphit-metadata convert --from-vitek-pdf vitek_pdfs/ -o antibiogram.csv

# External clinical XLSX → enrich metadata
python bin/staphit-metadata convert --from-external-xlsx clinical.xlsx \
    --metadata sample_metadata.csv -o supplement.csv

# SRA search → download public data
python bin/staphit-fetch search --organism "Staphylococcus aureus" --country "Saudi Arabia" -o results.tsv
python bin/staphit-fetch download results.tsv --output-dir fetched/
```

### Interactive dashboard

```bash
pip install dash plotly pandas dash-bootstrap-components
python bin/staphit-dashboard results/ --host 0.0.0.0 --port 8050
```

Two modes: **Investigation** (linked tree + map + table) and **Analysis** (resistance, trends, virulence, clusters, plasmids, QC tabs).

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Input/Output** | | |
| `--input` | | Samplesheet CSV (required) |
| `--outdir` | `results` | Output directory |
| `--metadata` | | Sample metadata CSV (PHA4GE schema) |
| `--antibiogram` | | Antibiogram CSV (NCBI long format) |
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
| **Clustering** | | |
| `--cluster_snp_tiers` | `5,15,40` | SNP thresholds: direct, outbreak, related |
| `--cluster_cgmlst_tiers` | `10,24,50` | cgMLST thresholds |

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
| `bin/staphit-fetch` | SRA/ENA search with metadata mapping |
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
