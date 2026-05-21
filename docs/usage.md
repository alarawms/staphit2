# alarawms/staphit2: Usage

## Introduction

**alarawms/staphit2** is an nf-core-compatible Nextflow pipeline for comprehensive MRSA (*Staphylococcus aureus*) genomic surveillance. Starting from raw Illumina paired-end reads, it performs quality control, *de novo* assembly, multi-locus typing (MLST, spa, SCCmec, agr), antimicrobial resistance gene detection and point mutation calling, structured virulence profiling, plasmid reconstruction, core genome phylogenetics, pairwise SNP distances, tiered outbreak clustering, and automated reporting with publication-ready figures.

The pipeline is built on a species-agnostic core architecture. Species-specific analysis steps (typing tools, virulence categories, outbreak thresholds) are driven by a descriptor file (`species/s_aureus.yml`), enabling future extension to additional pathogens without modifying the core workflow.

## Samplesheet input

You will need to create a samplesheet CSV describing the samples you want to analyse. Specify it with the `--input` parameter:

```bash
--input samplesheet.csv
```

The file must be comma-separated with a header row and at least three columns. Each row represents one sample (paired-end reads).

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
MRSA_001,/data/MRSA_001_R1.fastq.gz,/data/MRSA_001_R2.fastq.gz
MRSA_002,/data/MRSA_002_R1.fastq.gz,/data/MRSA_002_R2.fastq.gz
MRSA_003,/data/MRSA_003_R1.fastq.gz,/data/MRSA_003_R2.fastq.gz
```

| Column   | Description |
|----------|-------------|
| `sample` | Unique sample identifier. Must not contain spaces. If the same identifier appears on multiple rows the pipeline will concatenate the raw reads before downstream analysis (useful for samples sequenced across multiple lanes). |
| `fastq_1` | Full path to the forward-read FASTQ file. Must be gzipped (`.fastq.gz` or `.fq.gz`). |
| `fastq_2` | Full path to the reverse-read FASTQ file. Must be gzipped (`.fastq.gz` or `.fq.gz`). |

An [example samplesheet](../assets/samplesheet.csv) is provided with the pipeline.

### Multiple runs of the same sample

If a sample was sequenced more than once (e.g. to increase depth), use the same `sample` name on each row. The pipeline will concatenate reads before any downstream processing:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
MRSA_001,MRSA_001_L001_R1.fastq.gz,MRSA_001_L001_R2.fastq.gz
MRSA_001,MRSA_001_L002_R1.fastq.gz,MRSA_001_L002_R2.fastq.gz
```

## Metadata input (optional)

The pipeline accepts two optional metadata files that enable genotype--phenotype concordance analysis, epidemiological annotation of outbreak clusters, and richer surveillance reports.

### Sample metadata (`--metadata`)

A CSV following the [PHA4GE](https://pha4ge.org/) contextual data schema. At minimum it must contain a `sample_id` column that matches the `sample` column in the samplesheet.

```csv title="sample_metadata.csv"
sample_id,organism,collection_date,geo_loc_country,geo_loc_region,host,isolation_source,purpose_of_sampling
MRSA_001,Staphylococcus aureus,2025-11-01,Saudi Arabia,Riyadh,Homo sapiens,nasal,diagnostic testing
MRSA_002,Staphylococcus aureus,2025-11-03,Saudi Arabia,Riyadh,Homo sapiens,blood,diagnostic testing
```

Commonly used fields include:

| Field | Description |
|-------|-------------|
| `sample_id` | Must match samplesheet `sample` column |
| `organism` | Species name |
| `collection_date` | ISO 8601 date (YYYY-MM-DD) |
| `geo_loc_country` | Country of collection |
| `geo_loc_region` | Sub-national region |
| `host` | Host organism (e.g. `Homo sapiens`) |
| `isolation_source` | Body site or environmental source |
| `purpose_of_sampling` | Reason for sampling |
| `infection_origin` | HA-MRSA, CA-MRSA, or LA-MRSA |

### Antibiogram (`--antibiogram`)

A CSV in NCBI BioSample Antibiogram long format, with one row per sample--antibiotic combination:

```csv title="antibiogram.csv"
sample_id,antibiotic,resistance_phenotype,measurement,measurement_units,measurement_sign,laboratory_typing_method,testing_standard
MRSA_001,Vancomycin,susceptible,1,mg/L,<=,MIC,CLSI
MRSA_001,Oxacillin,resistant,4,mg/L,>=,MIC,CLSI
MRSA_002,Vancomycin,susceptible,0.5,mg/L,<=,MIC,CLSI
```

| Field | Description |
|-------|-------------|
| `sample_id` | Must match samplesheet `sample` column |
| `antibiotic` | Antimicrobial agent name |
| `resistance_phenotype` | `susceptible`, `intermediate`, or `resistant` |
| `measurement` | MIC value or zone diameter |
| `measurement_units` | `mg/L` (MIC) or `mm` (disk diffusion) |
| `measurement_sign` | `<=`, `=`, or `>=` |
| `laboratory_typing_method` | `MIC`, `DISK`, or `ETEST` |
| `testing_standard` | `CLSI` or `EUCAST` |

### Generating metadata from lab data

The pipeline ships helper scripts in `bin/` for converting common lab formats into the required CSVs:

```bash
# Convert Vitek 2 PDF reports to NCBI antibiogram format
python bin/staphit-metadata convert --from-vitek-pdf vitek_pdfs/ -o antibiogram.csv

# Convert a clinical Excel spreadsheet to PHA4GE metadata
python bin/staphit-metadata convert --from-external-xlsx clinical.xlsx \
    --metadata sample_metadata.csv -o enriched_metadata.csv

# Search SRA for public S. aureus data and download
python bin/staphit-fetch search --organism "Staphylococcus aureus" --country "Saudi Arabia" -o results.tsv
python bin/staphit-fetch download results.tsv --output-dir fetched/
```

An [external metadata template](../assets/external_metadata_template.csv) is provided for reference.

## Running the pipeline

### Basic run

```bash
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

### Full surveillance run (with metadata and antibiogram)

```bash
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --metadata sample_metadata.csv \
    --antibiogram antibiogram.csv
```

### Resume after failure or adding samples

```bash
nextflow run alarawms/staphit2 \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    -resume
```

### Using a params file

Rather than specifying every flag on the command line, you can place parameters in a YAML file:

```yaml title="params.yaml"
input: 'samplesheet.csv'
outdir: 'results'
metadata: 'sample_metadata.csv'
antibiogram: 'antibiogram.csv'
phylo_method: 'panaroo'
tree_builder: 'iqtree'
target_depth: 100
```

```bash
nextflow run alarawms/staphit2 -profile docker -params-file params.yaml
```

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

Note that the pipeline will create the following files in your working directory:

```
work/               # Nextflow working files
<OUTDIR>/           # Finished results (defined with --outdir)
.nextflow_log       # Log file from Nextflow
```

## Pipeline parameters

### Input/Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Path to samplesheet CSV |
| `--outdir` | (required) | Path to output directory |
| `--metadata` | `null` | Path to sample metadata CSV (PHA4GE schema) |
| `--antibiogram` | `null` | Path to antibiogram CSV (NCBI long format) |
| `--genome` | `null` | iGenomes genome key (not typically used) |
| `--reference` | `null` | Reference genome in GenBank format (required when `--phylo_method snippy`) |

### Assembly and read processing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genome_size` | `2800000` | Expected genome size in bp (used by Rasusa for subsampling) |
| `--target_depth` | `100` | Target coverage depth for Rasusa read subsampling |

### Quality control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_completeness` | `90` | Minimum CheckM2 completeness (%) to pass QC gate |
| `--max_contamination` | `5` | Maximum CheckM2 contamination (%) to pass QC gate |
| `--skip_qc_gate` | `false` | Skip the CheckM2 quality gate and process all samples |
| `--species_ani_threshold` | `95.0` | Minimum ANI (%) to NCTC 8325 for *S. aureus* species confirmation |
| `--skip_species_qc` | `false` | Skip fastANI species confirmation step |
| `--assembler` | `skesa` | Assembly tool: `skesa` (default) or `spades` |

### Typing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sccmec_viz` | `false` | Generate SVG/HTML visual maps of SCCmec cassette elements |

### Phylogenetics

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--phylo_method` | `panaroo` | Core genome method: `panaroo` (pangenome-based) or `snippy` (reference-based) |
| `--tree_builder` | `iqtree` | Tree inference: `iqtree` (ML, slower, more accurate) or `fasttree` (approximate ML, faster) |
| `--panaroo_clean` | `moderate` | Panaroo graph-cleaning stringency: `strict`, `moderate`, or `sensitive` |
| `--panaroo_threshold` | `0.95` | Fraction of samples a gene must appear in to be considered core |
| `--panaroo_aligner` | `mafft` | Alignment tool used by Panaroo |
| `--snippy_mincov` | `10` | Minimum read depth for Snippy variant calls |
| `--snippy_minqual` | `100` | Minimum mapping quality for Snippy variant calls |

### Clustering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cluster_snp_tiers` | `15,50,200` | Core gene SNP thresholds for direct transmission, outbreak, and related tiers (calibrated for Panaroo alignment; use `5,15,40` for WGS reference-mapped distances) |
| `--cluster_cgmlst_tiers` | `10,24,50` | Comma-separated cgMLST allelic distance thresholds for the same three tiers |

### MultiQC

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--multiqc_config` | `null` | Path to a custom MultiQC config YAML |
| `--multiqc_title` | `null` | Custom title for the MultiQC report |
| `--multiqc_logo` | `null` | Path to a custom logo for the MultiQC report |

## Profiles

Use `-profile` to select a software packaging method. Multiple profiles can be combined (e.g. `-profile test,docker`):

| Profile | Description |
|---------|-------------|
| `docker` | Run all tools via [Docker](https://docker.com/) containers (recommended) |
| `singularity` | Run via [Singularity](https://sylabs.io/docs/) containers (recommended for HPC) |
| `apptainer` | Run via [Apptainer](https://apptainer.org/) containers |
| `conda` | Run via [Conda](https://conda.io/) environments (last resort) |
| `mamba` | Conda with [Mamba](https://mamba.readthedocs.io/) solver |
| `podman` | Run via [Podman](https://podman.io/) containers |
| `shifter` | Run via [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/) containers |
| `charliecloud` | Run via [Charliecloud](https://charliecloud.io/) containers |
| `wave` | Enable [Wave](https://seqera.io/wave/) containers (use with another profile) |
| `test` | Minimal test dataset; runs automatically without additional input |
| `test_full` | Full-size test dataset for complete validation |

> [!IMPORTANT]
> We highly recommend Docker or Singularity for full reproducibility. Use Conda only when containers are not available.

The pipeline also dynamically loads institutional profiles from [nf-core/configs](https://github.com/nf-core/configs) at runtime.

## Updating the pipeline

When you first run the pipeline, Nextflow caches the code locally. To pull the latest version:

```bash
nextflow pull alarawms/staphit2
```

## Reproducibility

Pin a specific release when running production analyses:

```bash
nextflow run alarawms/staphit2 -r 1.0.0 -profile docker --input samplesheet.csv --outdir results
```

The version number is recorded in the MultiQC report and in `results/pipeline_info/` execution reports, ensuring traceability.

To repeat a run with identical settings, reuse a [params file](#using-a-params-file). When sharing params files (e.g. as supplementary material), remove any cluster-specific paths or institutional profile references.

## Clonal complex surveillance

This section describes the full workflow for building a curated phylogenetic tree from a mixture of locally sequenced isolates and public data from a target clonal complex (e.g. CC97).

### 1. Curate confirmed accessions

Start from a list of ENA/SRA accessions and cross-reference them against a metadata table that contains sequence-type (ST) information. Only accessions with a confirmed CC97 sequence type (ST97 or a known single-locus variant) should enter the pipeline. The script `bin/filter_cc97_samplesheet.py` handles the post-run filtering step, but manual pre-curation of the accession list is also recommended to avoid downloading large volumes of non-CC97 data.

CC97 includes the following sequence types (STs):

| ST | Relationship to ST97 |
|----|---------------------|
| 97 | Founder |
| 1153 | Single-locus variant (SLV) |
| 1465 | SLV |
| 3187 | SLV |
| 8064 | SLV |
| 2996 | SLV |
| 5435 | SLV |
| 3528 | SLV |
| 2458 | SLV |
| 3009 | SLV |
| 7570 | SLV |
| 3636 | SLV |
| 5264 | SLV |

Accessions from general surveillance studies that do not report ST metadata cannot be confirmed as CC97 and should be excluded until after typing.

### 2. Download public reads with fetchngs

Create a CSV with one ENA/SRA accession per line (header: `id`) and run nf-core/fetchngs:

```bash
nextflow run nf-core/fetchngs \
    --input cc97_accessions_to_fetch.csv \
    --outdir cc97_fetchngs \
    -profile docker \
    -resume
```

Downloaded FASTQ files appear in `cc97_fetchngs/fastq/`. Each accession also produces a `*.runinfo_ftp.tsv` metadata file in `cc97_fetchngs/metadata/` containing study information and available fields from ENA.

### 3. Build the combined samplesheet

Combine locally sequenced samples with the downloaded public samples into a single samplesheet. Local samples use internal IDs (e.g. `ID00001`); public samples use their ENA run accession (e.g. `ERR1234567`).

```csv
sample,fastq_1,fastq_2
ID00001,/path/to/ID00001_R1.fastq.gz,/path/to/ID00001_R2.fastq.gz
ERR1234567,cc97_fetchngs/fastq/ERR1234567_1.fastq.gz,cc97_fetchngs/fastq/ERR1234567_2.fastq.gz
```

### 4. Run the pipeline

Use `--phylo_method snippy` for CC97 collections. Panaroo (the default) can fail with `KeyError` on highly fragmented assemblies from older public datasets; Snippy-based alignment is more robust for mixed-quality inputs.

```bash
nextflow run alarawms/staphit2 \
    -profile docker \
    --input cc97_combined_samplesheet.csv \
    --outdir results \
    --phylo_method snippy \
    --reference assets/nctc8325.fasta \
    --metadata assets/sample_metadata.csv \
    -resume
```

> [!NOTE]
> The Snippy reference genome node (`Reference`) will appear as a leaf in the tree produced by `snippy_core`. It is automatically excluded from iTOL annotations by `export_itol.py` when `--tree` is provided.

### 5. Filter to confirmed CC97 by MLST

After the run completes, use `bin/filter_cc97_samplesheet.py` to remove any samples that did not type as CC97:

```bash
python bin/filter_cc97_samplesheet.py \
    results/ \
    cc97_combined_samplesheet.csv \
    cc97_filtered_samplesheet.csv
```

This script reads `results/summary/combined_summary.tsv`, keeps all samples with a CC97 sequence type, keeps local IDs (`ID*`) that had no MLST result (QC failures within the local cohort are reviewed separately), and drops public accessions that failed assembly QC or typed as non-CC97. A summary is printed showing counts and the ST distribution of retained samples.

### 6. Rebuild the tree on the filtered set

Re-run the pipeline on the filtered samplesheet. Nextflow `-resume` will reuse all cached typing results; only the alignment and tree steps will rerun:

```bash
nextflow run alarawms/staphit2 \
    -profile docker \
    --input cc97_filtered_samplesheet.csv \
    --outdir results \
    --phylo_method snippy \
    --reference assets/nctc8325.fasta \
    --metadata assets/sample_metadata.csv \
    -resume
```

### 7. Generate iTOL annotations

Export annotation tracks for the filtered tree:

```bash
python bin/export_itol.py \
    results/ \
    cc97 \
    --tree results/iqtree/core.treefile \
    cc97_metadata.tsv \
    assets/sample_metadata.csv
```

Annotation files are written to `results/itol/`. Upload the treefile and annotation files to [iTOL](https://itol.embl.de/) for interactive visualization. See the [itol/ output section](output.md#itol) for a description of each track.

**Annotation data sources for mixed local + public trees:**

| Annotation track | Local samples | Public (ENA) samples |
|-----------------|---------------|----------------------|
| ST, SCCmec, spa, agr, AMR, virulence, origin | `combined_summary.tsv` | `combined_summary.tsv` (if passed QC) |
| Country, year | `combined_summary.tsv` | `cc97_metadata.tsv` |
| Hospital, city, region, gender, patient type | `sample_metadata.csv` | Not available |

Public samples that failed assembly QC (CheckM2 completeness < 90%) pass Snippy alignment and therefore appear in the tree but have no typing data. Their annotation tracks will be blank for all fields except country and year (if present in the metadata TSV).

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-resume`

Reuse cached results from previous runs. Nextflow will only re-execute steps whose inputs have changed.

```bash
nextflow run alarawms/staphit2 -profile docker --input samplesheet.csv --outdir results -resume
```

### `-c`

Specify an additional Nextflow config file. Use this for resource tuning or infrastructure settings, **not** for pipeline parameters (use `-params-file` instead). See the [nf-core configuration docs](https://nf-co.re/docs/usage/configuration) for details.

## Custom configuration

### Resource requests

Default resource requests are defined in `conf/base.config`. Failed jobs are automatically retried with increased resources (up to 3 attempts). To override defaults, see the [nf-core resource tuning guide](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).

### Custom tool arguments

To pass additional arguments to a specific tool, override `ext.args` in a custom config. See the [nf-core tool arguments guide](https://nf-co.re/docs/usage/configuration#customising-tool-arguments).

## Running in the background

Use `screen`, `tmux`, or the Nextflow `-bg` flag to keep the pipeline running after you disconnect:

```bash
nextflow run alarawms/staphit2 -profile docker --input samplesheet.csv --outdir results -bg
```

## Nextflow memory requirements

If Nextflow consumes too much memory, add the following to your environment (e.g. `~/.bashrc`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
