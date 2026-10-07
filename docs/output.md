# alarawms/staphit2: Output

## Introduction

This document describes the output produced by the pipeline. All directories listed below are created under the results directory specified with `--outdir`. Paths shown are relative to that top-level directory.

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data through the following stages:

1. [Read QC and preprocessing](#read-qc-and-preprocessing) -- TrimGalore, Rasusa, FastQC
2. [Assembly and quality assessment](#assembly-and-quality-assessment) -- SKESA/SPAdes, QUAST, CheckM2, QC gate
3. [Typing](#typing) -- MLST, spaTyper, SCCmec, agr, Mash
4. [Antimicrobial resistance](#antimicrobial-resistance) -- AMRFinderPlus, ABRicate, KMA
5. [Plasmid analysis](#plasmid-analysis) -- MOB-suite, plasmid summary
6. [Genome annotation and phylogenetics](#genome-annotation-and-phylogenetics) -- Prokka, Panaroo/Snippy, IQ-TREE/FastTree, SNP-dists
7. [Aggregation and reporting](#aggregation-and-reporting) -- per-sample reports, combined summary, clusters, run report, figures
8. [Quality reports](#quality-reports) -- MultiQC, pipeline information

---

## Read QC and preprocessing

### fastp/

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `<sample>.fastp.json`: Trimming and filtering statistics (also read by MultiQC).
  - `<sample>.fastp.html`: Per-sample HTML report.

</details>

[fastp](https://github.com/OpenGene/fastp) removes adapters (auto-detected for paired-end reads), trims 3' bases below Q20 and drops reads shorter than 20 bp. It runs on short reads of `short` and `hybrid` samples; the trimmed reads go to Rasusa and every short-read step downstream.

### nanoplot/

<details markdown="1">
<summary>Output files</summary>

- `nanoplot/<sample>/`
  - `<sample>_NanoStats.txt`: Read count, N50, mean/median length and quality.
  - `<sample>_NanoPlot-report.html`: Interactive length/quality plots.

</details>

[NanoPlot](https://github.com/wdecoster/NanoPlot) summarises Oxford Nanopore reads of `hybrid` and `long` samples. The stats files are included in the MultiQC report.

### rasusa/

<details markdown="1">
<summary>Output files</summary>

- `rasusa/`
  - `*_1.fastq.gz`: Subsampled forward reads.
  - `*_2.fastq.gz`: Subsampled reverse reads.

</details>

[Rasusa](https://github.com/mbhall88/rasusa) randomly subsamples reads to a target coverage depth (default: 100x for a 2.8 Mb genome, controlled by `--target_depth` and `--genome_size`). Normalizing coverage prevents highly-sequenced samples from inflating assembly size or biasing downstream analysis.

### fastqc/

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: Interactive HTML quality report.
  - `*_fastqc.zip`: Zip archive containing the report, tab-delimited data, and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides per-base quality scores, GC content, adapter content, sequence duplication levels, and overrepresented sequences. Reports are generated on the trimmed and subsampled reads. For further reading see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

---

## Assembly and quality assessment

### skesa/

<details markdown="1">
<summary>Output files</summary>

- `skesa/`
  - `*.fasta`: Assembled scaffolds in FASTA format.

</details>

[SKESA](https://github.com/ncbi/SKESA) (Strategic K-mer Extension for Scrupulous Assemblies) is the default *de novo* assembler. It is optimized for bacterial genomes and produces conservative, high-quality assemblies. If SPAdes is selected instead, output appears in `spades/` with similar file structure.

### dragonflye/

<details markdown="1">
<summary>Output files</summary>

- `dragonflye/<sample>/`
  - `<sample>.dragonflye.log`: Full assembly and polishing log.
  - `<sample>.gfa`: Flye assembly graph (view in [Bandage](https://rrwick.github.io/Bandage/) to check circular chromosome and plasmids).

</details>

[Dragonflye](https://github.com/rpetit3/dragonflye) assembles `long` and `hybrid` samples with [Flye](https://github.com/fenderglass/Flye), polishes with [Racon](https://github.com/lbcb-sci/racon), then with the Illumina reads for hybrid samples ([Polypolish](https://github.com/rrwick/Polypolish)) or with [Medaka](https://github.com/nanoporetech/medaka) for long-only samples when `--medaka_model` is set. The contigs go into the same downstream steps as short-read assemblies. Samples whose assembly is below 500 kb are dropped before QC.

### quast/

<details markdown="1">
<summary>Output files</summary>

- `quast/`
  - `report.html`: Interactive HTML report with assembly metrics.
  - `report.tsv`: Tab-separated assembly statistics.
  - `transposed_report.tsv`: Transposed version (one sample per row) for easy downstream parsing.

</details>

[QUAST](http://quast.sourceforge.net/) evaluates assembly quality. Key metrics include total assembly length, number of contigs, largest contig, N50 (weighted median contig length), and GC content. The `transposed_report.tsv` file is particularly useful for multi-sample comparisons.

### checkm2/

<details markdown="1">
<summary>Output files</summary>

- `checkm2/`
  - `*_quality_report.tsv`: Per-sample completeness and contamination estimates.

</details>

[CheckM2](https://github.com/chklovski/CheckM2) uses machine-learning models to estimate genome completeness and contamination from assembled contigs. These values drive the QC gate: samples below `--min_completeness` (default 90%) or above `--max_contamination` (default 5%) are flagged as failed.

### qc_gate/

<details markdown="1">
<summary>Output files</summary>

- `qc_gate/`
  - `qc_report.tsv`: Per-sample QC metrics and pass/fail status.
  - `passed_samples.txt`: List of sample IDs that passed all quality thresholds.
  - `failed_samples.txt`: List of sample IDs that failed one or more thresholds.

</details>

The QC gate integrates metrics from QUAST and CheckM2 to produce a pass/fail decision for each sample. Only samples in `passed_samples.txt` proceed to typing, AMR detection, phylogenetics, and clustering. Skipping the gate is possible with `--skip_qc_gate true`, in which case all samples proceed regardless of quality.

---

## Typing

### mlst/

<details markdown="1">
<summary>Output files</summary>

- `mlst/`
  - `*.tsv`: Tab-separated MLST results with scheme, sequence type (ST), and allelic profile.

</details>

[mlst](https://github.com/tseemann/mlst) assigns a 7-locus multi-locus sequence type (ST) by scanning assembled contigs against the PubMLST *S. aureus* scheme. The ST is a key epidemiological marker for lineage classification.

### spatyper/

<details markdown="1">
<summary>Output files</summary>

- `spatyper/`
  - `*.tsv`: spa type assignment for each sample.

</details>

[spaTyper](https://github.com/HCGB-IGTP/spaTyper) determines the *spa* type based on the short-sequence repeat region of the staphylococcal protein A gene. spa typing provides finer resolution than MLST for outbreak investigations and is the most widely used single-locus typing method for *S. aureus*.

### sccmec/

<details markdown="1">
<summary>Output files</summary>

- `sccmec/<sample>/`
  - `<sample>_sccmec.tsv`: SCCmec type call.
  - `<sample>_sccmec.json`: Full typer output (ccr/mec complexes, candidates, scores, contigs).
  - `<sample>_sccmec_elements.csv`: Every detected element (mec, ccr, IS, orfX) with coordinates.
  - `<sample>_sccmec_map.svg`, `<sample>_sccmec_report.html` (optional): Element maps (`--sccmec_viz true`).

</details>

The SCCmec typer classifies the staphylococcal cassette chromosome *mec* element, which carries the *mecA*/*mecC* gene conferring methicillin resistance. SCCmec types (I--XIII) provide insight into the evolutionary origin of resistance (hospital-associated vs. community-associated lineages). Enable visual maps with `--sccmec_viz true`.

SCCmec columns in `summary/combined_summary.tsv`:

| Column | Meaning |
|--------|---------|
| `sccmec_type` | Typer call as reported (may include `(best-fit)` with `--sccmec_best_fit`) |
| `sccmec_iwg` | IWG-SCC designation: ccr complex + mec class (e.g. `IV(2B)`, `V(5C2&5)`); composite and tandem elements keep both ccr complexes |
| `sccmec_group` | Harmonized group used for figure colours |
| `sccmec_type_candidates` | Other types consistent with the detected elements |
| `sccmec_assembly_limited` | `true` when the assembly splits the cassette across contigs and the type could not be closed from the assembly alone |
| `sccmec_typing_mode` | `assembly` or `reads` (the read fallback rescued a split cassette) |
| `sccmec_cassette` | Ordered elements of the cassette (orfX → mec → ccr) |

### agr_typing/

<details markdown="1">
<summary>Output files</summary>

- `agr_typing/`
  - `*.tsv`: agr group (I--IV) assignment per sample.

</details>

The accessory gene regulator (*agr*) system controls quorum-sensing and virulence factor expression. The agr group (I through IV) correlates with the spectrum of secreted toxins and has been associated with clinical outcomes.

### mash/

<details markdown="1">
<summary>Output files</summary>

- `mash/`
  - `*.tsv`: Top hits from the Mash reference sketch database with genomic distances.

</details>

[Mash](https://github.com/marbl/Mash) provides rapid species confirmation by computing MinHash-based genomic distances against a reference database. Samples with unexpectedly large distances to *S. aureus* references may indicate contamination or mis-labelled isolates.

---

## Antimicrobial resistance

### amrfinderplus/

<details markdown="1">
<summary>Output files</summary>

- `amrfinderplus/`
  - `*.tsv`: AMR gene and point mutation results in NCBI AMRFinderPlus tabular format.

</details>

[AMRFinderPlus](https://github.com/ncbi/amr) identifies acquired AMR genes, stress response genes, and chromosomal point mutations from assembled contigs. When run with the `--organism Staphylococcus_aureus` flag (set automatically by the pipeline), it additionally reports clinically relevant point mutations in genes such as *gyrA*, *parC* (fluoroquinolone resistance), *rpoB* (rifampicin resistance), and *fusA* (fusidic acid resistance). Rows with `Subtype=POINT` in the output indicate point mutations rather than acquired genes.

### abricate/

<details markdown="1">
<summary>Output files</summary>

- `abricate/`
  - `*_resfinder.tsv`: Acquired AMR genes detected against the ResFinder database.
  - `*_vfdb.tsv`: Virulence factor genes detected against the VFDB database.
  - `*_plasmidfinder.tsv`: Plasmid replicon sequences detected against the PlasmidFinder database.

</details>

[ABRicate](https://github.com/tseemann/abricate) performs mass screening of assembled contigs against multiple curated databases. The pipeline runs three database searches in parallel:

- **ResFinder**: Acquired resistance genes, complementing AMRFinderPlus with an alternative gene catalogue.
- **VFDB**: Virulence factor genes including toxins, adhesins, and immune evasion factors.
- **PlasmidFinder**: Plasmid replicon typing to identify the incompatibility group of plasmids present.

### kma/

<details markdown="1">
<summary>Output files</summary>

- `kma/`
  - `*.res`: Read-based AMR detection results (template coverage, depth, identity).
  - `*.aln`: Consensus alignments to matched resistance gene templates.

</details>

[KMA](https://bitbucket.org/genomicepidemiology/kma/) maps raw reads directly against resistance gene databases, providing an assembly-independent confirmation of AMR determinants. This is particularly useful for detecting genes present on low-abundance plasmids or in mixed cultures that may not assemble well.

---

## Plasmid analysis

### mob_recon/

<details markdown="1">
<summary>Output files</summary>

- `mob_recon/`
  - `contig_report.txt`: Assignment of each assembled contig to chromosome or a reconstructed plasmid.
  - `mobtyper_results.txt`: Plasmid typing results including predicted mobility (conjugative, mobilizable, non-mobilizable), replicon type, and relaxase type.
  - `plasmid_*.fasta`: Reconstructed plasmid sequences (one FASTA per predicted plasmid).
  - `chromosome.fasta`: Contigs assigned to the chromosome.

</details>

[MOB-suite](https://github.com/phac-nml/mob-suite) reconstructs individual plasmids from assembled genomes, assigns replicon and mobility types, and classifies each contig as chromosomal or plasmid-borne. This enables tracking of plasmid-mediated resistance transfer between isolates.

### plasmids/

<details markdown="1">
<summary>Output files</summary>

- `plasmids/`
  - `*.tsv`: Merged plasmid summary per sample (plasmid count, replicon types, mobility classification).

</details>

A pipeline-specific summary that consolidates MOB-suite results into a per-sample table reporting the number of predicted plasmids, their replicon types, and mobility status.

---

## Genome annotation and phylogenetics

### prokka/

<details markdown="1">
<summary>Output files</summary>

- `prokka/`
  - `*.gff`: Annotated genome in GFF3 format.
  - `*.gbk`: Annotated genome in GenBank format.
  - `*.faa`: Predicted protein sequences (FASTA).
  - `*.ffn`: Predicted nucleotide coding sequences (FASTA).
  - `*.txt`: Annotation summary statistics.

</details>

[Prokka](https://github.com/tseemann/prokka) performs rapid whole-genome annotation. The GFF files are consumed by Panaroo for pangenome analysis. GenBank and protein FASTA files are useful for downstream comparative analyses.

### panaroo/

<details markdown="1">
<summary>Output files</summary>

- `panaroo/`
  - `core_gene_alignment.aln`: Concatenated core genome alignment (FASTA) used for tree building.
  - `gene_presence_absence.csv`: Gene presence/absence matrix across all samples.
  - `summary_statistics.txt`: Pangenome summary (core, soft-core, shell, cloud gene counts).

</details>

[Panaroo](https://github.com/gtonkinhill/panaroo) builds a pangenome graph from Prokka GFF annotations, identifies the core genome, and produces a concatenated core gene alignment. This is the default phylogenetic input method (`--phylo_method panaroo`). The `--panaroo_clean` parameter controls graph-cleaning stringency and `--panaroo_threshold` sets the fraction of samples required for a gene to be considered core.

When `--phylo_method snippy` is used instead, the output appears in `snippy_core/` and contains a reference-based core SNP alignment produced by [Snippy](https://github.com/tseemann/snippy).

### iqtree/

<details markdown="1">
<summary>Output files</summary>

- `iqtree/`
  - `*.treefile`: Maximum-likelihood phylogenetic tree in Newick format.
  - `*.log`: IQ-TREE log with model selection and bootstrap details.
  - `*.iqtree`: Full IQ-TREE report.

</details>

[IQ-TREE](http://www.iqtree.org/) infers a maximum-likelihood phylogenetic tree from the core genome alignment with automatic model selection (ModelFinder) and ultrafast bootstrap support values. This is the default tree builder (`--tree_builder iqtree`).

When `--tree_builder fasttree` is selected, output appears in `fasttree/` containing a Newick tree file produced by [FastTree](http://www.microbesonline.org/fasttree/) (approximate maximum-likelihood, faster but less accurate).

### snpdists/

<details markdown="1">
<summary>Output files</summary>

- `snpdists/`
  - `*.tsv`: Pairwise SNP distance matrix (tab-separated, symmetric).

</details>

[SNP-dists](https://github.com/tseemann/snp-dists) computes pairwise SNP distances from the core genome alignment. This matrix is consumed by the clustering module and is also useful for manual inspection of isolate relatedness.

### itol/

<details markdown="1">
<summary>Output files</summary>

- `itol/`
  - `01_sccmec.txt`: SCCmec type — colour strip.
  - `02_mlst.txt`: MLST sequence type — colour strip.
  - `03_agr.txt`: *agr* group (I–IV) — colour strip.
  - `04_pvl.txt`: PVL status (positive/negative) — binary symbol.
  - `05_spa.txt`: *spa* type label — text dataset.
  - `06_st_label.txt`: ST label — text dataset.
  - `07_amr_class.txt`: AMR drug-class presence — binary dataset.
  - `08_virulence.txt`: Selected virulence genes — binary dataset.
  - `09_origin.txt`: Infection origin (HA/CA/LA) — colour strip.
  - `10_country.txt`: Country of isolation — colour strip.
  - `11_host.txt`: Host (human/animal/environmental) — colour strip.
  - `12_source.txt`: Sample source (blood/wound/screen, etc.) — colour strip.
  - `13_hospital.txt`: Hospital/facility — text dataset.
  - `14_year.txt`: Collection year — colour strip (gradient).
  - `15_gender.txt`: Patient gender — colour strip.
  - `16_patient_type.txt`: Patient type (inpatient/outpatient) — colour strip.
  - `17_city.txt`: City of collection — text dataset.
  - `18_region.txt`: Administrative region — colour strip.

</details>

The iTOL annotation export (`bin/export_itol.py`) generates ready-to-upload annotation files for [iTOL (Interactive Tree of Life)](https://itol.embl.de/). Each file corresponds to one annotation track and is formatted according to the iTOL dataset specification for its type (DATASET_COLORSTRIP, DATASET_TEXT, DATASET_BINARY, DATASET_SYMBOL).

The script draws annotation data from three sources, checked in priority order:

| Source | Contents | Applies to |
|--------|----------|-----------|
| `results/summary/combined_summary.tsv` | Pipeline typing output (ST, SCCmec, spa, agr, AMR, virulence, origin) | All samples that passed QC |
| `--pub_meta` TSV | Manually curated public metadata (country, year, ST, spa) | Public ENA/SRA accessions |
| `--local_meta` CSV | Local cohort metadata (hospital, city, region, gender, patient type) | Local sequenced samples |

Only tracks with at least one data value are written to disk. Tracks whose data is entirely absent (e.g. hospital for a public-data-only tree) produce no file.

**Generating iTOL annotations:**

```bash
python bin/export_itol.py \
    results/ \
    my_run \
    --tree results/iqtree/core.treefile \
    cc97_metadata.tsv \
    assets/sample_metadata.csv
```

The `--tree` flag is strongly recommended: it uses the treefile as the authoritative source of node IDs, resolves mismatches between metadata keys and tree labels, and automatically excludes the Snippy reference genome node (`Reference`) which would otherwise appear as an unannotated leaf.

---

## Databases

<details markdown="1">
<summary>Output files</summary>

- `databases/`
  - `checkm2/uniref100.KO.1.dmnd`: CheckM2 database (reuse with `--checkm2_db`).
  - `amrfinderplus/amrfinderdb.tar.gz`: AMRFinderPlus database (reuse with `--amrfinder_db`).
  - `mash_refseq/refseq.genomes.msh`: Mash RefSeq sketch for species QC (reuse with `--mash_db`).
  - `mob_suite/mob_db/`: MOB-suite databases.
  - `resfinder/resfinder_db/`: ResFinder FASTA files at commit `--resfinder_db_commit`.

</details>

Database versions or checksums are written to `pipeline_info/staphit2_software_mqc_versions.yml` with the tool versions.

---

## Aggregation and reporting

### aggregated/

<details markdown="1">
<summary>Output files</summary>

- `aggregated/`
  - `*_report.json`: Comprehensive per-sample JSON report.
  - `*_summary.tsv`: Tabular per-sample summary.

</details>

The aggregator (`bin/staphit-aggregate`) merges all upstream results into a single structured JSON report per sample. Each report contains:

- **Typing**: ST, spa type, SCCmec type, agr group.
- **Resistance**: Acquired AMR genes, drug class predictions, and resistance profile summary.
- **Point mutations**: Chromosomal mutations (gyrA, parC, rpoB, fusA, etc.) with predicted drug class.
- **Virulence profile**: PVL status, TSST-1, IEC type, enterotoxin count, exfoliatins, biofilm/hemolysin operon completeness.
- **Plasmids**: Count, replicon types, mobility.
- **AST profile**: Antibiogram concordance (when `--antibiogram` is provided).
- **Infection origin**: HA-/CA-/LA-MRSA classification (when `--metadata` is provided).

### summary/

<details markdown="1">
<summary>Output files</summary>

- `summary/`
  - `combined_summary.tsv`: Merged tab-separated summary of all samples in the run.

</details>

A single wide-format table with one row per sample and columns for every typing result, resistance gene, virulence marker, and QC metric. This file is designed for direct import into spreadsheet software or downstream statistical analysis.

### clusters/

<details markdown="1">
<summary>Output files</summary>

- `clusters/`
  - `clusters.tsv`: Cluster assignment for each sample at three tiers (direct transmission, outbreak, related).
  - `cluster_report.json`: Structured JSON report with cluster membership, sizes, and metadata annotations.
  - `transmission_pairs.tsv`: Pairs of samples within the direct transmission threshold with their SNP distances.

</details>

The clustering module (`bin/staphit-cluster`) applies single-linkage clustering to the SNP distance matrix at three configurable thresholds (default: 5, 15, and 40 SNPs for *S. aureus*):

| Tier | Default SNP threshold | Interpretation |
|------|----------------------|----------------|
| Tier 1 | <=5 SNPs | Direct/recent transmission |
| Tier 2 | <=15 SNPs | Part of the same outbreak |
| Tier 3 | <=40 SNPs | Epidemiologically related |

When metadata is provided, cluster reports are annotated with epidemiological context (collection dates, locations, infection origin).

### report/

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `run_report.md`: Automated Markdown surveillance report.

</details>

The report generator (`bin/staphit-report`) produces a narrative Markdown report summarizing the run: sample counts, QC pass rates, lineage distribution, resistance patterns, cluster summary, and key findings. This report is suitable for direct sharing with infection control teams.

### figures/

<details markdown="1">
<summary>Output files</summary>

- `figures/`
  - `*.png`: Publication-ready figures in PNG format.

</details>

The visualization module (`bin/staphit-visualize`) generates figures including:

- Phylogenetic tree annotated with typing and resistance data.
- Resistance heatmap across samples and drug classes.
- Cluster network diagrams.
- SNP distance distributions.

All figures are rendered at publication resolution (300 DPI).

---

## Quality reports

### multiqc/

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Standalone HTML report aggregating QC metrics from all tools.
  - `multiqc_data/`: Directory containing parsed statistics in machine-readable format.
  - `multiqc_plots/`: Directory containing static plot images from the report.

</details>

[MultiQC](http://multiqc.info) aggregates quality metrics from FastQC, TrimGalore, QUAST, and other supported tools into a single interactive HTML report. This provides a rapid overview of sequencing quality and assembly statistics across all samples. For more information about interpreting MultiQC reports, see the [MultiQC documentation](http://multiqc.info).

### pipeline_info/

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`: Nextflow execution report with task-level resource usage.
  - `execution_timeline_*.html`: Timeline visualization of task execution.
  - `execution_trace_*.txt`: Tab-separated trace file with per-task metrics (CPU, memory, duration).
  - `pipeline_dag_*.html`: Directed acyclic graph (DAG) of the pipeline workflow.
  - `params.json`: Parameters used for the pipeline run.
  - `software_versions.yml`: Versions of all software tools used.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) automatically generates execution reports, timelines, and trace files for every run. These are essential for troubleshooting failed tasks, benchmarking resource usage, and ensuring reproducibility. The `params.json` file records the exact parameters used, and `software_versions.yml` lists all tool versions for citation purposes.
