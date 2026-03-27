# Staphit v2 — Migration & Build Plan

## Big Picture

Staphit v2 is an nf-core MRSA surveillance pipeline built on a **species-agnostic core**. The architecture separates generic bioinformatics (QC, assembly, AMR, phylogeny) from species-specific intelligence (typing, virulence, mutations, thresholds). Adding another pathogen later means writing a species descriptor + adapter scripts — the core doesn't change.

```
staphit2/
├── main.nf                          # Entry point
├── workflows/staphit2.nf            # Main workflow orchestrator
│
├── subworkflows/
│   ├── local/
│   │   ├── input_check.nf           # Samplesheet validation
│   │   ├── prepare_genome.nf        # Reference genome handling
│   │   ├── core/                    # ── GENERIC (any pathogen) ──
│   │   │   ├── qc.nf               # Trim + FastQC + CheckM2
│   │   │   ├── assembly.nf         # SPAdes/SKESA + QUAST
│   │   │   ├── amr.nf              # AMRFinderPlus + ABRicate + KMA
│   │   │   ├── phylogeny.nf        # Panaroo/Snippy → IQ-TREE/FastTree
│   │   │   ├── clustering.nf       # Distance-based outbreak detection
│   │   │   └── plasmids.nf         # MOB-suite
│   │   │
│   │   └── species/                 # ── SPECIES-SPECIFIC ──
│   │       └── s_aureus/
│   │           ├── typing.nf        # spa, SCCmec, agr (+ MLST, cgMLST)
│   │           ├── virulence.nf     # PVL, TSST, IEC, operons
│   │           └── mutations.nf     # gyrA, parC, rpoB point mutations
│   │
│   └── nf-core/                     # Imported nf-core subworkflows
│
├── modules/
│   ├── local/                       # Custom modules
│   │   ├── aggregator.nf           # Thin wrapper → bin/staphit-aggregate
│   │   ├── report.nf               # → bin/staphit-report
│   │   ├── outbreak_cluster.nf     # → bin/staphit-cluster
│   │   ├── qc_gate.nf             # → bin/staphit-qc
│   │   ├── visualization.nf       # → bin/staphit-visualize
│   │   └── ...
│   └── nf-core/                    # Standard nf-core modules (fastqc, quast, etc.)
│
├── bin/                            # All proven staphit-* scripts
├── species/                        # Species descriptors (future multi-species)
│   └── s_aureus.yml
├── conf/                           # nf-core configs
├── assets/                         # Schemas, templates, test data
└── tests/
    ├── python/                     # pytest for bin/ scripts
    └── nf-test/                    # Process-level tests
```

The species descriptor (`species/s_aureus.yml`) declares:
- Genome size, GC range, assembly thresholds
- Which typing tools to activate
- Virulence gene sets and interpretation rules
- Point mutation → drug class mappings
- Outbreak SNP/cgMLST thresholds
- Reference genome
- Vitek card definitions

For v2, only *S. aureus* is implemented. But the structure ensures adding *E. coli* or *K. pneumoniae* later is **additive** — new files in `subworkflows/local/species/` and `species/`, zero changes to the core.

---

## The Plan — Ordered by User Experience

### Phase 1: Foundation

**What the user gets:** `nf-core pipelines launch alarawms/staphit2` shows a working param form.

#### Step 1.1: nf-core scaffold + species descriptor

- The nf-core template is already created
- Create `species/s_aureus.yml` with all organism-specific parameters
- Create `subworkflows/local/core/` and `subworkflows/local/species/s_aureus/` directory structure
- Define `nextflow_schema.json` with all params from v1 organized into groups
- Create `conf/test.config` with paths to bundled test data

#### Step 1.2: Bundled test dataset

- Create 5 synthetic *S. aureus* read sets (subsample from real data or simulate)
- Include: samplesheet, sample_metadata.csv, antibiogram.csv
- Store in `assets/test-data/` (small enough for GitHub) or nf-core/test-datasets
- `-profile test,docker` should complete in <10 minutes

#### Step 1.3: Copy bin/ scripts from v1

- All 12 `bin/staphit-*` scripts copied verbatim (they're standalone, they work)
- Verify: each script's `--help` works
- Copy all 216 pytest tests into `tests/python/`
- Verify: `pytest tests/python/ -v` passes

---

### Phase 2: Data Acquisition

**What the user gets:** "Find all MRSA from Saudi Arabia on SRA, download them, get a ready-to-run samplesheet."

#### Step 2.1: Create bin/staphit-fetch

Search SRA/ENA and bridge to nf-core/fetchngs:

```bash
# Search: what's available?
staphit-fetch search \
    --organism "Staphylococcus aureus" \
    --country "Saudi Arabia" \
    --date-range 2020-2024 \
    --limit 500

# Output: search_results.tsv with accessions + available metadata + completeness scores

# Preview: metadata quality
staphit-fetch preview search_results.tsv
# "Found 342 samples: 310 have dates, 280 have country, 45 have host_disease"

# Fetch: download + map metadata
staphit-fetch download search_results.tsv \
    --output-dir fetched/ \
    --max-samples 200

# Uses nf-core/fetchngs under the hood or direct fasterq-dump
# Outputs: samplesheet.csv + sample_metadata.csv (mapped to PHA4GE schema)
# Reports: metadata_completeness.tsv showing what's filled vs missing
```

Functions:
- `search()` — queries NCBI Entrez with organism + filters, returns accession list with BioSample metadata
- `preview()` — shows metadata completeness stats
- `download()` — fetches reads (fasterq-dump or fetchngs), maps BioSample attributes to PHA4GE fields
- `merge_local()` — combine fetched samplesheet with local samples

Tests: mock NCBI responses, verify metadata mapping.

#### Step 2.2: Wire SRA search into the pipeline (optional)

Keep `--species` param for direct search-and-run (v1 behavior), but recommend the 2-step staphit-fetch workflow for better control.

---

### Phase 3: Input Preparation

**What the user gets:** "Convert my Vitek PDFs, lab spreadsheets, and clinical data into pipeline-ready CSVs."

#### Step 3.1: Migrate staphit-metadata

Already in `bin/` from Step 1.3. Create nf-core-compatible input validation:

- `assets/schema_input.json` — validates samplesheet format
- `assets/schema_metadata.json` — validates sample_metadata.csv
- `assets/schema_antibiogram.json` — validates antibiogram.csv

The `staphit-metadata validate` command aligns with these schemas.

#### Step 3.2: Add VALIDATE_METADATA as nf-core local module

```
modules/local/validate_metadata.nf
```

Thin wrapper calling `bin/staphit-metadata validate` + `normalize`.

---

### Phase 4: QC & Assembly (Core)

**What the user gets:** "My reads are trimmed, assembled, and quality-checked. Bad samples are flagged."

#### Step 4.1: Install nf-core modules

```bash
nf-core modules install fastqc
nf-core modules install trimgalore  # or trimmomatic
nf-core modules install quast
```

#### Step 4.2: Create local assembly module

```
modules/local/skesa.nf        # SKESA assembly
```

nf-core has `spades` module. SKESA needs a local module.

#### Step 4.3: Create QC gate subworkflow

```
subworkflows/local/core/qc.nf
```

Wires: TRIMGALORE → FASTQC → ASSEMBLY → QUAST → CHECKM2 → QC_GATE

Assembly size filter (>500KB) built into the subworkflow. CheckM2 runs with `--skip_qc_gate` option.

#### Step 4.4: nf-test for QC processes

Test: SKESA produces assembly >500KB for test data. CheckM2 reports >90% completeness. QC_GATE correctly filters.

---

### Phase 5: Typing (Species-Specific)

**What the user gets:** "Each sample has MLST, spa type, SCCmec, agr group, and optionally cgMLST."

#### Step 5.1: Install nf-core modules

```bash
nf-core modules install mlst
```

#### Step 5.2: Create S. aureus typing subworkflow

```
subworkflows/local/species/s_aureus/typing.nf
```

Wires: MLST + SPATYPER + SCCMEC + AGR_TYPING (+ optional CHEWBBACA for cgMLST)

Local modules:
```
modules/local/spatyper.nf
modules/local/sccmec.nf
modules/local/agr_typing.nf
modules/local/chewbbaca.nf     # 4 processes: prep, allele, join, dists
```

#### Step 5.3: Species detection via Mash

```
modules/local/mash.nf
```

Auto-confirm species. Future: auto-select species adapter if multiple are available.

#### Step 5.4: nf-test for typing processes

Test: MLST returns valid ST for test data. SPATYPER returns a spa type. SCCmec returns a type.

---

### Phase 6: Resistance (Core + Species-Specific)

**What the user gets:** "Every resistance gene, every point mutation, predicted phenotype, and concordance with lab MIC data."

#### Step 6.1: Install nf-core modules

```bash
nf-core modules install amrfinderplus
nf-core modules install abricate
```

#### Step 6.2: Create AMR subworkflow

```
subworkflows/local/core/amr.nf
```

Wires: AMRFINDERPLUS + ABRICATE (resfinder + vfdb + plasmidfinder) + KMA

#### Step 6.3: Create S. aureus mutation subworkflow

```
subworkflows/local/species/s_aureus/mutations.nf
```

Uses `bin/staphit-mutations` — already handles AMRFinderPlus POINT rows with S. aureus-specific gene→drug mapping.

#### Step 6.4: nf-test

Test: AMRFinderPlus finds mecA in MRSA test sample. Point mutation profiler extracts gyrA_S84L.

---

### Phase 7: Virulence (Species-Specific)

**What the user gets:** "Structured virulence profile: PVL status, TSST, IEC type, enterotoxin count, operon completeness."

#### Step 7.1: Create S. aureus virulence subworkflow

```
subworkflows/local/species/s_aureus/virulence.nf
```

Uses `bin/staphit-virulence` — the profiler is already a standalone script. The Nextflow module is a thin wrapper.

#### Step 7.2: nf-test

Test: PVL+ sample correctly identified. IEC type matches expected.

---

### Phase 8: Plasmids (Core)

**What the user gets:** "Full plasmid reconstruction with mobility classification."

#### Step 8.1: Create plasmid subworkflow

```
subworkflows/local/core/plasmids.nf
```

Wires: MOB_RECON → PLASMID_SUMMARY

Local modules:
```
modules/local/mob_recon.nf
modules/local/plasmid_summary.nf
```

#### Step 8.2: nf-test

Test: MOB_RECON produces contig_report.txt. Handles samples with no plasmids gracefully.

---

### Phase 9: Phylogeny & Outbreak Detection (Core)

**What the user gets:** "Phylogenetic tree, SNP distance matrix, outbreak clusters with epi annotation."

#### Step 9.1: Create phylogeny subworkflow

```
subworkflows/local/core/phylogeny.nf
```

Wires: alignment (PANAROO or SNIPPY_CORE) → tree (IQTREE or FASTTREE) → SNP_DISTS

Alignment method and tree builder are independent choices (the v1 refactor we did).

#### Step 9.2: Create clustering subworkflow

```
subworkflows/local/core/clustering.nf
```

Wires: SNP_DISTS + optional CGMLST_DISTS → OUTBREAK_CLUSTER

Uses `bin/staphit-cluster` — thresholds come from species descriptor.

#### Step 9.3: Install nf-core modules

```bash
nf-core modules install iqtree
nf-core modules install snpdists
```

Local: `modules/local/panaroo.nf`, `modules/local/snippy.nf`, `modules/local/fasttree.nf`

#### Step 9.4: nf-test

Test: Panaroo produces alignment. IQ-TREE produces treefile. Clustering assigns cluster IDs.

---

### Phase 10: Aggregation & Reporting

**What the user gets:** "Per-sample JSON reports, merged summary TSV, Markdown run report, publication figures."

#### Step 10.1: Create reporting subworkflow

```
subworkflows/local/core/reporting.nf
```

Wires: AGGREGATOR → SUMMARY_MERGER → REPORT → VISUALIZATION → MULTIQC

All use `bin/staphit-*` scripts as thin wrappers.

#### Step 10.2: SQLite backend (new in v2)

Create `bin/staphit-db`:
- `staphit-db init` — creates schema
- `staphit-db import results/` — loads pipeline output
- `staphit-db query "SELECT ..."` — SQL access
- `staphit-db export --format tsv` — dump

Tables: samples, typing, resistance, virulence, clusters, antibiogram, qc, plasmids.

The aggregator writes to SQLite in addition to JSON/CSV (backwards compatible).

#### Step 10.3: nf-test

Test: Aggregator produces valid JSON. Report produces valid Markdown.

---

### Phase 11: Dashboard (Separate Package)

**What the user gets:** "Interactive surveillance dashboard with linked views."

#### Step 11.1: Create staphit2-dashboard repo

```bash
gh repo create alarawms/staphit2-dashboard --private
```

- Move `dashboard/` package from v1
- Add `pyproject.toml` for pip install
- Reads from SQLite (primary) or flat files (fallback)
- `pip install staphit2-dashboard && staphit2-dashboard results/`

#### Step 11.2: Apply Plotly template to all charts

Use `utils.PLOTLY_TEMPLATE` consistently across all analysis modules.

#### Step 11.3: Full UI polish

Use frontend design best practices for the Clinical Intelligence theme.

---

### Phase 12: Continuous Surveillance

**What the user gets:** "Point at an incoming directory, get automatic analysis and alerts."

#### Step 12.1: Migrate staphit-watch

Already in `bin/` from Step 1.3. Add integration with SQLite:
- New clusters detected by querying DB instead of comparing flat files
- Notification includes cluster details from SQL query

#### Step 12.2: Documentation

- `docs/usage.md` — complete usage guide following user journey
- `docs/output.md` — every output file documented
- `docs/surveillance.md` — watch mode, incremental runs, alerts

---

### Phase 13: Testing & CI

**What the user gets:** confidence that the pipeline works.

#### Step 13.1: nf-test for all local modules

Priority order (most bug-prone first from v1 experience):
1. checkm2 (DB path, symlinks, name collisions)
2. mob_recon (empty output, name collisions)
3. iqtree (seed tree flag)
4. aggregator (quoting — now eliminated)
5. All others

#### Step 13.2: GitHub Actions

Enable from nf-core template:
- `ci.yml` — lint + pytest + `nextflow run -profile test`
- `linting.yml` — nf-core lint
- Push to dev → CI → green badge

#### Step 13.3: Conda environment

`environment.yml` per module for users without Docker.

---

### Phase 14: Release

#### Step 14.1: Versioned release

- CHANGELOG.md
- `v2.0.0` tag
- GitHub Release with notes
- Make repo public (when ready)

#### Step 14.2: nf-core submission (optional)

- Submit to nf-core for community review
- Requires passing nf-core lint, tests, documentation
- Gets listed on nf-co.re/pipelines

---

## Execution order

| Phase | What | Depends on | Effort |
|-------|------|-----------|--------|
| 1 | Foundation (scaffold, test data, copy scripts) | Nothing | 3 days |
| 2 | Data acquisition (staphit-fetch) | Phase 1 | 3 days |
| 3 | Input preparation (metadata tools) | Phase 1 | 1 day |
| 4 | QC & Assembly | Phase 1 | 2 days |
| 5 | Typing | Phase 4 | 2 days |
| 6 | Resistance | Phase 4 | 2 days |
| 7 | Virulence | Phase 6 | 1 day |
| 8 | Plasmids | Phase 4 | 1 day |
| 9 | Phylogeny & clustering | Phase 4 | 2 days |
| 10 | Reporting & SQLite | Phase 5-9 | 3 days |
| 11 | Dashboard | Phase 10 | 3 days |
| 12 | Surveillance | Phase 10 | 1 day |
| 13 | Testing & CI | Phase 4+ | 3 days (parallel) |
| 14 | Release | All | 1 day |

**Critical path:** Phase 1 → 4 → 5/6/7/8 (parallel) → 9 → 10 → 14

**Total:** ~4 weeks focused work, or 6-8 weeks part-time.

---

## Multi-species scaffold summary

When ready to add *E. coli*:

1. Create `species/e_coli.yml` — genome size, thresholds, typing tools
2. Create `subworkflows/local/species/e_coli/typing.nf` — serotype, phylogroup
3. Create `bin/staphit-virulence-ecoli` (or make `staphit-virulence` species-aware) — stx, pathotype
4. Create `bin/staphit-mutations-ecoli` — E. coli-specific gyrA positions
5. Add Vitek card `AST-N395.yml` to `bin/staphit-metadata`
6. Zero changes to core subworkflows

The species descriptor pattern means the **framework is already built** — it's just serving *S. aureus* first.
