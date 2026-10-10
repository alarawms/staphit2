# staphit2 review — to-do

Review of staphit2 against nf-core standards (A) and scientifically (B), October 2026.
Tick items when done and verified; link the PR.

A1–A11: branch `chore/nf-core-standards`. Notes:

- A1: local modules emit `[process, tool, version]` on the `versions` topic (same style as the nf-core modules); FastQC, IQ-TREE and snp-dists `versions.yml` are now collected too.
- A3: databases are published to `<outdir>/databases/` and reused with `--checkm2_db`, `--amrfinder_db`, `--mash_db`; ResFinder pinned by `--resfinder_db_commit`; CheckM2/Mash checksums recorded. (`storeDir` cannot be used with tuple/`eval` outputs.)
- A2 found two broken containers: FastTree was missing from `staphb/lyveset:2.0.1`, and the BEAST2 tag did not exist.

## A. nf-core standards

- [x] **A1** `versions.yml` from every local module (MultiQC software-versions section complete)
- [x] **A2** Pin all containers (no `:latest`); record image tags
- [x] **A3** Databases: pinned versions, reusable via params, versions written to outputs
      (CheckM2 DB, AMRFinderPlus DB, ResFinder DB, Mash RefSeq sketch)
- [x] **A4** `process.resourceLimits`; resource labels aligned with nf-core defaults; heavy tools via `withName`
- [x] **A5** `publishDir` moved from modules to `conf/modules.config`
- [x] **A6** `errorStrategy` moved from modules to `conf/modules.config`
- [x] **A7** nf-test: tests for the new modules (fastp, NanoPlot, Dragonflye, ITOL_EXPORT) + pipeline stub test with short / hybrid / long-only samples
- [x] **A8** `stub:` blocks for all local modules
- [x] **A9** Pinned Python image (not EOL `python:3.9`); call `bin/` scripts by name, not `${projectDir}/bin/...`
- [x] **A10** Docs: `usage.md` (long_fastq, accessions/metadata helpers, new params), `output.md` (fastp, NanoPlot, Dragonflye, iTOL, SCCmec IWG columns), `CITATIONS.md` (all tools)
- [x] **A11** Conda/Singularity: document docker/singularity-only support for local modules
- [x] **A15** `PLOT_TREE` runs on `ghcr.io/alarawms/staphit2-r:1.0.0` (pinned R 4.5.3 + ggtree stack); no run-time package installs
- [x] **A17** `sccmec_container` default pinned to the sccmec_typer build from merged PR #4 (`sha256:f4ba10ae…`, has `--fallback-1/2`)
- [x] **A16** CI uses nf-test 0.9.5 (strict-syntax support); the `NXF_SYNTAX_PARSER=v1` workaround is removed
- [ ] **A18** Real minimal test dataset for `-profile test` (the old `assets/test-data` FASTQs were never committed); CI currently runs the stub pipeline test (`tests/default.nf.test`)
- [x] **A12** Remove dead code (unused Trim Galore module, unused `reporting.nf`)
- [x] **A13** Deprecated syntax (`Channel.` → `channel.` etc.) in local workflows/subworkflows
- [x] **A14** Lint: `.gitignore` template difference — accept and ignore in `.nf-core.yml`

## B. Scientific

- [x] **B1** Reproducibility: every tool and database version is in `pipeline_info/staphit2_software_mqc_versions.yml` and MultiQC (A1–A3); all containers pinned; `tests/test_profile.nf.test` snapshots the versions so any change shows up in CI
- [ ] **B2** Within-species contamination check (mixed _S. aureus_ cultures; e.g. ConFindr or allele-ratio from read mapping) — e.g. SAMEA114860523: 1,192x, 579 contigs, passed QC
- [ ] **B3** Lineage-aware outbreak clustering (cgMLST on by default, or within-lineage mapping + recombination masking) instead of one cross-lineage core alignment / NCTC8325 mapping
- [x] **B4** Assembly gate: drop outside 2.5–3.2 Mb (`--min/max_assembly_size`, replaces the 500 kb filter); flag (keep) low depth (<30x Illumina, <20x ONT), >500 contigs, N50 <10 kb: `qc_flags` in the summary, `flags` in `sample_status.tsv`
- [x] **B5** `report/sample_status.tsv` + _Sample Status_ report section: every input sample with the stage and reason it was dropped (on the 703-sample run: 17 QC gate, 1 non-_S. aureus_, 1 empty assembly)
- [ ] **B6** ONT-only polishing: Medaka by default with model from the FASTQ header; flag truncated genes
- [ ] **B7** Assembler choice: SKESA fragments repeat regions (246/684 split SCCmec cassettes); evaluate SPAdes/shovill
- [x] **B8** `--use_beast`: Gubbins masks recombination on snippy-core's whole-genome alignment; BEAST2 gets the masked SNPs plus invariant-site counts (`constantSiteWeights`); XML rewritten for BEAST 2.7 (the old one never parsed); prep runs in the Python image; BEAST runs under Docker `-u` (user.home)
- [ ] **B9** Consensus _S. aureus_ resistance calls across AMRFinderPlus / ResFinder / KMA (down-weight KMA on raw ONT)
- [x] **B10** Species QC reports best and second-best bacterial species (Mash screen -w on every assembly, phage/plasmid hits ignored) and flags `possible_contamination` (≥0.90 identity, ≥100/1000 hashes); within-species mixtures remain for B2
- [ ] **B11** Cluster thresholds: document scheme per threshold (PubMLST scheme 20, 1,716 loci); recalibrate on known outbreaks
