# staphit2 review — to-do

Review of staphit2 against nf-core standards (A) and scientifically (B), October 2026.
Tick items when done and verified; link the PR.

## A. nf-core standards

- [ ] **A1** `versions.yml` from every local module (MultiQC software-versions section complete)
- [ ] **A2** Pin all containers (no `:latest`); record image tags
- [ ] **A3** Databases: pinned versions, cached (`storeDir`), versions written to outputs
      (CheckM2 DB, AMRFinderPlus DB, ResFinder DB, Mash RefSeq sketch)
- [ ] **A4** `process.resourceLimits`; resource labels aligned with nf-core defaults; heavy tools via `withName`
- [ ] **A5** `publishDir` moved from modules to `conf/modules.config`
- [ ] **A6** `errorStrategy` moved from modules to `conf/modules.config`
- [ ] **A7** nf-test: tests for the new modules (fastp, NanoPlot, Dragonflye, ITOL_EXPORT) + pipeline stub test with short / hybrid / long-only samples
- [ ] **A8** `stub:` blocks for all local modules
- [ ] **A9** Pinned Python image (not EOL `python:3.9`); call `bin/` scripts by name, not `${projectDir}/bin/...`
- [ ] **A10** Docs: `usage.md` (long_fastq, accessions/metadata helpers, new params), `output.md` (fastp, NanoPlot, Dragonflye, iTOL, SCCmec IWG columns), `CITATIONS.md` (all tools)
- [ ] **A11** Conda/Singularity: document docker/singularity-only support for local modules
- [x] **A12** Remove dead code (unused Trim Galore module, unused `reporting.nf`)
- [x] **A13** Deprecated syntax (`Channel.` → `channel.` etc.) in local workflows/subworkflows
- [x] **A14** Lint: `.gitignore` template difference — accept and ignore in `.nf-core.yml`

## B. Scientific

- [ ] **B1** Reproducibility over time: tool + database versions recorded with every result (depends on A1–A3)
- [ ] **B2** Within-species contamination check (mixed *S. aureus* cultures; e.g. ConFindr or allele-ratio from read mapping) — e.g. SAMEA114860523: 1,192x, 579 contigs, passed QC
- [ ] **B3** Lineage-aware outbreak clustering (cgMLST on by default, or within-lineage mapping + recombination masking) instead of one cross-lineage core alignment / NCTC8325 mapping
- [ ] **B4** Stricter assembly gate: genome size ~2.5–3.2 Mb, minimum depth (short ~30x, ONT ~20–30x), contig/N50 flags
- [ ] **B5** Report every dropped sample (ignored assembler/Snippy failures) with a reason in the QC report
- [ ] **B6** ONT-only polishing: Medaka by default with model from the FASTQ header; flag truncated genes
- [ ] **B7** Assembler choice: SKESA fragments repeat regions (246/684 split SCCmec cassettes); evaluate SPAdes/shovill
- [ ] **B8** Recombination masking (Gubbins/ClonalFrameML) before BEAST2 dating
- [ ] **B9** Consensus *S. aureus* resistance calls across AMRFinderPlus / ResFinder / KMA (down-weight KMA on raw ONT)
- [ ] **B10** Species QC: report second-best hit alongside NCTC8325 ANI
- [ ] **B11** Cluster thresholds: document scheme per threshold (PubMLST scheme 20, 1,716 loci); recalibrate on known outbreaks
