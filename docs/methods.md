# alarawms/staphit2: Methods and validation

This page explains how staphit2 decides what to keep, what to flag and how to group samples, and why each default is what it is. Every number quoted under "Validation" comes from one run of the pipeline on 703 public and Saudi _S. aureus_ samples (684 passed to typing: 639 Illumina, 45 ONT-only), unless stated otherwise. Parameters are described in [usage](usage.md); files in [output](output.md).

- [Sample fate: dropped or flagged](#sample-fate-dropped-or-flagged)
- [Assembly](#assembly)
- [Species and contamination](#species-and-contamination)
- [Mixed _S. aureus_ cultures](#mixed-s-aureus-cultures)
- [Long-read polishing and gene integrity](#long-read-polishing-and-gene-integrity)
- [Resistance genes: consensus across three callers](#resistance-genes-consensus-across-three-callers)
- [cgMLST](#cgmlst)
- [Outbreak clustering](#outbreak-clustering)
- [Recombination (Gubbins)](#recombination-gubbins)
- [Time-scaled phylogeny (BEAST2)](#time-scaled-phylogeny-beast2)
- [Reproducibility](#reproducibility)

## Sample fate: dropped or flagged

A sample is **dropped** only when its data cannot support typing: no assembly, an assembly of the wrong size, failed completeness/contamination, or not _S. aureus_. Everything else that weakens a result is a **flag**: the sample stays in every output and the reason is written next to it.

`report/sample_status.tsv` lists every input sample with its status (`included`, `not_in_tree`, `dropped`), the stage at which it was dropped and why, and its flags. The workflow records which stages each sample reached; the first missing stage is the reason it was dropped. The same flags appear in the `qc_flags` column of `summary/combined_summary.tsv`.

| Stage / flag                         | Drops or flags | Rule (default)                                                                                  |
| ------------------------------------ | -------------- | ----------------------------------------------------------------------------------------------- |
| `assembly_failed`                    | drop           | assembler produced no contigs                                                                   |
| `assembly_size_out_of_range`         | drop           | assembly length outside 2.5–3.2 Mb                                                              |
| `qc_gate_failed`                     | drop           | CheckM2 completeness < 90% or contamination > 5%                                                |
| `not_s_aureus`                       | drop           | ANI to NCTC 8325 < 95%                                                                          |
| `typing_incomplete`                  | drop           | a typing or AMR step failed for the sample                                                      |
| `not_in_tree`                        | kept           | in the summary but not in the phylogeny                                                         |
| `low_depth_short` / `low_depth_long` | flag           | trimmed Illumina < 30x / ONT < 20x                                                              |
| `contigs`, `n50`                     | flag           | > 500 contigs / N50 < 10 kb                                                                     |
| `mixed_strains`                      | flag           | > 800 mixed sites ([below](#mixed-s-aureus-cultures))                                           |
| `cgmlst_missing`                     | flag           | > 86 of 1,716 cgMLST loci without a call                                                        |
| `truncated_genes`                    | flag           | long-only assembly with > 80 truncated genes ([below](#long-read-polishing-and-gene-integrity)) |

**Validation.** All 19 samples dropped from the 703-sample run are explained: 17 by the QC gate, 1 not _S. aureus_ (_Enterobacter cloacae_) and 1 empty assembly. The size gate replaces an earlier 500 kb file-size filter and drops the same 19; the depth/contig/N50 flags mark 12 kept samples, including 4 ONT samples under 20x.

## Assembly

Read mode is detected per sample from the samplesheet: Illumina only (SKESA, or SPAdes with `--assembler spades`), ONT only (Flye + Racon + Medaka via Dragonflye) or hybrid (Flye + Racon + Polypolish with the Illumina reads). Illumina reads are trimmed with fastp and subsampled to `--target_depth` (100x) with Rasusa before assembly.

_S. aureus_ genomes are 2.7–2.9 Mb, so the gate drops assemblies outside 2.5–3.2 Mb (sequence length, not file size). Too small means a failed or partial assembly; too large usually means two genomes.

Known limitation: SKESA breaks assemblies at repeats, and 246 of 684 SCCmec cassettes were split across contigs. The SCCmec typer falls back to the reads for split cassettes; a comparison with SPAdes/shovill is open (review item B7).

## Species and contamination

fastANI against NCTC 8325 confirms the species (≥ 95% ANI). Mash screen (`-w`, winner-takes-all) runs on every assembly against RefSeq and reports the best and second-best bacterial species; phage, prophage and plasmid hits are ignored, and unnamed `sp.` hits rank below named species. A second species at ≥ 0.90 identity with ≥ 100/1000 shared hashes is reported as `possible_contamination` in `species_qc/species_confirmed.tsv`.

**Validation.** Flags SAMEA114040995 (_Lysinibacillus fusiformis_ alongside _S. aureus_). Mash only sees other species; a second _S. aureus_ strain needs the next check.

## Mixed _S. aureus_ cultures

**Method.** The trimmed Illumina reads of each sample are mapped back to the sample's own assembly (minimap2 `-ax sr`), and bcftools counts allele depths at every position (mapping quality ≥ 30, base quality ≥ 20). A position is a _mixed site_ when it has ≥ 20 reads and the second most common base has ≥ 5 reads and ≥ 10% of them. A single strain agrees with its own assembly almost everywhere; two _S. aureus_ strains disagree at every position where they differ, which is thousands for different lineages. The median minor-base fraction across mixed sites estimates the minor strain's share. This is the idea behind ConFindr, applied genome-wide and without the rMLST database. Long-only samples are not checked: ONT errors would swamp the signal.

**Calibration.** On the 655 Illumina samples of the run:

| Mixed sites | Samples |
| ----------- | ------- |
| < 200       | 456     |
| 200–399     | 169     |
| 400–599     | 5       |
| 600–919     | 0       |
| ≥ 920       | 25      |

The background in single-strain samples (median 157) comes from repeats and sequencing errors. The 25 samples above the gap have a minor strain at 11–35%; 6 of them had already failed the CheckM2 QC gate. Of the 19 that passed it, the 11 above 4,000 mixed sites also miss 47–531 cgMLST loci (typical sample: 34), an independent sign of a mixed assembly; the 8 between 920 and 4,000 mostly look normal to cgMLST (32–38 missing, two at 48 and 124). A low-level mixture leaves the assembly almost intact, which is why it needs this check. `--max_mixed_sites` 800 sits in the empty interval. The most extreme, SAMEA114860523, has 20,920 mixed sites with a 15% minor strain; it passed every other QC check (1,192x depth, 579 contigs).

**Interpretation.** A flagged sample's typing describes the majority strain, but its SNP and cgMLST distances to other samples are unreliable, and its AMR genes may come from either strain. Re-isolate if it matters.

## Long-read polishing and gene integrity

**Method.** Long-only assemblies are polished with Medaka by default (`--medaka_model auto`):

1. If the read headers name the basecaller model (`basecall_model_version_id=`, written by MinKNOW/Dorado), Medaka uses the matching model, downloading it if it is not bundled.
2. If they name a model that cannot be downloaded, the sample is left unpolished: the wrong model can add errors.
3. If they name none, Medaka's bundled R10.4.1 SUP model is used. Public reads usually land here: ENA rewrites FASTQ headers.
4. Samples with a NanoPlot median read quality below `--medaka_min_read_q` (10) are not polished.

Hybrid assemblies are polished with the Illumina reads (Polypolish) instead. The model used is the last line of `dragonflye/<sample>/<sample>.dragonflye.log`.

Uncorrected homopolymer errors shift the reading frame and truncate genes. `TRUNCATED_GENES` runs blastp of the Prokka proteins against the NCTC 8325 proteins (from the GenBank reference) and counts proteins whose best hit is the same gene (≥ 95% identity) but which are shorter than 90% of it.

**Validation.** Illumina assemblies have a median of 46 truncated genes (60 samples; real pseudogenes and contig ends). The 45 long-only assemblies (median read Q 7–22, R10.4.1-like) assembled with Flye + Racon alone had 33–75, except one Q8 sample with 97–108. Medaka on five of them:

| Sample         | Median read Q | Truncated, Racon only | After Medaka |
| -------------- | ------------- | --------------------- | ------------ |
| SAMEA117686580 | 11.0          | 44                    | 44           |
| SAMEA117686603 | 10.0          | 44                    | 44           |
| SAMEA117686582 | 20.4          | 65                    | 52           |
| SAMEA117686615 | 21.7          | 74                    | 55           |
| SAMEA117686597 | 8.0           | 108                   | 136          |

Medaka helps where there is something to fix and does nothing where there is not, but on very poor reads the default model made the assembly worse, hence the read-quality gate. `--max_truncated_genes` 80 is above every long-only assembly except the Q8 one.

## Resistance genes: consensus across three callers

**Method.** AMRFinderPlus and ResFinder (via ABRicate) search the assembly; KMA maps the reads to ResFinder. Gene names are reduced to families so the tools can be compared (allele numbers and accessions removed; known synonyms such as `blaPC1` = `blaZ` merged).

- `amr_consensus`: every AMRFinderPlus gene, plus genes ResFinder and KMA both find. AMRFinderPlus is the primary caller because its _S. aureus_ curation covers families ResFinder misses: `fosB` was called in 428 samples by AMRFinderPlus and in none by ResFinder, so a "two tools agree" rule would drop it.
- `amr_discordant`: genes only one secondary tool found, tagged `[resfinder]`, `[kma]`, `[kma-trace]` or `[kma-ont]`.

KMA depth is compared with the sample's read depth (bases after trimming / genome size, capped at the subsampling target). A gene at < 25% of the read depth cannot be in every cell: it is tagged `kma-trace` (contamination or index-hopping) and never counts. KMA on raw ONT reads of long-only samples is tagged `kma-ont` and never confirms a gene on its own.

**Validation.** 134 KMA calls were traces at 1–13% of the read depth, including _mecA_ in 6 SCCmec-negative samples at ~10%: counted naively, they would have been reported as MRSA. 58 `kma-ont` calls, mostly _erm(C)_ at 3,500–8,900x, look like small high-copy plasmids that Flye did not assemble.

## cgMLST

**Method.** The 1,716 loci of PubMLST _S. aureus_ cgMLST scheme 20 are downloaded once (`bin/staphit-cgmlst-schema`, ~840 MB) and prepared for chewBBACA (`PrepExternalSchema`). All samples of a run are called in a single `AlleleCall` (`--mode 4`, `--no-inferred`): new alleles are numbered consistently only within one run, so calling samples one by one would give two different new alleles the same number. The allele distance between two samples is the number of loci called in both whose alleles differ (pairwise deletion, as in SeqSphere+ and PubMLST); calls such as LNF, PLOT, ASM/ALM and NIPH count as missing.

**Validation.** Allele calling for 684 samples took 7 minutes on 32 CPUs. A typical sample leaves 34 loci without a call (95th percentile 44); samples with hundreds missing are fragmented, contaminated or mixed, hence the `cgmlst_missing` flag above 86 (5%).

## Outbreak clustering

**Method.** Samples are linked when their distance is at or below a threshold, and clusters are the connected groups (single linkage), at three tiers: `direct_transmission`, `outbreak`, `related`. By default the distance is cgMLST alleles (`--cluster_by cgmlst`). Unlike a core-genome alignment, a cgMLST distance between two samples does not change when unrelated lineages are added to the run: the Panaroo core alignment shrinks as the collection diversifies, so SNP thresholds on it mean different things in different runs. SNP distances are still reported (`max_snp_distance`, `transmission_pairs.tsv`) and can be used with `--cluster_by snp` or `both`. A cluster is `chain_linked` when its diameter exceeds twice the tier threshold: single linkage joined distant samples through intermediates.

**Thresholds.** The defaults and their sources are in [usage: outbreak clustering thresholds](usage.md#outbreak-clustering-thresholds). In short: cgMLST `10,24,50` alleles (24 is the SeqSphere+ _S. aureus_ cluster threshold; PubMLST's LIN-code bins for scheme 20 are at 27, 14, 7, 4, 2 and 1); SNP `15,50,200` on the Panaroo core alignment (15 after Coll et al. 2020, who set it on whole-genome SNPs).

**Validation.** For 13,527 sample pairs within 200 alleles, whole-genome SNPs (Snippy alignment against NCTC 8325, compared pair by pair) were set against cgMLST distances. Pairs within 15 SNPs differed by at most 9 alleles (95th percentile 8), and one allele corresponded to a median of 2.8 SNPs, so the 10-allele tier sits where 15–25 SNPs would. This shows the tiers are consistent with SNP practice; it does not calibrate them. Only epidemiologically confirmed outbreaks can: give `bin/staphit-calibrate-thresholds` the distance matrix and a `sample_id<TAB>outbreak` list (`sporadic` for unlinked samples) and it reports, for every threshold, the share of linked pairs caught (sensitivity), the share of unlinked pairs kept apart (specificity) and the threshold with the highest Youden's J.

## Recombination (Gubbins)

Gubbins (3.4) runs by default on Snippy's whole-genome alignment (`--skip_gubbins` turns it off). Imported segments carry many SNPs at once, so they lengthen branches and inflate SNP distances without reflecting time since a common ancestor. Gubbins reports the recombinant regions and r/m per branch, builds a recombination-aware tree, and masks those regions; `snpdists/gubbins_masked.tsv` holds SNP distances after masking. Clustering is not switched to the masked distances: the default cgMLST clustering is already robust to recombination (an imported segment changes at most one allele per locus), and the SNP thresholds were set on unmasked alignments.

## Time-scaled phylogeny (BEAST2)

With `--use_beast`, BEAST2 (2.7) dates the Gubbins-masked alignment, since recombination would inflate the clock rate. It gets the variable sites of the masked alignment plus the counts of invariant A/C/G/T sites (`constantSiteWeights`), so rates are not overestimated from analysing variable sites only. Model: HKY + Γ4, strict clock, constant-size coalescent. With a collection date for every sample, the clock rate is estimated under a log-normal prior centred on ~1 × 10⁻⁶ substitutions/site/year; otherwise the clock is fixed and branch lengths are substitutions/site. Check convergence (ESS > 200 in Tracer on `beast2/beast_run.log`) before using dates; the default chain is 10 million steps.

Before this was fixed, `--use_beast` could not run: the generated XML used BEAST 2.6 class names that 2.7 rejects, the container could not load its packages under Docker's `-u`, and TreeAnnotator's arguments had changed.

## Reproducibility

Every tool and database version (AMRFinderPlus database date, ResFinder commit, Mash sketch checksum, container digests) is written to `pipeline_info/staphit2_software_mqc_versions.yml` and the MultiQC report. All containers are pinned. The real-data test (`-profile test`: two Illumina, one hybrid and one ONT-only sample from [staphit2-test-data](https://github.com/alarawms/staphit2-test-data)) checks the expected STs and SCCmec types and snapshots the versions file, so any version change shows up in CI.
