#!/usr/bin/env python3
"""
beast2_prep.py  —  Generate a BEAST2 2.7 XML from a core-SNP FASTA alignment.

Usage:
    python beast2_prep.py \
        --alignment core.aln \
        --metadata metadata.tsv \
        --output beast_run.xml \
        --chain-length 10000000 \
        --log-every 1000
"""

import argparse
import os
import sys
from datetime import date, datetime


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Prepare a BEAST2 2.7 XML from a FASTA alignment + metadata TSV."
    )
    p.add_argument("--alignment",    required=True,  help="Core SNP FASTA alignment (core.aln)")
    p.add_argument("--metadata",     default=None,   help="Sample metadata TSV (collection_date column)")
    p.add_argument("--output",       required=True,  help="Output XML path")
    p.add_argument("--chain-length", type=int, default=10_000_000, help="MCMC chain length")
    p.add_argument("--log-every",    type=int, default=1_000,      help="Log frequency")
    p.add_argument("--constant-sites", default=None,
                   help="File with 'A C G T' counts of invariant sites (from beast_alignment.py); "
                        "corrects rates when the alignment holds variable sites only")
    return p.parse_args()


# ---------------------------------------------------------------------------
# FASTA parser
# ---------------------------------------------------------------------------

def read_fasta(path):
    """Return list of (name, seq) tuples."""
    sequences = []
    name = None
    seqbuf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, "".join(seqbuf)))
                name = line[1:].split()[0]
                seqbuf = []
            else:
                seqbuf.append(line)
    if name is not None:
        sequences.append((name, "".join(seqbuf)))
    if not sequences:
        sys.exit(f"ERROR: No sequences found in {path}")
    return sequences


# ---------------------------------------------------------------------------
# Date conversion helpers
# ---------------------------------------------------------------------------

def decimal_year(date_str):
    """Convert 'YYYY-MM-DD' or 'YYYY' to decimal year (float)."""
    date_str = date_str.strip()
    if len(date_str) == 4 and date_str.isdigit():
        return float(date_str) + 0.5
    try:
        d = datetime.strptime(date_str, "%Y-%m-%d").date()
        year_start = date(d.year, 1, 1)
        year_end   = date(d.year + 1, 1, 1)
        fraction   = (d - year_start).days / (year_end - year_start).days
        return d.year + fraction
    except ValueError:
        # fallback: try just the year prefix
        year = int(date_str[:4])
        return float(year) + 0.5


# ---------------------------------------------------------------------------
# Metadata parser
# ---------------------------------------------------------------------------

def load_dates(metadata_path, sample_names):
    """
    Read metadata TSV and return a dict {sample_name: decimal_year}.
    Returns empty dict if metadata absent, not parseable, or column missing.
    """
    if metadata_path is None:
        return {}
    sentinel_name = os.path.basename(metadata_path)
    if sentinel_name == "NO_METADATA" or not os.path.isfile(metadata_path):
        return {}

    try:
        with open(metadata_path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
    except Exception as e:
        print(f"WARNING: Could not open metadata file: {e}", file=sys.stderr)
        return {}

    # Support both tab and comma delimiters
    if len(header) == 1:
        with open(metadata_path) as fh:
            header = fh.readline().rstrip("\n").split(",")
        delimiter = ","
    else:
        delimiter = "\t"

    # Find required columns (case-insensitive)
    col_lower = [c.lower().strip() for c in header]
    try:
        sample_col = col_lower.index("sample")
    except ValueError:
        # try common alternatives
        for alt in ("sample_id", "isolate", "strain", "id"):
            if alt in col_lower:
                sample_col = col_lower.index(alt)
                break
        else:
            print("WARNING: No 'sample' column in metadata; running without dates.", file=sys.stderr)
            return {}

    try:
        date_col = col_lower.index("collection_date")
    except ValueError:
        for alt in ("date", "collection date", "isolation_date", "isolation date"):
            if alt in col_lower:
                date_col = col_lower.index(alt)
                break
        else:
            print("WARNING: No 'collection_date' column in metadata; running without dates.", file=sys.stderr)
            return {}

    dates = {}
    sample_set = set(sample_names)
    with open(metadata_path) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split(delimiter)
            if len(parts) <= max(sample_col, date_col):
                continue
            sname = parts[sample_col].strip()
            dval  = parts[date_col].strip()
            if sname in sample_set and dval:
                try:
                    dates[sname] = decimal_year(dval)
                except Exception:
                    pass  # silently skip unparseable dates

    if not dates:
        print("WARNING: Could not parse any dates from metadata; running without date-calibration.", file=sys.stderr)
    return dates


# ---------------------------------------------------------------------------
# BEAST 2.7 XML (BEAST.base class names)
# ---------------------------------------------------------------------------

NAMESPACE = ("beast.base.core:beast.base.inference:beast.base.evolution.alignment:"
             "beast.base.evolution.tree.coalescent:beast.base.util:beast.base.math:"
             "beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:"
             "beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:"
             "beast.base.evolution.likelihood")

# S. aureus: ~1e-6 substitutions/site/year (2-4 SNPs/genome/year); lognormal prior in log space
CLOCK_PRIOR_M, CLOCK_PRIOR_S = -13.8, 1.0

CLOCK_PRIOR = """                <prior id="ClockPrior" spec="beast.base.inference.distribution.Prior" name="distribution" x="@clockRate">
                    <distr spec="beast.base.inference.distribution.LogNormalDistributionModel" M="{m}" S="{s}"/>
                </prior>
"""

CLOCK_OPERATORS = """        <operator id="clockRateScaler" spec="kernel.BactrianScaleOperator" parameter="@clockRate" scaleFactor="0.5" weight="3.0"/>
        <operator id="clockUpDown" spec="beast.base.inference.operator.kernel.BactrianUpDownOperator" scaleFactor="0.75" weight="3.0">
            <up idref="clockRate"/>
            <down idref="Tree"/>
        </operator>
"""

TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- {run_id}: generated by staphit2 bin/beast2_prep.py for BEAST 2.7 -->
<beast version="2.7" namespace="{namespace}" required="">

    <data id="alignment" spec="Alignment" dataType="nucleotide">
{seqs}
    </data>
{filtered}
    <run id="mcmc" spec="MCMC" chainLength="{chain_length}">
        <state id="state" spec="State" storeEvery="5000">
            <tree id="Tree" spec="beast.base.evolution.tree.Tree" name="stateNode">
{taxa}
            </tree>
            <parameter id="kappa" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
            <parameter id="gammaShape" spec="parameter.RealParameter" lower="0.1" name="stateNode">1.0</parameter>
            <parameter id="freqParameter" spec="parameter.RealParameter" dimension="4" lower="0.0" upper="1.0" name="stateNode">0.25</parameter>
            <parameter id="clockRate" spec="parameter.RealParameter" lower="0.0" estimate="{clock_est}" name="stateNode">{clock_value}</parameter>
            <parameter id="popSize" spec="parameter.RealParameter" lower="0.0" name="stateNode">10.0</parameter>
        </state>

        <init id="RandomTree" spec="RandomTree" estimate="false" initial="@Tree" taxa="@alignment">
            <populationModel spec="ConstantPopulation">
                <parameter spec="parameter.RealParameter" name="popSize">10.0</parameter>
            </populationModel>
        </init>

        <distribution id="posterior" spec="CompoundDistribution">
            <distribution id="prior" spec="CompoundDistribution">
                <distribution id="Coalescent" spec="Coalescent">
                    <populationModel spec="ConstantPopulation" popSize="@popSize"/>
                    <treeIntervals spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree"/>
                </distribution>
                <prior id="KappaPrior" spec="beast.base.inference.distribution.Prior" name="distribution" x="@kappa">
                    <distr spec="beast.base.inference.distribution.LogNormalDistributionModel" M="1.0" S="1.25"/>
                </prior>
                <prior id="GammaShapePrior" spec="beast.base.inference.distribution.Prior" name="distribution" x="@gammaShape">
                    <distr spec="beast.base.inference.distribution.Exponential" mean="1.0"/>
                </prior>
                <prior id="PopSizePrior" spec="beast.base.inference.distribution.Prior" name="distribution" x="@popSize">
                    <distr spec="beast.base.inference.distribution.OneOnX"/>
                </prior>
{clock_prior}            </distribution>
            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood" spec="ThreadedTreeLikelihood" data="{data}" tree="@Tree">
                    <siteModel id="SiteModel" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape">
                        <parameter spec="parameter.RealParameter" estimate="false" name="mutationRate">1.0</parameter>
                        <parameter spec="parameter.RealParameter" estimate="false" name="proportionInvariant">0.0</parameter>
                        <substModel id="hky" spec="HKY" kappa="@kappa">
                            <frequencies id="freqs" spec="Frequencies" frequencies="@freqParameter"/>
                        </substModel>
                    </siteModel>
                    <branchRateModel id="StrictClock" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate"/>
                </distribution>
            </distribution>
        </distribution>

        <operator id="kappaScaler" spec="kernel.BactrianScaleOperator" parameter="@kappa" scaleFactor="0.1" weight="0.1"/>
        <operator id="gammaShapeScaler" spec="kernel.BactrianScaleOperator" parameter="@gammaShape" scaleFactor="0.5" weight="0.1"/>
        <operator id="freqExchanger" spec="kernel.BactrianDeltaExchangeOperator" delta="0.01" weight="0.1">
            <parameter idref="freqParameter"/>
        </operator>
        <operator id="popSizeScaler" spec="kernel.BactrianScaleOperator" parameter="@popSize" scaleFactor="0.75" weight="3.0"/>
{clock_ops}        <operator id="treeScaler" spec="beast.base.evolution.operator.kernel.BactrianScaleOperator" scaleFactor="0.5" tree="@Tree" weight="3.0"/>
        <operator id="treeRootScaler" spec="beast.base.evolution.operator.kernel.BactrianScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree" weight="3.0"/>
        <operator id="nodeUniform" spec="beast.base.evolution.operator.kernel.BactrianNodeOperator" tree="@Tree" weight="30.0"/>
        <operator id="subtreeSlide" spec="beast.base.evolution.operator.kernel.BactrianSubtreeSlide" tree="@Tree" weight="15.0"/>
        <operator id="narrow" spec="Exchange" tree="@Tree" weight="15.0"/>
        <operator id="wide" spec="Exchange" isNarrow="false" tree="@Tree" weight="3.0"/>
        <operator id="wilsonBalding" spec="WilsonBalding" tree="@Tree" weight="3.0"/>

        <logger id="tracelog" spec="Logger" fileName="{run_id}.log" logEvery="{log_every}" sanitiseHeaders="true">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood"/>
            <log id="treeHeight" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree"/>
            <log idref="kappa"/>
            <log idref="gammaShape"/>
            <log idref="clockRate"/>
            <log idref="popSize"/>
        </logger>
        <logger id="treelog" spec="Logger" fileName="{run_id}.trees" logEvery="{log_every}" mode="tree">
            <log id="treeWithMetaData" spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree"/>
        </logger>
        <logger id="screenlog" spec="Logger" logEvery="{screen_every}">
            <log idref="posterior"/>
        </logger>
    </run>
</beast>
"""


def build_xml(sequences, dates, chain_length, log_every, run_id, constant_sites=None):
    """HKY+G4, strict clock, constant-size coalescent. With a tip date for every sample the
    clock rate is estimated (time-scaled tree); otherwise it is fixed at 1 and branch
    lengths are substitutions/site."""
    names = [n for n, _ in sequences]
    dated = bool(dates) and all(n in dates for n in names)
    if dates and not dated:
        missing = [n for n in names if n not in dates]
        print(f"WARNING: {len(missing)} of {len(names)} samples have no date; running undated "
              f"(clock fixed). First missing: {', '.join(missing[:5])}", file=sys.stderr)

    seqs = "\n".join(f'        <sequence id="seq_{n}" spec="Sequence" taxon="{n}" totalcount="4" value="{q}"/>'
                     for n, q in sequences)
    filtered = (f'    <data id="alignmentFiltered" spec="FilteredAlignment" filter="-" data="@alignment" '
                f'constantSiteWeights="{constant_sites}"/>\n') if constant_sites else ""
    if dated:
        values = ",".join(f"{n}={dates[n]:.4f}" for n in names)
        taxa = (f'                <trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet" traitname="date" value="{values}">\n'
                f'                    <taxa id="TaxonSet" spec="TaxonSet"><alignment idref="alignment"/></taxa>\n'
                f'                </trait>\n'
                f'                <taxonset idref="TaxonSet"/>')
    else:
        taxa = '                <taxonset id="TaxonSet" spec="TaxonSet"><alignment idref="alignment"/></taxonset>'
    return TEMPLATE.format(
        run_id=run_id, namespace=NAMESPACE, seqs=seqs, filtered=filtered, chain_length=chain_length,
        taxa=taxa, clock_est="true" if dated else "false", clock_value="1e-6" if dated else "1.0",
        clock_prior=CLOCK_PRIOR.format(m=CLOCK_PRIOR_M, s=CLOCK_PRIOR_S) if dated else "",
        clock_ops=CLOCK_OPERATORS if dated else "",
        data="@alignmentFiltered" if constant_sites else "@alignment",
        log_every=log_every, screen_every=max(log_every * 100, 10000))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    sequences = read_fasta(args.alignment)
    sample_names = [name for name, _ in sequences]

    dates = load_dates(args.metadata, sample_names)

    # Derive run_id from output filename stem
    run_id = os.path.splitext(os.path.basename(args.output))[0]

    constant = None
    if args.constant_sites:
        with open(args.constant_sites) as fh:
            constant = ' '.join(fh.read().split())
    xml_content = build_xml(sequences, dates, args.chain_length, args.log_every, run_id, constant)

    with open(args.output, "w") as fh:
        fh.write(xml_content)

    n_dated = len([s for s in sample_names if s in dates])
    print(f"INFO: Wrote {args.output} | sequences={len(sequences)} | dated={n_dated}/{len(sequences)} | chain={args.chain_length} | logEvery={args.log_every}", file=sys.stderr)


if __name__ == "__main__":
    main()
