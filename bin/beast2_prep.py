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
# XML builder
# ---------------------------------------------------------------------------

INDENT = "    "


def seq_elements(sequences):
    lines = []
    for name, seq in sequences:
        lines.append(f'{INDENT}    <sequence id="seq_{name}" taxon="{name}" totalcount="4" value="{seq}"/>')
    return "\n".join(lines)


def date_trait_block(dates, sample_names):
    """Return the dateTrait XML block (empty string if no dates)."""
    if not dates:
        return ""
    # Only include samples that have dates
    pairs = []
    for name in sample_names:
        if name in dates:
            pairs.append(f"{name}={dates[name]:.6f}")
    if not pairs:
        return ""
    value_str = ",".join(pairs)
    return (
        f'\n    <trait id="dateTrait" spec="beast.evolution.tree.TraitSet" traitname="date-forward"\n'
        f'           value="{value_str}">\n'
        f'        <taxa id="TaxonSet.alignment" spec="TaxonSet">\n'
        f'            <alignment idref="alignment"/>\n'
        f'        </taxa>\n'
        f'    </trait>'
    )


def tree_block(has_dates):
    if has_dates:
        return (
            '\n    <tree id="Tree.t:alignment" spec="beast.evolution.tree.Tree" name="stateNode">\n'
            '        <trait idref="dateTrait"/>\n'
            '        <taxonset idref="TaxonSet.alignment"/>\n'
            '    </tree>'
        )
    else:
        return (
            '\n    <tree id="Tree.t:alignment" spec="beast.evolution.tree.Tree" name="stateNode">\n'
            '        <taxonset idref="TaxonSet.alignment"/>\n'
            '    </tree>'
        )


def build_xml(sequences, dates, chain_length, log_every, run_id):
    sample_names = [name for name, _ in sequences]
    has_dates = bool(dates)

    seq_block   = seq_elements(sequences)
    trait_block = date_trait_block(dates, sample_names) if has_dates else ""
    t_block     = tree_block(has_dates)

    xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="" version="2.7">

    <data id="alignment" dataType="nucleotide" name="alignment">
{seq_block}
    </data>
{trait_block}{t_block}

    <treeModel id="treeModel" spec="TreeModel" tree="@Tree.t:alignment"/>

    <siteModel id="SiteModel.s:alignment" spec="SiteModel">
        <shape id="gammaShape.s:alignment" spec="parameter.RealParameter" value="0.5"/>
        <gammaCategoryCount spec="parameter.IntegerParameter" value="4"/>
        <substModel id="hky.s:alignment" spec="HKY" kappa="@kappa.s:alignment">
            <frequencies id="freqParameter.s:alignment" spec="parameter.RealParameter"
                         dimension="4" lower="0.0" upper="1.0" value="0.25"/>
        </substModel>
    </siteModel>

    <parameter id="kappa.s:alignment" spec="parameter.RealParameter" lower="0.0" value="2.0"/>
    <parameter id="clockRate.c:alignment" spec="parameter.RealParameter" value="1e-3"/>
    <parameter id="popSize.t:alignment" spec="parameter.RealParameter" lower="0.0" value="0.3"/>

    <treeLikelihood id="treeLikelihood.alignment" spec="TreeLikelihood"
                    data="@alignment" tree="@Tree.t:alignment" siteModel="@SiteModel.s:alignment">
        <branchRateModel id="StrictClock.c:alignment" spec="beast.evolution.branchratemodel.StrictClockModel"
                         clock.rate="@clockRate.c:alignment"/>
    </treeLikelihood>

    <coalescentLikelihood id="CoalescentConstant.t:alignment"
                          spec="beast.evolution.tree.coalescent.ConstantPopulation"
                          popSize="@popSize.t:alignment" tree="@Tree.t:alignment"/>

    <run id="mcmc" spec="MCMC" chainLength="{chain_length}">
        <state id="state" storeEvery="5000">
            <tree idref="Tree.t:alignment"/>
            <parameter idref="kappa.s:alignment"/>
            <parameter idref="clockRate.c:alignment"/>
            <parameter idref="popSize.t:alignment"/>
        </state>

        <distribution id="posterior" spec="util.CompoundDistribution">
            <distribution id="prior" spec="util.CompoundDistribution">
                <prior id="KappaPrior" name="distribution" x="@kappa.s:alignment">
                    <LogNormal meanInRealSpace="false" M="1.0" S="1.25"/>
                </prior>
                <prior id="ClockRatePrior" name="distribution" x="@clockRate.c:alignment">
                    <Uniform lower="0.0" upper="1e-2"/>
                </prior>
                <prior id="PopSizePrior" name="distribution" x="@popSize.t:alignment">
                    <OneOnX/>
                </prior>
                <distribution idref="CoalescentConstant.t:alignment"/>
            </distribution>
            <distribution id="likelihood" idref="treeLikelihood.alignment"/>
        </distribution>

        <!-- Operators -->
        <operator id="treeScaler" spec="ScaleOperator" scaleFactor="0.5"
                  tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="treeRootScaler" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5"
                  tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="UniformOperator" spec="Uniform" tree="@Tree.t:alignment" weight="30.0"/>
        <operator id="subtreeSlide" spec="SubtreeSlide" tree="@Tree.t:alignment" weight="15.0"/>
        <operator id="narrow" spec="Exchange" tree="@Tree.t:alignment" weight="15.0"/>
        <operator id="wide" spec="Exchange" isNarrow="false" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="WilsonBalding" spec="WilsonBalding" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="kappaScaler" spec="ScaleOperator" parameter="@kappa.s:alignment"
                  scaleFactor="0.5" weight="0.1"/>
        <operator id="clockRateScaler" spec="ScaleOperator" parameter="@clockRate.c:alignment"
                  scaleFactor="0.5" weight="3.0"/>
        <operator id="popSizeScaler" spec="ScaleOperator" parameter="@popSize.t:alignment"
                  scaleFactor="0.75" weight="3.0"/>

        <!-- Loggers -->
        <logger id="tracelog" spec="Logger" fileName="beast_{run_id}.log" logEvery="{log_every}">
            <log idref="posterior"/>
            <log idref="prior"/>
            <log idref="likelihood"/>
            <parameter idref="kappa.s:alignment"/>
            <parameter idref="clockRate.c:alignment"/>
            <parameter idref="popSize.t:alignment"/>
            <log spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:alignment"/>
        </logger>
        <logger id="treelog" spec="Logger" fileName="beast_{run_id}.trees" logEvery="{log_every}"
                mode="tree">
            <log spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:alignment"/>
        </logger>
        <logger id="screenlog" spec="Logger" logEvery="100000">
            <log idref="posterior"/>
            <log spec="beast.util.ESS" arg="@posterior"/>
        </logger>
    </run>
</beast>
"""
    return xml


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

    xml_content = build_xml(sequences, dates, args.chain_length, args.log_every, run_id)

    with open(args.output, "w") as fh:
        fh.write(xml_content)

    n_dated = len([s for s in sample_names if s in dates])
    print(f"INFO: Wrote {args.output} | sequences={len(sequences)} | dated={n_dated}/{len(sequences)} | chain={args.chain_length} | logEvery={args.log_every}", file=sys.stderr)


if __name__ == "__main__":
    main()
