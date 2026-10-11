"""Tests for bin/beast2_prep.py (BEAST 2.7 XML)."""
import importlib.machinery
import importlib.util
import xml.etree.ElementTree as ET
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'beast2_prep.py'
_loader = importlib.machinery.SourceFileLoader('beast2_prep', str(BIN))
_spec = importlib.util.spec_from_loader('beast2_prep', _loader)
bp = importlib.util.module_from_spec(_spec)
_loader.exec_module(bp)

SEQS = [('s1', 'ACGT'), ('s2', 'ACGA'), ('s3', 'TCGA')]


def test_dated_xml_is_time_scaled_with_constant_sites():
    xml = bp.build_xml(SEQS, {'s1': 2019.5, 's2': 2020.1, 's3': 2021.9}, 1000, 10, 'beast_run', '700 400 400 700')
    root = ET.fromstring(xml)                          # well-formed
    assert 'beast.base.inference' in root.get('namespace')   # BEAST 2.7 class names
    filt = root.find(".//data[@id='alignmentFiltered']")
    assert filt.get('constantSiteWeights') == '700 400 400 700'
    assert root.find(".//*[@id='treeLikelihood']").get('data') == '@alignmentFiltered'
    assert root.find(".//trait[@id='dateTrait']") is not None
    assert root.find(".//parameter[@id='clockRate']").get('estimate') == 'true'
    assert root.find(".//operator[@id='clockUpDown']") is not None
    assert {l.get('fileName') for l in root.iter('logger') if l.get('fileName')} == {'beast_run.log', 'beast_run.trees'}
    # every element names its class (no BEAUti <map> shortcuts are defined)
    for el in root.iter():
        if el.tag not in ('beast', 'sequence') and el.get('idref') is None:
            assert el.get('spec'), el.tag


def test_undated_or_partly_dated_fixes_the_clock():
    for dates in ({}, {'s1': 2019.5}):
        root = ET.fromstring(bp.build_xml(SEQS, dates, 1000, 10, 'beast_run'))
        assert root.find(".//trait") is None
        assert root.find(".//parameter[@id='clockRate']").get('estimate') == 'false'
        assert root.find(".//data[@id='alignmentFiltered']") is None
