#!/usr/bin/env python3
"""Fetch standardized sample metadata from NCBI BioSample for staphit2.

Input: any file containing BioSample accessions (SAMN/SAMEA/SAMD) — a staphit2
samplesheet, a manifest, a plain list. For a samplesheet (has a 'sample' column),
one metadata row is written per sample (merged samples like SAMEA1_SAMEA2 combine
their BioSamples); otherwise one row per BioSample.

NCBI BioSample mirrors all INSDC BioSamples (incl. ENA's SAMEA) and maps submitter
attributes to harmonized names (collection_date, geo_loc_name, host, ...). These are
mapped onto the pipeline's PHA4GE-style columns (see bin/staphit-metadata):
  - dates normalized to ISO 8601 (YYYY, YYYY-MM or YYYY-MM-DD)
  - geo_loc_name "Country: region" split into geo_loc_country / geo_loc_region
  - human host names normalized to "Homo sapiens"
  - INSDC null terms ("missing", "not collected", ...) left empty

Outputs:
  <out.csv>                    metadata, ready for --metadata
  <out>.completeness.tsv       per-field fill rate vs. the minimum required fields

Usage: fetch_metadata.py <file with accessions> <metadata.csv>
Set NCBI_API_KEY to raise NCBI's rate limit (3 -> 10 requests/s).
"""
import csv
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSAMPLE = re.compile(r"(?<![A-Za-z0-9])SAM(?:N|EA|D)\d+(?![0-9])")
RUN = re.compile(r"(?<![A-Za-z0-9])[EDS]RR\d{6,}(?![0-9])")

# Pipeline metadata columns (bin/staphit-metadata METADATA_COLUMNS / REQUIRED_METADATA)
COLUMNS = [
    'sample_id', 'organism', 'collection_date', 'geo_loc_country', 'geo_loc_region',
    'host', 'isolation_source', 'collected_by', 'host_age', 'host_sex',
    'host_disease', 'host_body_site', 'patient_status', 'infection_origin',
    'mrsa_status', 'strain', 'lat_lon', 'sequenced_by', 'sequencing_platform',
    'pvl_status', 'biosample_accession', 'sra_accession', 'outbreak_id',
    'travel_history', 'icu_admission',
]
REQUIRED = ['sample_id', 'organism', 'collection_date', 'geo_loc_country',
            'geo_loc_region', 'host', 'isolation_source']

# pipeline column <- NCBI harmonized names, in order of preference
FIELD_MAP = {
    'collection_date':     ['collection_date'],
    'host':                ['host'],
    'isolation_source':    ['isolation_source', 'host_tissue_sampled', 'env_medium'],
    'collected_by':        ['collected_by'],
    'host_age':            ['host_age'],
    'host_sex':            ['host_sex'],
    'host_disease':        ['host_disease'],
    'host_body_site':      ['host_body_site', 'host_tissue_sampled'],
    'strain':              ['strain', 'isolate'],
    'lat_lon':             ['lat_lon'],
    'sequencing_platform': ['platform_model', 'seq_methods', 'sequencing_meth'],
}
NULL_TERMS = {'', 'missing', 'not collected', 'not provided', 'not applicable',
              'restricted access', 'unknown', 'na', 'n/a', 'none', '-', 'null',
              'not available', 'not determined', 'unspecified', 'missing: not provided',
              'missing: not collected', 'missing: control sample',
              'not avaliable', 'not availble', 'n.a.'}          # common typos seen in INSDC
HUMAN = {'homo sapiens', 'human', 'humans', 'homo sapien', 'h. sapiens', 'patient',
         'homo sapinens', 'homo sapines', 'homo spaiens'}
# Values that say nothing about the field they are in
UNINFORMATIVE = {'isolation_source': {'host', 'sample', 'isolate', 'bacteria', 'culture'}}


# ── Normalization ────────────────────────────────────────────────────────────
def clean(value):
    v = (value or '').strip()
    return '' if v.lower() in NULL_TERMS else v


def norm_date(value):
    """ISO 8601 date (YYYY, YYYY-MM, YYYY-MM-DD) or '' if unparseable."""
    v = clean(value)
    if not v:
        return ''
    v = v.split('T')[0].split('/')[0].strip()          # drop time; first date of a range
    if re.fullmatch(r'\d{4}(-\d{2}(-\d{2})?)?', v):
        return v
    for fmt, out in (('%d-%b-%Y', '%Y-%m-%d'), ('%b-%Y', '%Y-%m'), ('%d/%m/%Y', '%Y-%m-%d'),
                     ('%Y/%m/%d', '%Y-%m-%d'), ('%B %Y', '%Y-%m'), ('%d %B %Y', '%Y-%m-%d'),
                     ('%Y%m%d', '%Y-%m-%d')):
        try:
            return datetime.strptime(v, fmt).strftime(out)
        except ValueError:
            pass
    m = re.search(r'\b(19|20)\d{2}\b', v)               # last resort: a year
    return m.group(0) if m else ''


def split_geo(value):
    """'Saudi Arabia: Riyadh, X' -> ('Saudi Arabia', 'Riyadh, X')."""
    v = clean(value)
    if not v:
        return '', ''
    country, _, region = v.partition(':')
    return country.strip(), clean(region)


def norm_host(value):
    v = clean(value)
    return 'Homo sapiens' if v.lower() in HUMAN else v


def record_to_row(rec):
    """One BioSample record (dict) -> pipeline metadata row (without sample_id)."""
    attrs = rec['attributes']
    row = {c: '' for c in COLUMNS}
    for col, names in FIELD_MAP.items():
        row[col] = next((clean(attrs[n]) for n in names if clean(attrs.get(n))), '')
        if row[col].lower() in UNINFORMATIVE.get(col, ()):
            row[col] = ''
    row['collection_date'] = norm_date(row['collection_date'])
    row['host'] = norm_host(row['host'])
    row['geo_loc_country'], row['geo_loc_region'] = split_geo(attrs.get('geo_loc_name'))
    row['organism'] = rec['organism']
    row['biosample_accession'] = rec['accession']
    return row


def merge_rows(rows):
    """Combine rows from several BioSamples of one sample: first non-empty value wins;
    accession fields are joined."""
    out = {c: '' for c in COLUMNS}
    for c in COLUMNS:
        vals = [r[c] for r in rows if r[c]]
        if c in ('biosample_accession', 'sra_accession'):
            out[c] = ';'.join(dict.fromkeys(vals))
        elif vals:
            out[c] = vals[0]
    return out


# ── NCBI fetch ───────────────────────────────────────────────────────────────
def _get(url, params, tries=6):
    if os.environ.get('NCBI_API_KEY'):
        params = dict(params, api_key=os.environ['NCBI_API_KEY'])
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, data=data, timeout=120) as r:
                return r.read()
        except Exception as e:
            if attempt == tries:
                raise RuntimeError(f"NCBI request failed: {e}")
            time.sleep(min(60, 3 * 2 ** attempt))


def parse_biosample_xml(xml_bytes):
    """BioSampleSet XML -> list of {accession, organism, attributes{harmonized: value}}."""
    out = []
    for bs in ET.fromstring(xml_bytes).iter('BioSample'):
        acc = bs.get('accession', '')
        org = bs.find('Description/Organism')
        attrs = {}
        for a in bs.iter('Attribute'):
            key = a.get('harmonized_name') or a.get('attribute_name', '').strip().lower().replace(' ', '_')
            attrs.setdefault(key, (a.text or '').strip())
        out.append({'accession': acc,
                    'organism': org.get('taxonomy_name', '') if org is not None else '',
                    'attributes': attrs})
    return out


def fetch_biosamples(accessions, batch=150):
    """{accession: record} for all accessions (missing ones are simply absent)."""
    found = {}
    accs = sorted(set(accessions))
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        term = ' OR '.join(f'{a}[accn]' for a in chunk)
        ids = re.findall(rb'<Id>(\d+)</Id>', _get(f'{EUTILS}/esearch.fcgi',
                         {'db': 'biosample', 'term': term, 'retmax': len(chunk) * 2}))
        if ids:
            xml = _get(f'{EUTILS}/efetch.fcgi',
                       {'db': 'biosample', 'id': ','.join(i.decode() for i in ids), 'retmode': 'xml'})
            for rec in parse_biosample_xml(xml):
                found[rec['accession']] = rec
        print(f'\rfetched {min(i + batch, len(accs))}/{len(accs)} BioSamples', end='', file=sys.stderr, flush=True)
        time.sleep(0.4)
    print(file=sys.stderr)
    return found


# ── Input & reporting ────────────────────────────────────────────────────────
def samples_from_file(path):
    """[(sample_id, [biosamples], [runs])]. Samplesheet -> per 'sample'; else per BioSample."""
    with open(path, newline='') as f:
        text = f.read()
    first = text.splitlines()[0] if text else ''
    if first.split(',')[0].strip() == 'sample':
        out = []
        for r in csv.DictReader(text.splitlines()):
            out.append((r['sample'], BIOSAMPLE.findall(r['sample']),
                        RUN.findall(' '.join(r.get(c) or '' for c in ('fastq_1', 'fastq_2', 'long_fastq')))))
        return out
    return [(a, [a], []) for a in dict.fromkeys(BIOSAMPLE.findall(text))]


def completeness(rows):
    n = len(rows) or 1
    lines = ['field\trequired\tfilled\tpercent']
    for c in COLUMNS:
        k = sum(1 for r in rows if r[c])
        lines.append(f"{c}\t{'yes' if c in REQUIRED else 'no'}\t{k}\t{100 * k / n:.1f}")
    full = sum(1 for r in rows if all(r[c] for c in REQUIRED))
    return '\n'.join(lines) + '\n', full


def main(src, dst):
    samples = samples_from_file(src)
    records = fetch_biosamples([b for _, bs, _ in samples for b in bs])
    rows, not_found = [], []
    for sid, bss, runs in samples:
        parts = [record_to_row(records[b]) for b in bss if b in records]
        not_found += [b for b in bss if b not in records]
        row = merge_rows(parts) if parts else {c: '' for c in COLUMNS}
        row['sample_id'] = sid
        row['sra_accession'] = row['sra_accession'] or ';'.join(dict.fromkeys(runs))
        rows.append(row)
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    report, full = completeness(rows)
    rep_path = re.sub(r'\.csv$', '', dst) + '.completeness.tsv'
    open(rep_path, 'w').write(report)
    for b in not_found:
        print(f'NOT FOUND in BioSample: {b}', file=sys.stderr)
    print(report, end='')
    print(f'{len(rows)} samples -> {dst}; {full} have all {len(REQUIRED)} required fields; report: {rep_path}')


def _selftest():
    assert norm_date('2020-03-01T10:00:00Z') == '2020-03-01'
    assert norm_date('Mar-2019') == '2019-03'
    assert norm_date('01-Mar-2019') == '2019-03-01'
    assert norm_date('2018/2019') == '2018'
    assert norm_date('not collected') == '' and norm_date('missing') == ''
    assert split_geo('Saudi Arabia: Riyadh, King Fahad') == ('Saudi Arabia', 'Riyadh, King Fahad')
    assert split_geo('Saudi Arabia') == ('Saudi Arabia', '')
    assert split_geo('not provided') == ('', '')
    assert norm_host('human') == 'Homo sapiens' and norm_host('Homo sapinens') == 'Homo sapiens'
    assert clean('Not Avaliable') == ''
    xml = b'''<BioSampleSet><BioSample accession="SAMEA1"><Description>
      <Organism taxonomy_name="Staphylococcus aureus"/></Description><Attributes>
      <Attribute attribute_name="collection date" harmonized_name="collection_date">2020-03-01</Attribute>
      <Attribute attribute_name="geographic location (country and/or sea)" harmonized_name="geo_loc_name">Saudi Arabia: Qatif</Attribute>
      <Attribute attribute_name="host scientific name" harmonized_name="host">human</Attribute>
      <Attribute attribute_name="isolation_source" harmonized_name="isolation_source">blood</Attribute>
      <Attribute attribute_name="host disease">missing</Attribute>
    </Attributes></BioSample></BioSampleSet>'''
    rec = parse_biosample_xml(xml)[0]
    row = record_to_row(rec)
    assert row['collection_date'] == '2020-03-01' and row['geo_loc_country'] == 'Saudi Arabia'
    assert row['geo_loc_region'] == 'Qatif' and row['host'] == 'Homo sapiens'
    assert row['isolation_source'] == 'blood' and row['host_disease'] == ''
    rec['attributes']['isolation_source'] = 'Host'
    assert record_to_row(rec)['isolation_source'] == ''
    assert row['organism'] == 'Staphylococcus aureus' and row['biosample_accession'] == 'SAMEA1'
    other = dict(row, collection_date='', host_disease='sepsis', biosample_accession='SAMEA2')
    m = merge_rows([dict(row, host_disease=''), other])
    assert m['host_disease'] == 'sepsis' and m['biosample_accession'] == 'SAMEA1;SAMEA2'
    print('selftest ok')


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        _selftest()
    elif len(sys.argv) == 3:
        main(*sys.argv[1:])
    else:
        sys.exit(__doc__)
