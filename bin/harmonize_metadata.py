#!/usr/bin/env python3
"""Single source of truth for harmonized metadata labels and colours.

Used by staphit-aggregate (adds year / host_group / source_group to the summary),
which bin/plot_tree.R and bin/itol_annotations.py then read — so every figure and
table shows the same categories. Edit the rules here, nowhere else.
"""
import re

UNKNOWN = 'Unknown'

# NCBI scientific name (lower case) -> common label
HOST_NAMES = {
    'homo sapiens': 'Human', 'camelus dromedarius': 'Camel', 'bos taurus': 'Cattle',
    'ovis aries': 'Sheep', 'capra hircus': 'Goat', 'gallus gallus': 'Chicken',
    'sus scrofa': 'Pig', 'equus caballus': 'Horse', 'canis lupus familiaris': 'Dog',
    'felis catus': 'Cat', 'vertebrata': 'Animal (unspecified)',
}

# Isolation source: first matching rule wins (case-insensitive regex).
# Food is kept as separate classes (meat / fish / milk-dairy), fish and milk before meat.
SOURCE_RULES = [
    ('Blood',                      r'blood|peripheral|bacteremia|septic'),
    ('Wound / skin & soft tissue', r'wound|tissue|skin|abscess|pus|ulcer|toe|burn'),
    ('Nasal swab',                 r'nasal|nose|nares|nasopharyn'),
    ('Respiratory',                r'sputum|respir|bronch|trache|lung|\bbal\b'),
    ('Urine',                      r'urin'),
    ('Other clinical site',        r'eye|ear|fluid|csf|bone|joint|catheter'),
    ('Clinical, unspecified',      r'clinical|diagnostic|hospital|patient'),
    ('Food: fish',                 r'fish|seafood|shrimp'),
    ('Food: milk/dairy',           r'milk|dairy|cheese|yogh?urt'),
    ('Food: meat',                 r'meat|beef|lamb|mutton|pork|chi[c]?ken|poultry|camel|retail'),
    ('Environment',                r'^env|environment|soil|water|surface|sewage|dust'),
]

HOST_COLORS = {
    'Human': '#636363', 'Camel': '#a6761d', 'Cattle': '#1b9e77', 'Sheep': '#7570b3',
    'Goat': '#e7298a', 'Chicken': '#e6ab02', 'Pig': '#fb9a99', 'Horse': '#8c510a',
    'Dog': '#66a61e', 'Cat': '#bf812d', 'Animal (unspecified)': '#bdbdbd', UNKNOWN: '#f0f0f0',
}
SOURCE_COLORS = {
    'Blood': '#e41a1c', 'Wound / skin & soft tissue': '#ff7f00', 'Nasal swab': '#4daf4a',
    'Respiratory': '#377eb8', 'Urine': '#ffff33', 'Other clinical site': '#a65628',
    'Clinical, unspecified': '#f781bf', 'Food: meat': '#984ea3', 'Food: fish': '#1f78b4',
    'Food: milk/dairy': '#cab2d6', 'Environment': '#66c2a5', UNKNOWN: '#f0f0f0',
}
OTHER_COLOR = '#969696'

# SCCmec display groups: best-fit '(est.)' calls collapse into their type; composite
# elements into one class; no element = MSSA; element of unclear type = Untypeable.
SCCMEC_COLORS = {
    'Type I': '#a6cee3', 'Type II': '#1f78b4', 'Type III': '#b2df8a', 'Type IV': '#4daf4a',
    'Type V': '#377eb8', 'Type VI': '#ffd92f', 'Type VII': '#ff7f00', 'Type VIII': '#fb9a99',
    'Type IX': '#e31a1c', 'Type X': '#fdbf6f', 'Type XI': '#984ea3', 'Type XII': '#cab2d6',
    'Type XIII': '#6a3d9a', 'Type XIV': '#b15928', 'Composite': '#a8d978',
    'Untypeable': '#f781bf', 'MSSA': '#d9d9d9',
}


def year(collection_date):
    m = re.match(r'\s*(\d{4})', collection_date or '')
    return m.group(1) if m else UNKNOWN


def host_group(host):
    h = (host or '').strip().lower()
    if not h:
        return UNKNOWN
    return HOST_NAMES.get(h, h[:1].upper() + h[1:])


def source_group(source):
    s = (source or '').strip().lower()
    for label, pattern in SOURCE_RULES:
        if re.search(pattern, s):
            return label
    return UNKNOWN


def sccmec_group(raw):
    v = re.sub(r'\s*\(est\.\)\s*$', '', (raw or '').strip())
    if v in ('', '-', 'Negative', 'ND', 'None'):
        return 'MSSA'
    if v.startswith('Composite'):
        return 'Composite'
    if v.lower() in ('unknown', 'untypeable', 'nt'):
        return 'Untypeable'
    return v


def year_colors(years):
    """Sequential light-yellow -> dark-red ramp over the sorted known years; Unknown grey."""
    known = sorted(y for y in set(years) if y != UNKNOWN)
    stops = [(0xff, 0xf5, 0xb1), (0xfd, 0x8d, 0x3c), (0x7f, 0x00, 0x00)]
    out = {}
    for i, y in enumerate(known):
        t = i / max(1, len(known) - 1) * (len(stops) - 1)
        a, b = stops[int(t)], stops[min(int(t) + 1, len(stops) - 1)]
        f = t - int(t)
        out[y] = '#%02x%02x%02x' % tuple(round(a[k] + (b[k] - a[k]) * f) for k in range(3))
    out[UNKNOWN] = '#f0f0f0'
    return out


if __name__ == '__main__':  # self-test
    assert year('2020-03-01') == '2020' and year('') == UNKNOWN
    assert host_group('Homo sapiens') == 'Human' and host_group('Vertebrata') == 'Animal (unspecified)'
    assert host_group('') == UNKNOWN and host_group('mus musculus') == 'Mus musculus'
    cases = {'clinical sample': 'Clinical, unspecified', 'PERIPHERAL': 'Blood', 'Chiken': 'Food: meat',
             'goat milk': 'Food: milk/dairy', 'fish': 'Food: fish', 'Goat Nasal swabs': 'Nasal swab',
             'retail meat from Camelus dromedarius': 'Food: meat', 'env': 'Environment',
             'TOE/TISSUE': 'Wound / skin & soft tissue', 'RESPIRATORY CULTURE': 'Respiratory',
             'Riyadh': UNKNOWN, '': UNKNOWN}
    for raw, want in cases.items():
        assert source_group(raw) == want, (raw, source_group(raw), want)
    yc = year_colors(['2010', '2022', '2015', UNKNOWN])
    assert yc['2010'] == '#fff5b1' and yc['2022'] == '#7f0000' and yc[UNKNOWN] == '#f0f0f0'
    assert set(SOURCE_COLORS) == {l for l, _ in SOURCE_RULES} | {UNKNOWN}
    assert sccmec_group('Type V (est.)') == 'Type V' and sccmec_group('Negative') == 'MSSA'
    assert sccmec_group('Composite (Type III)') == 'Composite' and sccmec_group('Unknown') == 'Untypeable'
    assert sccmec_group('Type XIV') == 'Type XIV'
    print('selftest ok')
