"""Every boolean/numeric schema param must be typed in main.nf, else `--flag` arrives as a string."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NF_TYPE = {'boolean': 'Boolean', 'integer': 'Integer', 'number': 'Float'}


def test_non_string_params_are_typed_in_main_nf():
    block = re.search(r'^params \{(.*?)^\}', (ROOT / 'main.nf').read_text(), re.S | re.M).group(1)
    declared = dict(re.findall(r'^\s*(\w+):\s*(\w+)\??', block, re.M))
    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())
    for group in schema['$defs'].values():
        for name, spec in group['properties'].items():
            if isinstance(spec.get('type'), str) and spec['type'] in NF_TYPE:  # help is boolean|string
                assert declared.get(name) == NF_TYPE[spec['type']], name
