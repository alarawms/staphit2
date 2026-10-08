"""Tests for bin/staphit-sample-status."""
import csv
import importlib.machinery
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'staphit-sample-status'
_loader = importlib.machinery.SourceFileLoader('sample_status', str(BIN))
_spec = importlib.util.spec_from_loader('sample_status', _loader)
ss = importlib.util.module_from_spec(_spec)
_loader.exec_module(ss)


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_every_drop_reason(tmp_path):
    stages = _write(tmp_path / 'stages.tsv', '\n'.join([
        *[f'input\t{s}' for s in ['OK', 'NOASM', 'SMALL', 'LOWQC', 'NOTSA', 'NOTYPE', 'NOTREE']],
        *[f'assembled\t{s}' for s in ['OK', 'SMALL', 'LOWQC', 'NOTSA', 'NOTYPE', 'NOTREE']],
        'size\tSMALL\t120000',
        'size\tOK\t2800000', 'depth_short\tOK\t95.0', 'depth_short\tNOTREE\t23.4', 'depth_long\tNOTREE\t12.0',
        *[f'size_pass\t{s}' for s in ['OK', 'LOWQC', 'NOTSA', 'NOTYPE', 'NOTREE']],
        *[f'qc_pass\t{s}' for s in ['OK', 'NOTSA', 'NOTYPE', 'NOTREE']],
        *[f'species_pass\t{s}' for s in ['OK', 'NOTYPE', 'NOTREE']],
    ]) + '\n')
    qc = _write(tmp_path / 'qc.tsv', 'sample_id\tpass\treason\tcompleteness\tcontamination\nLOWQC\tFAIL\tcompleteness 71.2% < 90.0%\t71.2\t0.4\n')
    sp = _write(tmp_path / 'sp.tsv', 'sample_id\tani_percent\tidentified_species\tmash_identity\tmash_accession\nNOTSA\t0.00\tEnterobacter cloacae\t0.99\tGCF_1\n')
    summ = _write(tmp_path / 'summary.tsv', 'sample_id\tmlst_st\tcontigs\tn50\nOK\t8\t70\t95000\nNOTREE\t5\t1096\t5200\n')
    tree = _write(tmp_path / 'core.treefile', '(OK:0.1,(Reference:0.2,X:0.3)100:0.1);\n')
    out = tmp_path / 'status.tsv'

    flagged = tmp_path / 'summary_flagged.tsv'
    ss.main(['--stages', stages, '--qc', qc, '--species-excluded', sp, '--summary', summ, '--tree', tree,
             '--summary-out', str(flagged), '-o', str(out)])
    rows = {r['sample_id']: r for r in csv.DictReader(out.open(), delimiter='\t')}

    assert list(rows) == ['OK', 'NOASM', 'SMALL', 'LOWQC', 'NOTSA', 'NOTYPE', 'NOTREE']
    assert rows['OK']['status'] == 'included'
    assert rows['NOASM']['stage'] == 'assembly_failed'
    assert rows['SMALL']['stage'] == 'assembly_size_out_of_range' and '0.12 Mb outside 2.50-3.20 Mb' in rows['SMALL']['reason']
    assert rows['LOWQC']['stage'] == 'qc_gate_failed' and '71.2' in rows['LOWQC']['reason']
    assert rows['NOTSA']['stage'] == 'not_s_aureus' and 'Enterobacter' in rows['NOTSA']['reason']
    assert rows['NOTYPE']['stage'] == 'typing_incomplete'
    assert rows['NOTREE']['status'] == 'not_in_tree'
    # kept samples are flagged, dropped ones are not
    assert rows['OK']['flags'] == ''
    assert rows['NOTREE']['flags'] == 'low_depth_short:23x;low_depth_long:12x;contigs:1096;n50:5.2kb'
    assert rows['SMALL']['flags'] == ''
    summary = {r['sample_id']: r for r in csv.DictReader(flagged.open(), delimiter='\t')}
    assert summary['NOTREE']['qc_flags'].startswith('low_depth_short') and summary['OK']['qc_flags'] == ''
    assert summary['OK']['mlst_st'] == '8'


def test_placeholders_and_no_tree(tmp_path):
    stages = _write(tmp_path / 'stages.tsv', 'input\tA\nassembled\tA\nsize_pass\tA\nqc_pass\tA\nspecies_pass\tA\n')
    summ = _write(tmp_path / 'summary.tsv', 'sample_id\nA\n')
    no = _write(tmp_path / 'NO_QC', 'placeholder\n')
    out = tmp_path / 'status.tsv'
    ss.main(['--stages', stages, '--qc', no, '--species-excluded', no, '--summary', summ, '--tree', no, '-o', str(out)])
    rows = list(csv.DictReader(out.open(), delimiter='\t'))
    assert rows == [{'sample_id': 'A', 'status': 'included', 'stage': '', 'reason': '', 'flags': ''}]
