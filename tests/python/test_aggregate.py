"""Tests for bin/staphit-aggregate."""
import csv
import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

import pytest

# ---------------------------------------------------------------------------
# Import the script as a module (no .py extension)
# ---------------------------------------------------------------------------

TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-aggregate')

loader = importlib.machinery.SourceFileLoader('staphit_aggregate', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_aggregate', loader)
agg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agg)


# ===========================================================================
# TestParseTrimmomatic
# ===========================================================================

class TestParseTrimmomatic:

    def test_basic(self, tmp_path):
        log = tmp_path / 'trim.log'
        log.write_text(
            'TrimmomaticPE: Started with arguments:\n'
            'Input Read Pairs: 500000 Both Surviving: 480000 (96.00%) '
            'Forward Only Surviving: 10000 (2.00%) '
            'Reverse Only Surviving: 5000 (1.00%) Dropped: 5000 (1.00%)\n'
        )
        result = agg.parse_trimmomatic(str(log))
        assert result['raw_reads'] == 500000
        assert result['trimmed_reads'] == 480000
        assert abs(result['survival_rate'] - 96.0) < 0.01
        assert result['q30_rate'] == 0

    def test_no_match(self, tmp_path):
        log = tmp_path / 'empty.log'
        log.write_text('Some unrelated log output\n')
        result = agg.parse_trimmomatic(str(log))
        assert result['raw_reads'] == 0
        assert result['trimmed_reads'] == 0

    def test_zero_reads(self, tmp_path):
        log = tmp_path / 'zero.log'
        log.write_text('Input Read Pairs: 0 Both Surviving: 0 (0.00%)\n')
        result = agg.parse_trimmomatic(str(log))
        assert result['raw_reads'] == 0
        assert result['survival_rate'] == 0


# ===========================================================================
# TestParseQuast
# ===========================================================================

class TestParseQuast:

    def test_basic(self, tmp_path):
        report = tmp_path / 'transposed_report.tsv'
        report.write_text(
            'Assembly\t# contigs\tTotal length\tGC (%)\tN50\n'
            'sample1\t150\t2800000\t32.75\t45000\n'
        )
        result = agg.parse_quast(str(tmp_path))
        assert result['n50'] == 45000
        assert result['contigs'] == 150
        assert result['length'] == 2800000
        assert abs(result['gc'] - 32.75) < 0.01

    def test_missing_dir(self, tmp_path):
        result = agg.parse_quast(str(tmp_path / 'nonexistent'))
        assert result['n50'] == 0
        assert result['contigs'] == 0


# ===========================================================================
# TestParseMLST
# ===========================================================================

class TestParseMLST:

    def test_tab_delimited(self, tmp_path):
        tsv = tmp_path / 'mlst.tsv'
        tsv.write_text('sample.fasta\tsaureus\t398\tarc(2)\taroE(3)\n')
        result = agg.parse_mlst(str(tsv))
        assert result['scheme'] == 'saureus'
        assert result['st'] == '398'
        assert result['alleles'] == ['arc(2)', 'aroE(3)']

    def test_comma_delimited(self, tmp_path):
        csvf = tmp_path / 'mlst.csv'
        csvf.write_text('sample.fasta,saureus,8,arc(1),aroE(1)\n')
        result = agg.parse_mlst(str(csvf))
        assert result['scheme'] == 'saureus'
        assert result['st'] == '8'

    def test_short_row(self, tmp_path):
        tsv = tmp_path / 'mlst.tsv'
        tsv.write_text('x\ty\n')
        result = agg.parse_mlst(str(tsv))
        assert result['scheme'] == '-'


# ===========================================================================
# TestParseSpatyper
# ===========================================================================

class TestParseSpatyper:

    def test_with_header(self, tmp_path):
        tsv = tmp_path / 'spa.tsv'
        tsv.write_text('File\tRepeats\tType\nsample.fa\t08-16-02-25\tt008\n')
        result = agg.parse_spatyper(str(tsv))
        assert result['type'] == 't008'
        assert result['repeats'] == '08-16-02-25'

    def test_no_header(self, tmp_path):
        tsv = tmp_path / 'spa.tsv'
        tsv.write_text('sample.fa\t04-21-12\tt037\n')
        result = agg.parse_spatyper(str(tsv))
        assert result['type'] == 't037'
        assert result['repeats'] == '04-21-12'

    def test_empty_file(self, tmp_path):
        tsv = tmp_path / 'spa.tsv'
        tsv.write_text('')
        result = agg.parse_spatyper(str(tsv))
        assert result['type'] == '-'


# ===========================================================================
# TestParseSccmec
# ===========================================================================

class TestParseSccmec:

    def test_basic(self, tmp_path):
        tsv = tmp_path / 'sccmec.tsv'
        tsv.write_text('SCCmec_Type\tSubtype\nIVa\tIVa(2B)\n')
        result = agg.parse_sccmec(str(tsv))
        assert result['type'] == 'IVa'
        assert 'SCCmec_Type' in result['full_row']

    def test_alt_column(self, tmp_path):
        tsv = tmp_path / 'sccmec.tsv'
        tsv.write_text('type\tnotes\nII\tMRSA\n')
        result = agg.parse_sccmec(str(tsv))
        assert result['type'] == 'II'

    def test_empty(self, tmp_path):
        tsv = tmp_path / 'sccmec.tsv'
        tsv.write_text('SCCmec_Type\tSubtype\n')
        result = agg.parse_sccmec(str(tsv))
        assert result['type'] == 'ND'


# ===========================================================================
# TestParseAgr
# ===========================================================================

class TestParseAgr:

    def test_basic(self, tmp_path):
        jp = tmp_path / 'agr.json'
        jp.write_text(json.dumps({'agr_group': 'I', 'confidence': 0.98}))
        result = agg.parse_agr(str(jp))
        assert result['group'] == 'I'
        assert abs(result['confidence'] - 0.98) < 0.001

    def test_missing_file(self, tmp_path):
        result = agg.parse_agr(str(tmp_path / 'nofile.json'))
        assert result['group'] == 'ND'

    def test_empty_file(self, tmp_path):
        jp = tmp_path / 'agr.json'
        jp.write_text('')
        result = agg.parse_agr(str(jp))
        assert result['group'] == 'ND'


# ===========================================================================
# TestParseAbricate
# ===========================================================================

class TestParseAbricate:

    def _write_tab(self, path, db_name, rows):
        """Helper: write an ABRicate-style tab file."""
        fpath = path / f'sample_{db_name}.tab'
        header = '#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\n'
        lines = [header]
        for gene, cov, ident in rows:
            lines.append(f'sample.fa\tcontig1\t1\t100\t+\t{gene}\t{cov}\t{ident}\t{db_name}\tACC001\tProduct\n')
        fpath.write_text(''.join(lines))
        return fpath

    def test_all_databases(self, tmp_path):
        self._write_tab(tmp_path, 'resfinder', [('mecA', '100', '99.5')])
        self._write_tab(tmp_path, 'vfdb', [('hla', '100', '98.0'), ('lukF-PV', '95', '97.0')])
        self._write_tab(tmp_path, 'plasmidfinder', [('rep16', '100', '99.0')])

        resistance, virulence, plasmids = agg.parse_abricate(str(tmp_path))
        assert len(resistance) == 1
        assert resistance[0]['gene'] == 'mecA'
        assert len(virulence) == 2
        assert len(plasmids) == 1
        assert plasmids[0]['gene'] == 'rep16'

    def test_empty_dir(self, tmp_path):
        resistance, virulence, plasmids = agg.parse_abricate(str(tmp_path))
        assert resistance == []
        assert virulence == []
        assert plasmids == []

    def test_nonexistent_dir(self, tmp_path):
        resistance, virulence, plasmids = agg.parse_abricate(str(tmp_path / 'no_dir'))
        assert resistance == []


# ===========================================================================
# TestParseAmrfinder
# ===========================================================================

class TestParseAmrfinder:

    HEADER = (
        'Protein id\tGene symbol\tSequence name\tScope\tElement type\t'
        'Element subtype\tClass\tSubclass\tMethod\t'
        'Target length\tReference sequence length\t'
        '% Coverage of reference\t% Identity to reference\t'
        'Alignment length\tAccession of closest sequence\t'
        'Closest reference name\tHMM id\tHMM description'
    )

    def _make_row(self, gene, subtype='AMR', cls='', subcls='', cov='100', ident='99.5'):
        return (
            f'prot1\t{gene}\tSome name\tcore\tAMR\t{subtype}\t{cls}\t{subcls}\t'
            f'EXACTX\t300\t300\t{cov}\t{ident}\t300\tACC001\tRef name\t\t'
        )

    def test_amr_and_point_mutations(self, tmp_path):
        tsv = tmp_path / 'amrfinder.tsv'
        lines = [
            self.HEADER,
            self._make_row('mecA', subtype='AMR', cls='BETA-LACTAM'),
            self._make_row('gyrA_S84L', subtype='POINT', cls='QUINOLONE'),
        ]
        tsv.write_text('\n'.join(lines) + '\n')

        amr_items, point_muts = agg.parse_amrfinder(str(tsv))
        assert len(amr_items) == 2
        assert len(point_muts) == 1
        assert point_muts[0]['gene'] == 'gyrA'
        assert point_muts[0]['mutation'] == 'S84L'
        assert point_muts[0]['drug_class'] == 'Fluoroquinolones'

    def test_empty_file(self, tmp_path):
        tsv = tmp_path / 'amrfinder.tsv'
        tsv.write_text(self.HEADER + '\n')
        amr_items, point_muts = agg.parse_amrfinder(str(tsv))
        assert amr_items == []
        assert point_muts == []

    def test_multiple_point_mutations(self, tmp_path):
        tsv = tmp_path / 'amrfinder.tsv'
        lines = [
            self.HEADER,
            self._make_row('gyrA_S84L', subtype='POINT', cls='QUINOLONE'),
            self._make_row('rpoB_H481N', subtype='POINT', cls='RIFAMYCIN'),
        ]
        tsv.write_text('\n'.join(lines) + '\n')

        _, point_muts = agg.parse_amrfinder(str(tsv))
        assert len(point_muts) == 2
        classes = {m['drug_class'] for m in point_muts}
        assert 'Fluoroquinolones' in classes
        assert 'Rifampicin' in classes


# ===========================================================================
# TestBuildSummaryRow
# ===========================================================================

class TestBuildSummaryRow:

    def test_full_data(self):
        data = {
            'sample_id': 'S001',
            'qc': {'raw_reads': 1000, 'trimmed_reads': 950, 'survival_rate': 95.0, 'q30_rate': 0},
            'assembly': {'length': 2800000, 'contigs': 100, 'n50': 50000, 'gc': 32.8},
            'typing': {
                'mlst': {'scheme': 'saureus', 'st': '8', 'alleles': []},
                'spa': {'type': 't008', 'repeats': ''},
                'sccmec': {'type': 'IVa', 'full_row': {}},
                'agr': {'group': 'I', 'confidence': 0.9},
            },
            'resistance': {
                'amrfinder': [{'gene': 'mecA'}],
                'abricate': [{'gene': 'blaZ'}],
                'kma': [{'gene': 'tetM'}],
            },
            'virulence': [{'gene': 'hla'}],
            'virulence_profile': {'_summary': 'PVL-;TSST-;IEC-negative;ica:absent;hlg:absent'},
            'plasmids': [{'gene': 'rep16'}],
            'point_mutations': {
                '_summary': 'gyrA:S84L',
                'genotype_phenotype': {
                    'Fluoroquinolones': {'concordant': True, 'genotype': 'R', 'phenotype': 'Levofloxacin=R'},
                },
            },
            'metadata': {
                'infection_origin': 'blood',
                'antibiogram': [
                    {'antibiotic': 'Oxacillin', 'sir': 'R'},
                ],
            },
        }
        headers, values = agg.build_summary_row(data)
        assert headers == agg.SUMMARY_HEADERS
        assert len(headers) == len(values)
        # Check specific values by position
        d = dict(zip(headers, values))
        assert d['sample_id'] == 'S001'
        assert d['total_reads'] == 1000
        assert d['mlst_st'] == '8'
        assert d['spa_type'] == 't008'
        assert d['amrfinder_genes'] == 'mecA'
        assert d['abricate_genes'] == 'blaZ'
        assert d['kma_genes'] == 'tetM'
        assert d['point_mutations'] == 'gyrA:S84L'
        assert d['virulence_genes'] == 'hla'
        assert d['infection_origin'] == 'blood'
        assert d['ast_profile'] == 'OXA:R'

    def test_empty_data(self):
        data = {
            'sample_id': 'empty',
            'qc': {},
            'assembly': {},
            'typing': {},
            'resistance': {},
            'virulence': [],
            'virulence_profile': {},
            'plasmids': [],
            'point_mutations': {},
            'metadata': {},
        }
        headers, values = agg.build_summary_row(data)
        assert len(headers) == len(values)
        d = dict(zip(headers, values))
        assert d['sample_id'] == 'empty'
        assert d['total_reads'] == 0
        assert d['mlst_scheme'] == '-'

    def test_discordance_count(self):
        data = {
            'sample_id': 'disc',
            'qc': {}, 'assembly': {}, 'typing': {},
            'resistance': {}, 'virulence': [], 'virulence_profile': {},
            'plasmids': [], 'metadata': {},
            'point_mutations': {
                '_summary': 'gyrA:S84L',
                'genotype_phenotype': {
                    'Fluoroquinolones': {'concordant': False},
                    'Rifampicin': {'concordant': True},
                },
            },
        }
        headers, values = agg.build_summary_row(data)
        d = dict(zip(headers, values))
        assert d['geno_pheno_discordance'] == '1'


# ===========================================================================
# TestCLI
# ===========================================================================

class TestCLI:

    def test_minimal_run(self, tmp_path):
        """Run the CLI with only required args and a trim log."""
        trim_log = tmp_path / 'trim.log'
        trim_log.write_text(
            'Input Read Pairs: 100 Both Surviving: 90 (90.00%) '
            'Forward Only Surviving: 5 Reverse Only Surviving: 3 Dropped: 2\n'
        )
        outdir = tmp_path / 'out'

        result = subprocess.run(
            [
                sys.executable, TOOL_PATH,
                '--sample-id', 'TEST01',
                '--trim-log', str(trim_log),
                '--outdir', str(outdir),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert 'Aggregation complete' in result.stdout

        # Check JSON report
        report = json.loads((outdir / 'TEST01_report.json').read_text())
        assert report['sample_id'] == 'TEST01'
        assert report['qc']['raw_reads'] == 100
        assert report['qc']['trimmed_reads'] == 90

        # Check summary CSV
        lines = (outdir / 'TEST01_summary.csv').read_text().strip().split('\n')
        assert len(lines) == 2
        headers = lines[0].split('\t')
        assert 'sample_id' in headers

    def test_full_pipeline(self, tmp_path):
        """End-to-end with all input types."""
        # Trim log
        trim_log = tmp_path / 'trim.log'
        trim_log.write_text('Input Read Pairs: 5000 Both Surviving: 4500 (90.00%)\n')

        # QUAST
        quast_dir = tmp_path / 'quast'
        quast_dir.mkdir()
        (quast_dir / 'transposed_report.tsv').write_text(
            'Assembly\t# contigs\tTotal length\tGC (%)\tN50\n'
            'sample\t120\t2700000\t33.0\t60000\n'
        )

        # MLST
        mlst = tmp_path / 'mlst.tsv'
        mlst.write_text('sample.fa\tsaureus\t239\tarc(1)\taroE(3)\tglpF(1)\n')

        # Spa
        spa = tmp_path / 'spa.tsv'
        spa.write_text('sample.fa\t04-21-12\tt037\n')

        # SCCmec
        sccmec = tmp_path / 'sccmec.tsv'
        sccmec.write_text('SCCmec_Type\tSubtype\nIII\tIII(3A)\n')

        # Agr
        agr = tmp_path / 'agr.json'
        agr.write_text(json.dumps({'agr_group': 'I', 'confidence': 0.95}))

        # Abricate dir
        abr_dir = tmp_path / 'abricate'
        abr_dir.mkdir()
        (abr_dir / 'sample_resfinder.tab').write_text(
            '#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\n'
            'sample.fa\tc1\t1\t100\t+\tmecA\t100\t99.8\tresfinder\tACC\tPBP2a\n'
        )
        (abr_dir / 'sample_vfdb.tab').write_text(
            '#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\n'
            'sample.fa\tc1\t1\t100\t+\thla\t100\t98.5\tvfdb\tACC\tAlpha-hemolysin\n'
        )

        # AMRFinder
        amrf = tmp_path / 'amrfinder.tsv'
        amrf.write_text(
            'Protein id\tGene symbol\tSequence name\tScope\tElement type\t'
            'Element subtype\tClass\tSubclass\tMethod\t'
            'Target length\tReference sequence length\t'
            '% Coverage of reference\t% Identity to reference\t'
            'Alignment length\tAccession of closest sequence\t'
            'Closest reference name\tHMM id\tHMM description\n'
            'prot1\tgyrA_S84L\tGyrA\tcore\tAMR\tPOINT\tQUINOLONE\tFluoroquinolone\t'
            'POINTX\t300\t300\t100\t99.9\t300\tACC\tRef\t\t\n'
        )

        # KMA
        kma = tmp_path / 'kma.res'
        kma.write_text(
            '#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\t'
            'Query_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n'
            'tetM\t500\t100\t1000\t99.5\t100.0\t99.5\t100.0\t35.2\t480.0\t1.0e-26\n'
        )

        # Metadata
        meta = tmp_path / 'metadata.json'
        meta.write_text(json.dumps([
            {
                'sample_id': 'FULL01',
                'infection_origin': 'wound',
                'antibiogram': [
                    {'antibiotic': 'Levofloxacin', 'sir': 'R'},
                ],
            }
        ]))

        outdir = tmp_path / 'out'
        result = subprocess.run(
            [
                sys.executable, TOOL_PATH,
                '--sample-id', 'FULL01',
                '--trim-log', str(trim_log),
                '--quast-dir', str(quast_dir),
                '--mlst', str(mlst),
                '--spa', str(spa),
                '--sccmec', str(sccmec),
                '--agr', str(agr),
                '--abricate-dir', str(abr_dir),
                '--amrfinder', str(amrf),
                '--kma', str(kma),
                '--metadata', str(meta),
                '--outdir', str(outdir),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr

        report = json.loads((outdir / 'FULL01_report.json').read_text())
        assert report['sample_id'] == 'FULL01'
        assert report['qc']['raw_reads'] == 5000
        assert report['assembly']['n50'] == 60000
        assert report['typing']['mlst']['st'] == '239'
        assert report['typing']['spa']['type'] == 't037'
        assert report['typing']['sccmec']['type'] == 'III'
        assert report['typing']['agr']['group'] == 'I'
        assert len(report['resistance']['amrfinder']) == 1
        assert len(report['resistance']['abricate']) == 1
        assert len(report['resistance']['kma']) == 1
        assert report['point_mutations']['mutations'][0]['gene'] == 'gyrA'
        assert report['metadata']['infection_origin'] == 'wound'

    def test_missing_optional_inputs(self, tmp_path):
        """CLI should succeed even with no optional inputs."""
        outdir = tmp_path / 'out'
        result = subprocess.run(
            [
                sys.executable, TOOL_PATH,
                '--sample-id', 'BARE',
                '--outdir', str(outdir),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        report = json.loads((outdir / 'BARE_report.json').read_text())
        assert report['sample_id'] == 'BARE'
        assert report['qc'] == {}


# ===========================================================================
# TestAmrConsensus
# ===========================================================================

class TestAmrConsensus:

    def test_family_names_line_up_across_tools(self):
        assert agg.amr_family('blaZ_138', 'resfinder') == 'blaZ'
        assert agg.amr_family('mecA_8_NC_002745', 'kma') == 'mecA'
        assert agg.amr_family("aac(6')-aph(2'')_1_M13771", 'kma') == "aac(6')-Ie/aph(2'')-Ia"
        assert agg.amr_family("aph(3')-III_1_M26832", 'resfinder') == "aph(3')-IIIa"
        assert agg.amr_family('aadD_1_AF181950', 'resfinder') == 'aadD1'
        assert agg.amr_family('blaPC1', 'amrfinder') == 'blaZ'
        assert agg.amr_family('erm(C)', 'amrfinder') == 'erm(C)'

    def test_consensus_rules(self):
        amr = [
            {'gene': 'mecA', 'element_type': 'AMR', 'element_subtype': 'AMR'},
            {'gene': 'fosB', 'element_type': 'AMR', 'element_subtype': 'AMR'},       # AMRFinderPlus only: kept
            {'gene': 'glpT_A100V', 'element_type': 'AMR', 'element_subtype': 'POINT'},  # point mutation: not a gene
            {'gene': 'qacA', 'element_type': 'STRESS', 'element_subtype': 'BIOCIDE'},   # not AMR
        ]
        abr = [{'gene': 'mecA_1'}, {'gene': 'tet(K)_1'}, {'gene': 'erm(C)_13'}]
        kma = [{'gene': 'mecA_8_NC_002745'}, {'gene': 'tet(K)_4_U38428'}, {'gene': 'blaZ_130_AHKZ01000073'}]
        cons, disc, support = agg.amr_consensus(amr, abr, kma)
        assert cons == ['fosB', 'mecA', 'tet(K)']          # tet(K): ResFinder + KMA agree
        assert disc == ['blaZ[kma]', 'erm(C)[resfinder]']
        assert support['mecA'] == ['amrfinder', 'kma', 'resfinder']

    def test_raw_ont_kma_cannot_confirm(self):
        abr = [{'gene': 'tet(K)_1'}]
        kma = [{'gene': 'tet(K)_4_U38428'}]
        cons, disc, _ = agg.amr_consensus([], abr, kma, kma_raw_ont=True)
        assert cons == []
        assert disc == ['tet(K)[kma-ont,resfinder]']

    def test_trace_kma_hits_do_not_count(self):
        # KMA saw ~100x; mecA at 10x is trace (carry-over), tet(K) at 95x is real
        abr = [{'gene': 'tet(K)_1'}]
        kma = [{'gene': 'mecA_6_BX571856', 'depth': 10.0}, {'gene': 'tet(K)_4_U38428', 'depth': 95.0}]
        cons, disc, _ = agg.amr_consensus([], abr, kma, kma_ref_depth=100.0)
        assert cons == ['tet(K)']
        assert disc == ['mecA[kma-trace]']
