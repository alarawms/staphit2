"""Tests for staphit-fetch: SRA search and metadata bridge."""
import csv
import importlib.machinery
import importlib.util
import io
import os
import sys
import tempfile

import pytest

# Import the module directly from the bin script (no .py extension)
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'bin')
_loader = importlib.machinery.SourceFileLoader('staphit_fetch', os.path.join(BIN_DIR, 'staphit-fetch'))
_spec = importlib.util.spec_from_loader('staphit_fetch', _loader)
fetch = importlib.util.module_from_spec(_spec)
_loader.exec_module(fetch)


# ---------------------------------------------------------------------------
# Mock data
# ---------------------------------------------------------------------------

MOCK_RUNINFO_CSV = """\
Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash,Organism,strain,geo_loc_name_country_calc,collection_date_sam,Host,host_disease,Isolation_source
SRR12345678,2023-01-15,2023-01-14,5000000,1500000000,5000000,300,500,,https://example.com,SRX123,lib1,WGS,RANDOM,GENOMIC,PAIRED,350,50,ILLUMINA,Illumina NovaSeq 6000,SRP123,PRJNA123456,,123456,SRS123,SAMN12345678,simple,1280,Staphylococcus aureus,SA001,,,,,,,,,,,,CENTER1,SUB123,,public,abc123,def456,Staphylococcus aureus,MRSA-SA001,Saudi Arabia: Riyadh,2023-01-10,Homo sapiens,bacteremia,blood culture
SRR12345679,2023-02-20,2023-02-19,6000000,1800000000,6000000,300,600,,https://example.com,SRX124,lib2,WGS,RANDOM,GENOMIC,PAIRED,350,50,ILLUMINA,Illumina NovaSeq 6000,SRP123,PRJNA123456,,123456,SRS124,SAMN12345679,simple,1280,Staphylococcus aureus,SA002,,,,,,,,,,,,CENTER1,SUB123,,public,ghi789,jkl012,Staphylococcus aureus,,Saudi Arabia: Jeddah,2023-02-15,Homo sapiens,,wound swab
SRR12345680,2023-03-10,2023-03-09,4500000,1350000000,4500000,300,450,,https://example.com,SRX125,lib3,WGS,RANDOM,GENOMIC,PAIRED,350,50,ILLUMINA,Illumina MiSeq,SRP456,PRJNA789012,,789012,SRS125,SAMN12345680,simple,1280,Staphylococcus aureus,SA003,,,,,,,,,,,,CENTER2,SUB456,,public,mno345,pqr678,Staphylococcus aureus,ST239,,2023-03-05,,MRSA infection,nasal swab
"""

MOCK_RECORDS = [
    {
        'run_accession': 'SRR12345678',
        'biosample': 'SAMN12345678',
        'bioproject': 'PRJNA123456',
        'organism': 'Staphylococcus aureus',
        'geo_loc_name': 'Saudi Arabia: Riyadh',
        'collection_date': '2023-01-10',
        'host': 'Homo sapiens',
        'isolation_source': 'blood culture',
        'host_disease': 'bacteremia',
        'strain': 'MRSA-SA001',
        'spots': '5000000',
        'bases': '1500000000',
        'platform': 'ILLUMINA',
    },
    {
        'run_accession': 'SRR12345679',
        'biosample': 'SAMN12345679',
        'bioproject': 'PRJNA123456',
        'organism': 'Staphylococcus aureus',
        'geo_loc_name': 'Saudi Arabia: Jeddah',
        'collection_date': '2023-02-15',
        'host': 'Homo sapiens',
        'isolation_source': 'wound swab',
        'host_disease': '',
        'strain': '',
        'spots': '6000000',
        'bases': '1800000000',
        'platform': 'ILLUMINA',
    },
    {
        'run_accession': 'SRR12345680',
        'biosample': 'SAMN12345680',
        'bioproject': 'PRJNA789012',
        'organism': 'Staphylococcus aureus',
        'geo_loc_name': '',
        'collection_date': '2023-03-05',
        'host': '',
        'isolation_source': 'nasal swab',
        'host_disease': 'MRSA infection',
        'strain': 'ST239',
        'spots': '4500000',
        'bases': '1350000000',
        'platform': 'ILLUMINA',
    },
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestParseRuninfo:
    """Test parsing of SRA runinfo CSV into structured records."""

    def test_parse_runinfo_basic(self):
        """Parse mock runinfo CSV and check record count."""
        records = fetch.parse_runinfo_csv(MOCK_RUNINFO_CSV)
        assert len(records) == 3

    def test_parse_runinfo_fields(self):
        """Check that key fields are correctly extracted."""
        records = fetch.parse_runinfo_csv(MOCK_RUNINFO_CSV)
        rec = records[0]
        assert rec['run_accession'] == 'SRR12345678'
        assert rec['biosample'] == 'SAMN12345678'
        assert rec['bioproject'] == 'PRJNA123456'
        assert rec['organism'] == 'Staphylococcus aureus'
        assert rec['platform'] == 'ILLUMINA'
        assert rec['spots'] == '5000000'

    def test_parse_runinfo_geo_loc(self):
        """Check geo_loc_name extraction including country:region format."""
        records = fetch.parse_runinfo_csv(MOCK_RUNINFO_CSV)
        assert records[0]['geo_loc_name'] == 'Saudi Arabia: Riyadh'
        assert records[1]['geo_loc_name'] == 'Saudi Arabia: Jeddah'
        assert records[2]['geo_loc_name'] == ''

    def test_parse_runinfo_collection_date(self):
        """Check collection_date extraction."""
        records = fetch.parse_runinfo_csv(MOCK_RUNINFO_CSV)
        assert records[0]['collection_date'] == '2023-01-10'
        assert records[2]['collection_date'] == '2023-03-05'

    def test_parse_runinfo_empty(self):
        """Parsing empty CSV returns empty list."""
        records = fetch.parse_runinfo_csv('')
        assert records == []

    def test_parse_runinfo_header_only(self):
        """Parsing CSV with only headers returns empty list."""
        header = MOCK_RUNINFO_CSV.split('\n')[0] + '\n'
        records = fetch.parse_runinfo_csv(header)
        assert records == []

    def test_parse_runinfo_all_output_columns(self):
        """Every record has all OUTPUT_COLUMNS keys."""
        records = fetch.parse_runinfo_csv(MOCK_RUNINFO_CSV)
        for rec in records:
            for col in fetch.OUTPUT_COLUMNS:
                assert col in rec, f'Missing column: {col}'


class TestMetadataMapping:
    """Test mapping from SRA fields to PHA4GE schema."""

    def test_map_basic(self):
        """Basic PHA4GE mapping preserves key fields."""
        mapped = fetch.map_to_pha4ge(MOCK_RECORDS)
        assert len(mapped) == 3
        row = mapped[0]
        assert row['sra_accession'] == 'SRR12345678'
        assert row['biosample_accession'] == 'SAMN12345678'
        assert row['organism'] == 'Staphylococcus aureus'
        assert row['host'] == 'Homo sapiens'

    def test_map_geo_loc_parsed(self):
        """geo_loc_name is split into country and region."""
        mapped = fetch.map_to_pha4ge(MOCK_RECORDS)
        assert mapped[0]['geo_loc_country'] == 'Saudi Arabia'
        assert mapped[0]['geo_loc_region'] == 'Riyadh'

    def test_map_sample_id(self):
        """run_accession is used as sample_id."""
        mapped = fetch.map_to_pha4ge(MOCK_RECORDS)
        assert mapped[0]['sample_id'] == 'SRR12345678'
        assert mapped[2]['sample_id'] == 'SRR12345680'

    def test_map_empty_fields(self):
        """Empty SRA fields map to empty strings, not None."""
        mapped = fetch.map_to_pha4ge(MOCK_RECORDS)
        # Record 1 has no strain
        assert mapped[1]['strain'] == ''
        # Record 2 has no host
        assert mapped[2]['host'] == ''

    def test_map_all_records(self):
        """All records are mapped."""
        mapped = fetch.map_to_pha4ge(MOCK_RECORDS)
        accessions = [r['sra_accession'] for r in mapped]
        assert accessions == ['SRR12345678', 'SRR12345679', 'SRR12345680']


class TestPreviewCompleteness:
    """Test completeness statistics calculation."""

    def test_completeness_full_field(self):
        """Field present in all records shows 100%."""
        stats = fetch.calculate_completeness(MOCK_RECORDS)
        assert stats['run_accession']['pct'] == 100.0
        assert stats['run_accession']['count'] == 3

    def test_completeness_partial_field(self):
        """Field present in some records shows correct percentage."""
        stats = fetch.calculate_completeness(MOCK_RECORDS)
        # strain: 'MRSA-SA001', '', 'ST239' -> 2/3
        assert stats['strain']['count'] == 2
        assert stats['strain']['pct'] == pytest.approx(66.7, abs=0.1)

    def test_completeness_empty_field(self):
        """geo_loc_name has 2 of 3 filled."""
        stats = fetch.calculate_completeness(MOCK_RECORDS)
        assert stats['geo_loc_name']['count'] == 2

    def test_completeness_empty_records(self):
        """Empty records list returns empty stats."""
        stats = fetch.calculate_completeness([])
        assert stats == {}

    def test_completeness_specific_fields(self):
        """Can calculate completeness for a subset of fields."""
        stats = fetch.calculate_completeness(MOCK_RECORDS, fields=['run_accession', 'strain'])
        assert len(stats) == 2
        assert 'run_accession' in stats
        assert 'strain' in stats


class TestGenerateSamplesheet:
    """Test samplesheet CSV generation."""

    def test_generate_samplesheet_format(self):
        """Samplesheet has correct columns and accessions."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            tmp_path = f.name
        try:
            fetch.generate_samplesheet(MOCK_RECORDS, tmp_path)
            with open(tmp_path, newline='') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            assert len(rows) == 3
            assert rows[0]['sample'] == 'SRR12345678'
            assert rows[1]['sample'] == 'SRR12345679'
            assert rows[2]['sample'] == 'SRR12345680'
            assert 'fastq_1' in rows[0]
            assert 'fastq_2' in rows[0]
        finally:
            os.unlink(tmp_path)

    def test_generate_samplesheet_header(self):
        """Samplesheet has the expected header."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            tmp_path = f.name
        try:
            fetch.generate_samplesheet(MOCK_RECORDS, tmp_path)
            with open(tmp_path) as f:
                header = f.readline().strip()
            assert header == 'sample,fastq_1,fastq_2'
        finally:
            os.unlink(tmp_path)

    def test_generate_samplesheet_empty(self):
        """Empty records produces header-only samplesheet."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            tmp_path = f.name
        try:
            fetch.generate_samplesheet([], tmp_path)
            with open(tmp_path, newline='') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            assert len(rows) == 0
        finally:
            os.unlink(tmp_path)


class TestGeoLocParsing:
    """Test parsing of geo_loc_name field."""

    def test_country_region(self):
        """Standard 'Country: Region' format."""
        country, region = fetch.parse_geo_loc_name('Saudi Arabia: Riyadh')
        assert country == 'Saudi Arabia'
        assert region == 'Riyadh'

    def test_country_only(self):
        """Country without region."""
        country, region = fetch.parse_geo_loc_name('Germany')
        assert country == 'Germany'
        assert region == ''

    def test_country_region_spaces(self):
        """Extra whitespace is trimmed."""
        country, region = fetch.parse_geo_loc_name('  USA :  California  ')
        assert country == 'USA'
        assert region == 'California'

    def test_empty_string(self):
        """Empty string returns empty tuple."""
        country, region = fetch.parse_geo_loc_name('')
        assert country == ''
        assert region == ''

    def test_none_value(self):
        """None value returns empty tuple."""
        country, region = fetch.parse_geo_loc_name(None)
        assert country == ''
        assert region == ''

    def test_multiple_colons(self):
        """Multiple colons: only split on first."""
        country, region = fetch.parse_geo_loc_name('USA: California: Los Angeles')
        assert country == 'USA'
        assert region == 'California: Los Angeles'


class TestBuildSearchQuery:
    """Test NCBI search query construction."""

    def test_basic_query(self):
        """Basic query with organism only."""
        q = fetch.build_search_query('Staphylococcus aureus')
        assert '"Staphylococcus aureus"[Organism]' in q
        assert 'paired[Layout]' in q
        assert 'illumina[Platform]' in q

    def test_query_with_country(self):
        """Query with country filter."""
        q = fetch.build_search_query('Staphylococcus aureus', country='Saudi Arabia')
        assert '"Saudi Arabia"[geo_loc_name]' in q

    def test_query_with_date_range(self):
        """Query with date range."""
        q = fetch.build_search_query('Staphylococcus aureus', date_range='2020-2024')
        assert '2020/01/01' in q
        assert '2024/12/31' in q

    def test_query_with_bioproject(self):
        """Query with BioProject."""
        q = fetch.build_search_query('Staphylococcus aureus', bioproject='PRJNA123456')
        assert 'PRJNA123456[BioProject]' in q

    def test_query_all_filters(self):
        """Query with all filters combined."""
        q = fetch.build_search_query(
            'Staphylococcus aureus',
            country='Germany',
            date_range='2021-2023',
            bioproject='PRJNA999',
        )
        assert '"Staphylococcus aureus"[Organism]' in q
        assert '"Germany"[geo_loc_name]' in q
        assert 'PRJNA999[BioProject]' in q
        assert ' AND ' in q


class TestWriteReadRoundtrip:
    """Test TSV write/read roundtrip."""

    def test_roundtrip(self):
        """Records survive a write -> read roundtrip."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            tmp_path = f.name
        try:
            fetch.write_search_results(MOCK_RECORDS, tmp_path)
            loaded = fetch.read_search_results(tmp_path)
            assert len(loaded) == len(MOCK_RECORDS)
            for orig, loaded_rec in zip(MOCK_RECORDS, loaded):
                for col in fetch.OUTPUT_COLUMNS:
                    assert loaded_rec[col] == orig[col], f'Mismatch on {col}'
        finally:
            os.unlink(tmp_path)


class TestGenerateMetadataCSV:
    """Test PHA4GE metadata CSV generation."""

    def test_generate_metadata(self):
        """Mapped records produce a valid CSV with PHA4GE fields."""
        mapped = fetch.map_to_pha4ge(MOCK_RECORDS)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            tmp_path = f.name
        try:
            fetch.generate_metadata_csv(mapped, tmp_path)
            with open(tmp_path, newline='') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            assert len(rows) == 3
            assert rows[0]['sra_accession'] == 'SRR12345678'
            assert rows[0]['geo_loc_country'] == 'Saudi Arabia'
            assert rows[0]['geo_loc_region'] == 'Riyadh'
            assert rows[0]['sample_id'] == 'SRR12345678'
        finally:
            os.unlink(tmp_path)

    def test_generate_metadata_empty(self):
        """Empty records list produces no file content issues."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            tmp_path = f.name
        try:
            fetch.generate_metadata_csv([], tmp_path)
            # File should still exist (generate_metadata_csv returns early for empty)
            assert os.path.exists(tmp_path)
        finally:
            os.unlink(tmp_path)
