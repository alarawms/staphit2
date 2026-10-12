"""Tests for bin/staphit-cgmlst-schema."""
import importlib.machinery
import importlib.util
import io
import urllib.error
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'staphit-cgmlst-schema'
_loader = importlib.machinery.SourceFileLoader('cgmlst_schema', str(BIN))
_spec = importlib.util.spec_from_loader('cgmlst_schema', _loader)
cs = importlib.util.module_from_spec(_spec)
_loader.exec_module(cs)


def test_get_waits_out_rate_limit(monkeypatch):
    calls, sleeps = [], []

    def urlopen(url, timeout):
        calls.append(url)
        if len(calls) < 3:
            raise urllib.error.HTTPError(url, 429, 'Too Many Requests', {'Retry-After': '90'}, None)
        return io.BytesIO(b'>a_1\nACGT\n')

    monkeypatch.setattr(cs.urllib.request, 'urlopen', urlopen)
    monkeypatch.setattr(cs.time, 'sleep', sleeps.append)
    assert cs.get('u') == b'>a_1\nACGT\n'
    assert sleeps == [90, 90]  # Retry-After honoured, not the 5 s / 10 s backoff


def test_get_backs_off_on_network_errors(monkeypatch):
    sleeps = []

    def urlopen(url, timeout):
        raise urllib.error.URLError('down')

    monkeypatch.setattr(cs.urllib.request, 'urlopen', urlopen)
    monkeypatch.setattr(cs.time, 'sleep', sleeps.append)
    try:
        cs.get('u', tries=4)
    except RuntimeError as e:
        assert 'down' in str(e)
    else:
        raise AssertionError('expected RuntimeError')
    assert sleeps == [5, 10, 20]
