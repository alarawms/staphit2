#!/usr/bin/env python3
"""
pubmlst_fetch.py — Fetch public S. aureus WGS isolates by ST, CC, or location.
Optionally downloads the FASTQ files and writes a staphit2-ready samplesheet.

Sources used (in combination):
  1. PubMLST BIGSdb  — precise ST/CC/country typing queries (metadata, limited SRA links)
  2. NCBI Entrez SRA — text-based ST search; returns run accessions + metadata
  3. ENA Portal API  — geographic / host / year filtering of WGS reads

Outputs:
  * metadata TSV     — all isolates with provenance fields
  * samplesheet CSV  — staphit2-ready, with local paths if --download-dir is set

Usage:
  # Fetch metadata only (FTP URLs in samplesheet, no download)
  python pubmlst_fetch.py --st 97 --out st97.tsv --samplesheet st97_samplesheet.csv

  # Fetch metadata AND download FASTQs → ready to pass straight to staphit2
  python pubmlst_fetch.py --st 97,1153 \\
      --download-dir /path/to/reads/st97 \\
      --samplesheet st97_local.csv \\
      --sra-only

  # Limit download to N samples (useful for pilot runs)
  python pubmlst_fetch.py --st 97 --continent europe \\
      --download-dir ./reads --max-downloads 30 --sra-only

  # CC97 bovine Europe context genomes for a phylogeny run
  python pubmlst_fetch.py --st 398 --host "Bos taurus" --continent europe \\
      --download-dir ./reads/bovine_eu --samplesheet bovine_eu.csv --sra-only

After downloading, run staphit2:
  nextflow run /path/to/staphit2 --input st97_local.csv --outdir results/st97

Credentials (PubMLST BIGSdb OAuth1 personal key):
  Consumer key:    buqeJdSSd50NMsH4ixIP4mdz
  Consumer secret: mmDBv0FqiKPuuxzrLdI4Ne1khcwLYMGnLeKWKhyNDV
"""

import argparse
import csv
import io
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import parse_qs, urlencode
from xml.etree import ElementTree as ET

import requests
from oauthlib.oauth1 import Client as OAuth1Client

# ── Credentials ───────────────────────────────────────────────────────────────
CONSUMER_KEY    = "buqeJdSSd50NMsH4ixIP4mdz"
CONSUMER_SECRET = "mmDBv0FqiKPuuxzrLdI4Ne1khcwLYMGnLeKWKhyNDV"

# ── API constants ─────────────────────────────────────────────────────────────
PUBMLST_BASE   = "https://rest.pubmlst.org"
DB_ISOLATES    = "pubmlst_saureus_isolates"
DB_SEQDEF      = "pubmlst_saureus_seqdef"
MLST_SCHEME_ID = 1

ENA_PORTAL     = "https://www.ebi.ac.uk/ena/portal/api/search"
NCBI_ESEARCH   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESUMMARY  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

SESSION_CACHE  = Path.home() / ".pubmlst_session.json"

# Known CC → member STs fallback (if seqdef CC lookup fails)
CC_ST_MAP: dict[str, list[int]] = {
    "97":  [97, 1153, 1154, 2249, 2573, 2648],
    "8":   [8, 247, 250, 254, 612],
    "5":   [5, 225, 231, 1821],
    "22":  [22, 128, 217, 1294],
    "30":  [30, 34, 36, 45, 1374],
    "45":  [45, 574],
    "398": [398, 291, 1232],
    "239": [239, 241, 612],
}

# ── PubMLST OAuth1 helpers ────────────────────────────────────────────────────

def _qsign_get(url: str, consumer_key: str, consumer_secret: str,
               token: str = "", token_secret: str = "",
               verifier: str = "", callback: str = None) -> requests.Response:
    """
    OAuth1 signed GET using QUERY signature type (parameters in URL).
    PubMLST requires QUERY mode; HEADER mode is rejected with 403.
    """
    cli = OAuth1Client(
        consumer_key, client_secret=consumer_secret,
        resource_owner_key=token or None,
        resource_owner_secret=token_secret or None,
        verifier=verifier or None,
        callback_uri=callback,
        signature_type="QUERY",
    )
    uri, headers, _ = cli.sign(url, "GET")
    hdrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in headers.items()}
    r = requests.get(uri, headers=hdrs, timeout=30)
    return r


def _load_session() -> tuple[str, str] | None:
    if not SESSION_CACHE.exists():
        return None
    try:
        data = json.loads(SESSION_CACHE.read_text())
        if data.get("expires_at", 0) > time.time() + 60:
            return data["token"], data["secret"]
    except Exception:
        pass
    return None


def _save_session(token: str, secret: str) -> None:
    SESSION_CACHE.write_text(json.dumps({
        "token": token, "secret": secret,
        "expires_at": time.time() + 11 * 3600,
    }))
    SESSION_CACHE.chmod(0o600)


def get_pubmlst_session(db: str = DB_ISOLATES) -> tuple[str, str]:
    """
    Obtain a BIGSdb session token via OAuth1 personal-key flow.
    Personal keys return oauth_verifier inline in get_request_token,
    bypassing the browser redirect.
    """
    cached = _load_session()
    if cached:
        return cached

    print("Authenticating with PubMLST...", file=sys.stderr)

    # Step 1: request token (callback=oob for personal/pre-authorized key)
    url_rt = f"{PUBMLST_BASE}/db/{db}/oauth/get_request_token"
    r1 = _qsign_get(url_rt, CONSUMER_KEY, CONSUMER_SECRET, callback="oob")
    if r1.status_code != 200:
        raise RuntimeError(f"Request token failed {r1.status_code}: {r1.text}")
    d1 = r1.json()
    req_tok    = d1["oauth_token"]
    req_secret = d1["oauth_token_secret"]
    verifier   = d1.get("oauth_verifier")

    if not verifier:
        raise RuntimeError(
            "No oauth_verifier returned — this key requires browser authorization.\n"
            "Visit: https://pubmlst.org/bigsdb?db=pubmlst_saureus_isolates&page=authorizeClient\n"
            "Then re-run with --no-cache after completing the authorization."
        )

    # Step 2: access token
    url_at = f"{PUBMLST_BASE}/db/{db}/oauth/get_access_token"
    r2 = _qsign_get(url_at, CONSUMER_KEY, CONSUMER_SECRET,
                    token=req_tok, token_secret=req_secret, verifier=verifier)
    if r2.status_code != 200:
        raise RuntimeError(f"Access token failed {r2.status_code}: {r2.text}")
    d2 = r2.json()
    acc_tok    = d2["oauth_token"]
    acc_secret = d2["oauth_token_secret"]

    # Step 3: session token (12-hour validity)
    url_st = f"{PUBMLST_BASE}/db/{db}/oauth/get_session_token"
    r3 = _qsign_get(url_st, CONSUMER_KEY, CONSUMER_SECRET,
                    token=acc_tok, token_secret=acc_secret)
    if r3.status_code != 200:
        raise RuntimeError(f"Session token failed {r3.status_code}: {r3.text}")
    d3 = r3.json()
    sess_tok    = d3["oauth_token"]
    sess_secret = d3["oauth_token_secret"]

    _save_session(sess_tok, sess_secret)
    print("  PubMLST authentication OK.", file=sys.stderr)
    return sess_tok, sess_secret


def pubmlst_post(url: str, body: dict,
                 session: tuple[str, str] | None = None) -> dict:
    """POST JSON to PubMLST. Auth optional — many endpoints are public."""
    headers = {"Content-Type": "application/json"}
    if session:
        cli = OAuth1Client(CONSUMER_KEY, CONSUMER_SECRET, session[0], session[1],
                           signature_type="QUERY")
        url, _, _ = cli.sign(url, "GET")  # QUERY adds token params
    r = requests.post(url, json=body, headers=headers, timeout=30)
    r.raise_for_status()
    return r.json()


def pubmlst_get(url: str, session: tuple[str, str] | None = None) -> dict:
    """GET from PubMLST. Auth optional."""
    if session:
        cli = OAuth1Client(CONSUMER_KEY, CONSUMER_SECRET, session[0], session[1],
                           signature_type="QUERY")
        url, headers, _ = cli.sign(url, "GET")
        headers = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in headers.items()}
    else:
        headers = {}
    r = requests.get(url, headers=headers, timeout=30)
    r.raise_for_status()
    return r.json()


# ── PubMLST source ────────────────────────────────────────────────────────────

def pubmlst_get_cc_sts(cc: str, session) -> list[int]:
    """Return all STs in a clonal complex from seqdef, or fallback map."""
    # Try to look up from seqdef scheme profiles (may have clonal_complex field)
    try:
        url = f"{PUBMLST_BASE}/db/{DB_SEQDEF}/schemes/{MLST_SCHEME_ID}/profiles_csv"
        r = requests.get(url, timeout=60)
        if r.status_code == 200:
            reader = csv.DictReader(io.StringIO(r.text), delimiter="\t")
            sts = [int(row["ST"]) for row in reader
                   if row.get("clonal_complex", "").strip() == str(cc)]
            if sts:
                print(f"  CC{cc}: {len(sts)} member STs from seqdef.", file=sys.stderr)
                return sts
    except Exception as e:
        print(f"  Seqdef CC lookup failed ({e}), using fallback.", file=sys.stderr)

    sts = CC_ST_MAP.get(str(cc))
    if not sts:
        sys.exit(f"ERROR: CC{cc} not in fallback map. Use --st with explicit STs.")
    print(f"  CC{cc} fallback STs: {sts}", file=sys.stderr)
    return sts


def _fetch_pubmlst_isolate(uri: str) -> dict:
    """Fetch a single PubMLST isolate record."""
    try:
        r = requests.get(uri, timeout=20)
        if r.status_code != 200:
            return {}
        data = r.json()
        prov = data.get("provenance", {})

        # Extract ST from schemes block
        st = ""
        for sch in data.get("schemes", []):
            if "ST" in sch.get("fields", {}):
                st = str(sch["fields"]["ST"])
                break

        row = {
            "source":              "pubmlst",
            "isolate_name":        prov.get("isolate", ""),
            "species":             prov.get("species", "Staphylococcus aureus"),
            "ST":                  st,
            "country":             prov.get("country", ""),
            "continent":           prov.get("continent", ""),
            "year":                str(prov.get("year", "")),
            "disease":             prov.get("disease", ""),
            "source_category":     prov.get("source_category", ""),
            "host":                prov.get("host", ""),
            "spa_type":            prov.get("spa_type", ""),
            "run_accession":       prov.get("run_accession", "") or "",
            "biosample_accession": prov.get("biosample_accession", "") or "",
            "assembly_accession":  prov.get("NCBI_assembly_accession", "") or "",
            "pubmlst_id":          str(prov.get("id", "")),
        }
        _add_ena_ftp(row)
        return row
    except Exception:
        return {}


def fetch_from_pubmlst(sts: list[int],
                       country: str | None,
                       continent: str | None,
                       host: str | None,
                       year_from: int | None,
                       year_to: int | None,
                       workers: int = 8) -> list[dict]:
    """Query PubMLST isolates by ST + provenance filters."""
    print("\n[PubMLST] Searching...", file=sys.stderr)
    url = f"{PUBMLST_BASE}/db/{DB_ISOLATES}/isolates/search"

    all_uris: list[str] = []
    for st in sts:
        body: dict = {}
        if st:
            body[f"scheme.{MLST_SCHEME_ID}.ST"] = st
        if country:
            body["field.country"] = country
        if continent:
            body["field.continent"] = continent.lower()
        if host:
            body["field.host"] = host

        # return_all must be a query param, not in body
        try:
            r = requests.post(
                url + "?return_all=1",
                json=body,
                headers={"Content-Type": "application/json"},
                timeout=30,
            )
            if r.status_code != 200:
                print(f"  ST{st} search error {r.status_code}: {r.text[:100]}",
                      file=sys.stderr)
                continue
            data = r.json()
            uris = data.get("isolates", [])
            print(f"  ST{st}: {data.get('records', len(uris))} isolates", file=sys.stderr)
            all_uris.extend(uris)
        except Exception as e:
            print(f"  ST{st} search failed: {e}", file=sys.stderr)

    # Deduplicate
    seen: set[str] = set()
    uris_uniq = [u for u in all_uris if not (u in seen or seen.add(u))]
    print(f"  Total unique PubMLST isolates: {len(uris_uniq)}", file=sys.stderr)

    if not uris_uniq:
        return []

    print(f"  Fetching isolate details ({workers} workers)...", file=sys.stderr)
    results: list[dict] = []
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_fetch_pubmlst_isolate, u): u for u in uris_uniq}
        for fut in as_completed(futs):
            row = fut.result()
            if row:
                # Apply year filter (PubMLST doesn't support year range in POST)
                y = row.get("year", "")
                if y and y.isdigit():
                    yi = int(y)
                    if year_from and yi < year_from:
                        continue
                    if year_to and yi > year_to:
                        continue
                results.append(row)
            done += 1
            if done % 100 == 0:
                print(f"  {done}/{len(uris_uniq)} fetched...", file=sys.stderr)

    return results


# ── NCBI Entrez SRA source ────────────────────────────────────────────────────

def _build_ncbi_query(sts: list[int],
                      country: str | None,
                      continent: str | None,
                      host: str | None,
                      year_from: int | None,
                      year_to: int | None) -> str:
    """Build an Entrez search query for SRA."""
    terms = ['Staphylococcus aureus[Organism]', 'WGS[Strategy]']

    # ST terms: search for "ST97" OR "ST 97" OR "sequence type 97" OR "CC97"
    st_terms = []
    for st in sts:
        st_terms += [f'"ST{st}"', f'"ST {st}"', f'"sequence type {st}"']
    if st_terms:
        terms.append(f'({" OR ".join(st_terms)})')

    if country:
        terms.append(f'"{country}"[All Fields]')
    if continent:
        continent_countries = {
            "europe": ["UK", "Germany", "France", "Spain", "Belgium", "Netherlands",
                       "Italy", "Sweden", "Denmark", "Norway", "Switzerland"],
            "asia": ["Saudi Arabia", "Iran", "India", "China", "Japan", "Korea"],
            "americas": ["USA", "Canada", "Brazil"],
            "africa": ["Egypt", "South Africa", "Nigeria"],
        }
        if continent.lower() in continent_countries:
            countries = continent_countries[continent.lower()]
            terms.append(f'({" OR ".join(f"{c}[All Fields]" for c in countries)})')
    if host:
        safe_host = host.replace('"', '')
        terms.append(f'"{safe_host}"[All Fields]')
    if year_from or year_to:
        yf = year_from or 1990
        yt = year_to or 2030
        terms.append(f'{yf}:{yt}[Publication Date]')

    return " AND ".join(terms)


def _parse_sra_efetch(xml_text: str) -> list[dict]:
    """Parse SRA efetch XML into isolate rows."""
    rows = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return rows

    for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
        run_el = pkg.find(".//RUN")
        if run_el is None:
            continue
        run_acc = run_el.get("accession", "")

        # Sample attributes
        attrs: dict[str, str] = {}
        for attr in pkg.findall(".//SAMPLE_ATTRIBUTE"):
            tag = (attr.findtext("TAG") or "").lower().strip()
            val = (attr.findtext("VALUE") or "").strip()
            attrs[tag] = val

        # Resolve ST from attributes
        st = ""
        for key in ("mlst st", "sequence type", "st", "mlst_st", "sequence_type",
                    "mlst sequence type"):
            if key in attrs:
                st = attrs[key]
                break

        # Infer ST from run title / description if not in attrs
        sample = pkg.find(".//SAMPLE")
        title = ""
        if sample is not None:
            title = (sample.findtext(".//TITLE") or
                     sample.findtext(".//NAME") or "")
        if not st:
            m = re.search(r'\bST[-\s]?(\d+)\b', title, re.I)
            if m:
                st = m.group(1)

        row = {
            "source":              "ncbi_sra",
            "isolate_name":        attrs.get("isolate", attrs.get("strain", title)),
            "species":             "Staphylococcus aureus",
            "ST":                  st,
            "country":             attrs.get("geo_loc_name", attrs.get("country", "")),
            "continent":           "",
            "year":                attrs.get("collection_date", "")[:4],
            "disease":             attrs.get("disease", attrs.get("host_disease", "")),
            "source_category":     attrs.get("isolation_source", ""),
            "host":                attrs.get("host", ""),
            "spa_type":            attrs.get("spa_type", attrs.get("spa type", "")),
            "run_accession":       run_acc,
            "biosample_accession": pkg.findtext(".//IDENTIFIERS/PRIMARY_ID") or "",
            "assembly_accession":  "",
            "pubmlst_id":          "",
        }
        _add_ena_ftp(row)
        rows.append(row)
    return rows


def fetch_from_ncbi(sts: list[int],
                    country: str | None,
                    continent: str | None,
                    host: str | None,
                    year_from: int | None,
                    year_to: int | None,
                    ncbi_email: str = "alarawms@gmail.com") -> list[dict]:
    """Search NCBI SRA for WGS runs matching ST/location criteria."""
    print("\n[NCBI SRA] Searching...", file=sys.stderr)

    query = _build_ncbi_query(sts, country, continent, host, year_from, year_to)
    print(f"  Query: {query[:120]}...", file=sys.stderr)

    # esearch
    r = requests.get(NCBI_ESEARCH, params={
        "db": "sra", "term": query,
        "retmax": 10000, "retmode": "json",
        "email": ncbi_email,
    }, timeout=30)
    r.raise_for_status()
    result = r.json().get("esearchresult", {})
    ids = result.get("idlist", [])
    total = result.get("count", "?")
    print(f"  {total} SRA records found.", file=sys.stderr)

    if not ids:
        return []

    # Fetch in batches of 200
    all_rows: list[dict] = []
    batch_size = 200
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        r2 = requests.post(NCBI_EFETCH, data={
            "db": "sra", "id": ",".join(batch),
            "rettype": "full", "retmode": "xml",
            "email": ncbi_email,
        }, timeout=60)
        if r2.status_code == 200:
            rows = _parse_sra_efetch(r2.text)
            all_rows.extend(rows)
            print(f"  Fetched batch {i//batch_size + 1}: {len(rows)} records",
                  file=sys.stderr)
        time.sleep(0.4)   # NCBI rate limit: max 3 req/s without API key

    return all_rows


# ── ENA Portal API source ─────────────────────────────────────────────────────

ENA_FIELDS = (
    "run_accession,sample_accession,study_accession,"
    "country,collection_date,host,strain,sample_description"
)


def fetch_from_ena(sts: list[int],
                   country: str | None,
                   continent: str | None,
                   host: str | None,
                   year_from: int | None,
                   year_to: int | None) -> list[dict]:
    """
    Query ENA Portal API for S. aureus WGS reads with provenance filters.
    ENA does not expose MLST/ST as a searchable field; ST is inferred from
    strain/description text when available.
    """
    print("\n[ENA] Searching...", file=sys.stderr)

    query_parts = ['tax_eq(1280)']
    if country:
        query_parts.append(f'country="{country}"')
    if host:
        safe = host.replace('"', '')
        query_parts.append(f'host="{safe}"')

    query = " AND ".join(query_parts)

    params = {
        "result": "read_run",
        "query": query,
        "fields": ENA_FIELDS,
        "format": "json",
        "limit": 100000,
    }

    r = requests.get(ENA_PORTAL, params=params, timeout=60)
    if r.status_code != 200:
        print(f"  ENA search failed {r.status_code}: {r.text[:200]}", file=sys.stderr)
        return []

    records = r.json()
    print(f"  {len(records)} ENA records returned (pre-ST filter).", file=sys.stderr)

    # Build ST pattern for text matching
    st_patterns = [re.compile(rf'\bST[-\s]?{st}\b', re.I) for st in sts]

    rows: list[dict] = []
    for rec in records:
        # ST filtering: check strain + description fields
        text = " ".join(filter(None, [rec.get("strain", ""), rec.get("sample_description", "")]))
        if sts and not any(p.search(text) for p in st_patterns):
            continue

        # Extract year from collection_date
        year = ""
        cd = rec.get("collection_date", "")
        m = re.search(r'\b(19|20)\d{2}\b', cd)
        if m:
            year = m.group(0)

        # Year range filter
        if year and year.isdigit():
            yi = int(year)
            if year_from and yi < year_from:
                continue
            if year_to and yi > year_to:
                continue

        # Infer ST from text
        st_found = ""
        for pat in st_patterns:
            if pat.search(text):
                st_found = str(sts[st_patterns.index(pat)])
                break

        run_acc = rec.get("run_accession", "")
        row = {
            "source":              "ena",
            "isolate_name":        rec.get("strain", ""),
            "species":             "Staphylococcus aureus",
            "ST":                  st_found,
            "country":             rec.get("country", ""),
            "continent":           "",
            "year":                year,
            "disease":             "",
            "source_category":     "",
            "host":                rec.get("host", ""),
            "spa_type":            "",
            "run_accession":       run_acc,
            "biosample_accession": rec.get("sample_accession", ""),
            "assembly_accession":  "",
            "pubmlst_id":          "",
        }
        _add_ena_ftp(row)
        rows.append(row)

    print(f"  {len(rows)} records after ST filter.", file=sys.stderr)
    return rows


# ── Shared utilities ──────────────────────────────────────────────────────────

def _add_ena_ftp(row: dict) -> None:
    """Derive ENA FTP download URLs from run_accession."""
    run = row.get("run_accession", "")
    row["ena_ftp_1"] = ""
    row["ena_ftp_2"] = ""
    if run and re.match(r'^(ERR|SRR|DRR)\d+$', run):
        prefix = run[:6]
        sub    = run[:9]
        base   = f"ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{sub}/{run}"
        row["ena_ftp_1"] = f"{base}/{run}_1.fastq.gz"
        row["ena_ftp_2"] = f"{base}/{run}_2.fastq.gz"


def _dedup(rows: list[dict]) -> list[dict]:
    """Deduplicate by run_accession (prefer ncbi_sra > pubmlst > ena)."""
    source_rank = {"ncbi_sra": 0, "pubmlst": 1, "ena": 2}
    by_run: dict[str, dict] = {}
    no_run: list[dict] = []

    for row in rows:
        run = row.get("run_accession", "")
        if not run:
            no_run.append(row)
            continue
        existing = by_run.get(run)
        if existing is None or (source_rank.get(row["source"], 9) <
                                source_rank.get(existing["source"], 9)):
            by_run[run] = row

    result = list(by_run.values()) + no_run
    result.sort(key=lambda r: (r.get("country", ""), r.get("ST", ""),
                                r.get("year", "")))
    return result


# ── Output ────────────────────────────────────────────────────────────────────

OUTPUT_COLS = [
    "source", "isolate_name", "ST", "country", "continent", "year",
    "disease", "source_category", "host", "spa_type",
    "run_accession", "biosample_accession", "assembly_accession",
    "pubmlst_id", "ena_ftp_1", "ena_ftp_2",
]


def write_metadata(rows: list[dict], path: str, sra_only: bool) -> int:
    if sra_only:
        rows = [r for r in rows if r.get("run_accession")]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLS,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_samplesheet(rows: list[dict], path: str,
                      download_dir: Path | None = None) -> int:
    """
    Write samplesheet CSV.  If download_dir is set, uses local file paths;
    otherwise uses ENA FTP URLs (not valid for staphit2 without downloading).
    """
    runnable = [r for r in rows
                if r.get("run_accession") and r.get("ena_ftp_1")]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        for r in runnable:
            run = r["run_accession"]
            if download_dir:
                f1 = str(download_dir / f"{run}_1.fastq.gz")
                f2 = str(download_dir / f"{run}_2.fastq.gz")
            else:
                f1 = "https://" + r["ena_ftp_1"].removeprefix("ftp://").removeprefix("https://")
                f2 = "https://" + r["ena_ftp_2"].removeprefix("ftp://").removeprefix("https://") if r.get("ena_ftp_2") else ""
            writer.writerow([run, f1, f2])
    return len(runnable)


# ── Download ──────────────────────────────────────────────────────────────────

def _download_one(run: str, ftp_1: str, ftp_2: str,
                  dest: Path) -> tuple[str, bool, str]:
    """
    Download paired FASTQ files for one run via ENA FTP.
    Returns (run_accession, success, message).
    Skips files that already exist with non-zero size.
    """
    import subprocess
    dest.mkdir(parents=True, exist_ok=True)
    errors = []
    for ftp_url, suffix in [(ftp_1, "_1.fastq.gz"), (ftp_2, "_2.fastq.gz")]:
        if not ftp_url:
            errors.append(f"{suffix}: no FTP URL")
            continue
        local = dest / f"{run}{suffix}"
        if local.exists() and local.stat().st_size > 0:
            continue
        local.unlink(missing_ok=True)  # remove any zero-byte partial from a prior attempt
        bare = ftp_url.removeprefix("ftp://").removeprefix("https://").removeprefix("http://")
        url  = f"https://{bare}"
        try:
            result = subprocess.run(
                ["wget", "--tries=3", "--timeout=60", "-O", str(local), url],
                capture_output=True, timeout=3600,
            )
            stderr_lines = result.stderr.decode(errors="replace").strip().splitlines()
            stderr_tail  = " | ".join(stderr_lines[-3:]) if stderr_lines else ""
            if result.returncode != 0 or not local.exists() or local.stat().st_size == 0:
                local.unlink(missing_ok=True)
                errors.append(
                    f"{suffix}: wget rc={result.returncode}"
                    + (f" — {stderr_tail}" if stderr_tail else "")
                )
        except subprocess.TimeoutExpired:
            local.unlink(missing_ok=True)
            errors.append(f"{suffix}: timeout")
        except Exception as e:
            errors.append(f"{suffix}: {e}")

    if errors:
        return run, False, "; ".join(errors)
    return run, True, "ok"


def download_reads(rows: list[dict], download_dir: Path,
                   max_downloads: int | None = None,
                   workers: int = 4) -> list[dict]:
    """
    Download FASTQ files for all rows that have ena_ftp_1.
    Updates rows in-place with local_fastq_1/2 fields.
    Returns only rows that downloaded successfully.
    """
    import shutil
    if not shutil.which("wget"):
        sys.exit("ERROR: wget not found. Install with: sudo dnf install wget")

    download_dir.mkdir(parents=True, exist_ok=True)

    targets = [r for r in rows if r.get("run_accession") and r.get("ena_ftp_1")]
    if max_downloads:
        targets = targets[:max_downloads]

    total = len(targets)
    print(f"\nDownloading {total} samples → {download_dir}/", file=sys.stderr)
    print(f"  ({workers} parallel streams; skips existing files)\n", file=sys.stderr)

    success_runs: set[str] = set()
    done = 0

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {
            pool.submit(_download_one,
                        r["run_accession"], r.get("ena_ftp_1"), r.get("ena_ftp_2"),
                        download_dir): r
            for r in targets
        }
        for fut in as_completed(futs):
            run, ok, msg = fut.result()
            done += 1
            if ok:
                success_runs.add(run)
                print(f"  [{done}/{total}] ✓ {run}", file=sys.stderr)
            else:
                print(f"  [{done}/{total}] ✗ {run}: {msg}", file=sys.stderr)

    print(f"\n  Downloaded: {len(success_runs)}/{total}", file=sys.stderr)
    return [r for r in targets if r["run_accession"] in success_runs]


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Fetch public S. aureus WGS isolates by ST / CC / location.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Credentials")[0].strip(),
    )
    p.add_argument("--st",  metavar="ST[,ST2,...]",
                   help="Sequence type(s), comma-separated (e.g. 97,1153)")
    p.add_argument("--cc",  metavar="CC",
                   help="Clonal complex (e.g. 97). Expands to member STs.")
    p.add_argument("--country",    help="Country filter (e.g. 'Saudi Arabia')")
    p.add_argument("--continent",  help="Continent: europe / asia / africa / americas")
    p.add_argument("--host",       help="Host species (e.g. 'Homo sapiens', 'Bos taurus')")
    p.add_argument("--year-from",  type=int, metavar="YEAR")
    p.add_argument("--year-to",    type=int, metavar="YEAR")
    p.add_argument("--source",     default="all",
                   choices=["all", "pubmlst", "ncbi", "ena"],
                   help="Data source to query (default: all)")
    p.add_argument("--out",        default="pubmlst_metadata.tsv",
                   help="Output metadata TSV (default: pubmlst_metadata.tsv)")
    p.add_argument("--samplesheet", metavar="FILE",
                   help="Write staphit2 samplesheet CSV for SRA-linked isolates")
    p.add_argument("--accessions", metavar="FILE",
                   help="Write one SRA run accession per line (input for nf-core/fetchngs)")
    p.add_argument("--sra-only",   action="store_true",
                   help="Only output isolates with a run_accession")
    p.add_argument("--download-dir", metavar="DIR",
                   help="Download FASTQ files here and use local paths in the samplesheet. "
                        "Files already present are skipped. Requires wget.")
    p.add_argument("--max-downloads", type=int, metavar="N",
                   help="Cap the number of samples downloaded (for pilot runs)")
    p.add_argument("--download-workers", type=int, default=4,
                   help="Parallel download streams (default: 4; be polite to ENA FTP)")
    p.add_argument("--workers",    type=int, default=8,
                   help="Parallel metadata fetch workers (default: 8)")
    p.add_argument("--no-cache",   action="store_true",
                   help="Discard cached PubMLST session token")
    return p.parse_args()


def main():
    args = parse_args()

    if not args.st and not args.cc:
        sys.exit("ERROR: specify at least --st or --cc")

    if args.no_cache and SESSION_CACHE.exists():
        SESSION_CACHE.unlink()

    # Resolve STs
    sts: list[int] = []
    if args.cc:
        # PubMLST session needed for live CC lookup
        try:
            session = get_pubmlst_session()
        except Exception as e:
            print(f"  PubMLST auth failed ({e}), using CC fallback map.", file=sys.stderr)
            session = None
        sts.extend(pubmlst_get_cc_sts(args.cc, session))

    if args.st:
        for s in args.st.split(","):
            v = int(s.strip())
            if v not in sts:
                sts.append(v)

    print(f"Target STs: {sts}", file=sys.stderr)
    print(f"Filters — country: {args.country or 'any'}, "
          f"continent: {args.continent or 'any'}, "
          f"host: {args.host or 'any'}", file=sys.stderr)

    all_rows: list[dict] = []

    # Query selected sources
    use_pubmlst = args.source in ("all", "pubmlst")
    use_ncbi    = args.source in ("all", "ncbi")
    use_ena     = args.source in ("all", "ena")

    if use_pubmlst:
        try:
            session = get_pubmlst_session()
        except Exception as e:
            print(f"  PubMLST auth unavailable ({e}), skipping.", file=sys.stderr)
            session = None
        rows = fetch_from_pubmlst(
            sts, args.country, args.continent, args.host,
            args.year_from, args.year_to, workers=args.workers,
        )
        all_rows.extend(rows)

    if use_ncbi:
        rows = fetch_from_ncbi(
            sts, args.country, args.continent, args.host,
            args.year_from, args.year_to,
        )
        all_rows.extend(rows)

    if use_ena:
        rows = fetch_from_ena(
            sts, args.country, args.continent, args.host,
            args.year_from, args.year_to,
        )
        all_rows.extend(rows)

    # Merge + deduplicate
    merged = _dedup(all_rows)
    print(f"\nTotal unique records: {len(merged)}", file=sys.stderr)
    sra_n = sum(1 for r in merged if r.get("run_accession"))
    print(f"With SRA run_accession: {sra_n}", file=sys.stderr)

    if not merged:
        print("No records found.", file=sys.stderr)
        sys.exit(0)

    # Write outputs
    n = write_metadata(merged, args.out, sra_only=args.sra_only)
    print(f"Wrote {n} rows → {args.out}", file=sys.stderr)

    download_dir: Path | None = None
    if args.download_dir:
        download_dir = Path(args.download_dir)
        merged = download_reads(
            merged, download_dir,
            max_downloads=args.max_downloads,
            workers=args.download_workers,
        )

    if args.accessions:
        acc_rows = [r for r in merged if r.get("run_accession")]
        with open(args.accessions, "w") as fh:
            for r in acc_rows:
                fh.write(r["run_accession"] + "\n")
        print(f"Wrote {len(acc_rows)} accessions → {args.accessions}", file=sys.stderr)

    if args.samplesheet:
        n2 = write_samplesheet(merged, args.samplesheet, download_dir=download_dir)
        print(f"Wrote {n2} SRA entries → {args.samplesheet}", file=sys.stderr)

    # Summary
    from collections import Counter
    print("\n── Summary ──────────────────────────────────────────", file=sys.stderr)
    print(f"By source:   {dict(Counter(r['source'] for r in merged).most_common())}",
          file=sys.stderr)
    print(f"By ST:       {dict(Counter(r['ST'] for r in merged if r['ST']).most_common(8))}",
          file=sys.stderr)
    print(f"By country:  {dict(Counter(r['country'] for r in merged if r['country']).most_common(5))}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
