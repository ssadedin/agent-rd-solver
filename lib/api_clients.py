"""
Thin, cached clients for the external knowledge sources the workflow is allowed to use:
ClinVar + PubMed (NCBI E-utilities), Ensembl REST, UCSC REST, gnomAD GraphQL.

Design goals:
- Every response is cached to disk (audit/api_cache/) so a run is reproducible and the raw
  payload is available as provenance for an evidence_ref.
- Functions return parsed primitives AND record the accession/URL used, so the agent can cite it.
- Rate limiting is respected (NCBI: 10 rps with key, 3 without).

These are intentionally dependency-light (requests + stdlib). Network failures raise; the
orchestrator decides whether to retry or mark the fact UNKNOWN.
"""
from __future__ import annotations
import hashlib
import json
import os
import time
import urllib.parse
from pathlib import Path

import requests

CACHE_DIR = Path(os.environ.get("WF_API_CACHE", "audit/api_cache"))
CACHE_DIR.mkdir(parents=True, exist_ok=True)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_BASE = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
UCSC_BASE = "https://api.genome.ucsc.edu"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"

_NCBI_KEY = os.environ.get("NCBI_API_KEY")
_last_ncbi_call = [0.0]
_NCBI_MIN_INTERVAL = 0.11 if _NCBI_KEY else 0.34


def _cache_path(tag: str, key: str) -> Path:
    h = hashlib.sha256(key.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{tag}_{h}.json"


def _cached_get(tag: str, url: str, params: dict | None = None, ncbi: bool = False) -> dict:
    key = url + "?" + urllib.parse.urlencode(sorted((params or {}).items()))
    cp = _cache_path(tag, key)
    if cp.exists():
        return json.loads(cp.read_text())
    if ncbi:
        delta = time.time() - _last_ncbi_call[0]
        if delta < _NCBI_MIN_INTERVAL:
            time.sleep(_NCBI_MIN_INTERVAL - delta)
        _last_ncbi_call[0] = time.time()
    r = requests.get(url, params=params, timeout=30,
                     headers={"Content-Type": "application/json",
                              "User-Agent": "rare-disease-workflow/1.0"})
    r.raise_for_status()
    try:
        payload = r.json()
    except ValueError:
        payload = {"_text": r.text}
    payload["_provenance"] = {"url": r.url, "retrieved_at": _now()}
    cp.write_text(json.dumps(payload, indent=2))
    return payload


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# ----------------------------------------------------------------------------- ClinVar
def clinvar_by_variant(hgvs_or_rsid: str) -> dict:
    """
    Look up ClinVar for a variant (HGVS, rsID, or 'chr-pos-ref-alt' spdi-ish).
    Returns {variation_id, clinical_significance, review_status, star_rating,
             conditions[], last_evaluated, accession, provenance}.
    """
    params = {"db": "clinvar", "term": hgvs_or_rsid, "retmode": "json"}
    if _NCBI_KEY:
        params["api_key"] = _NCBI_KEY
    es = _cached_get("clinvar_search", f"{NCBI_BASE}/esearch.fcgi", params, ncbi=True)
    ids = es.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"variation_id": None, "clinical_significance": "not_in_clinvar",
                "_provenance": es.get("_provenance")}
    sp = {"db": "clinvar", "id": ids[0], "retmode": "json"}
    if _NCBI_KEY:
        sp["api_key"] = _NCBI_KEY
    summ = _cached_get("clinvar_summary", f"{NCBI_BASE}/esummary.fcgi", sp, ncbi=True)
    doc = summ.get("result", {}).get(ids[0], {})
    germ = doc.get("germline_classification", {})
    return {
        "variation_id": doc.get("accession") or ids[0],
        "clinical_significance": germ.get("description", doc.get("clinical_significance", {}).get("description", "UNKNOWN")),
        "review_status": germ.get("review_status", "UNKNOWN"),
        "star_rating": _stars(germ.get("review_status", "")),
        "conditions": [t.get("trait_name") for t in germ.get("trait_set", []) if t.get("trait_name")],
        "last_evaluated": germ.get("last_evaluated", ""),
        "_provenance": summ.get("_provenance"),
    }


def _stars(review_status: str) -> int:
    rs = (review_status or "").lower()
    if "practice guideline" in rs:
        return 4
    if "reviewed by expert panel" in rs:
        return 3
    if "multiple submitters" in rs and "no conflict" in rs:
        return 2
    if "single submitter" in rs or "conflicting" in rs:
        return 1
    return 0


def clinvar_same_residue(gene: str, protein_change: str) -> list[dict]:
    """Find ClinVar variants at the same protein residue (for PM5/PS1). Returns summaries."""
    term = f"{gene}[gene] AND {protein_change}"
    params = {"db": "clinvar", "term": term, "retmode": "json", "retmax": "50"}
    if _NCBI_KEY:
        params["api_key"] = _NCBI_KEY
    es = _cached_get("clinvar_residue", f"{NCBI_BASE}/esearch.fcgi", params, ncbi=True)
    return es.get("esearchresult", {}).get("idlist", [])


# ----------------------------------------------------------------------------- PubMed
def pubmed_search(query: str, retmax: int = 20) -> list[str]:
    params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": str(retmax),
              "sort": "relevance"}
    if _NCBI_KEY:
        params["api_key"] = _NCBI_KEY
    es = _cached_get("pubmed_search", f"{NCBI_BASE}/esearch.fcgi", params, ncbi=True)
    return es.get("esearchresult", {}).get("idlist", [])


def pubmed_summaries(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "json"}
    if _NCBI_KEY:
        params["api_key"] = _NCBI_KEY
    s = _cached_get("pubmed_summary", f"{NCBI_BASE}/esummary.fcgi", params, ncbi=True)
    out = []
    res = s.get("result", {})
    for pid in res.get("uids", []):
        d = res[pid]
        out.append({"pmid": pid, "title": d.get("title", ""),
                    "year": (d.get("pubdate", "")[:4] or None),
                    "journal": d.get("fulljournalname", "")})
    return out


# ----------------------------------------------------------------------------- Ensembl
def ensembl_vep(region_allele: str, assembly: str = "GRCh38") -> dict:
    """region_allele like '9:22125503-22125502:1/C' or HGVS via /vep/human/hgvs/."""
    base = ENSEMBL_BASE[assembly]
    return _cached_get("ensembl_vep", f"{base}/vep/human/region/{region_allele}",
                       {"content-type": "application/json"})


def ensembl_lookup_symbol(symbol: str, assembly: str = "GRCh38") -> dict:
    base = ENSEMBL_BASE[assembly]
    return _cached_get("ensembl_symbol", f"{base}/lookup/symbol/homo_sapiens/{symbol}",
                       {"expand": "1"})


def ensembl_phenotype_gene(symbol: str, assembly: str = "GRCh38") -> dict:
    base = ENSEMBL_BASE[assembly]
    return _cached_get("ensembl_pheno", f"{base}/phenotype/gene/homo_sapiens/{symbol}", {})


# ----------------------------------------------------------------------------- UCSC
def ucsc_track(track: str, chrom: str, start: int, end: int, genome: str = "hg38") -> dict:
    """Fetch a UCSC track slice, e.g. track='phyloP100way' or 'rmsk', 'genomicSuperDups'."""
    return _cached_get("ucsc", f"{UCSC_BASE}/getData/track",
                       {"genome": genome, "track": track, "chrom": chrom,
                        "start": str(start), "end": str(end)})


# ----------------------------------------------------------------------------- gnomAD
def gnomad_variant(variant_id: str, dataset: str = "gnomad_r4") -> dict:
    """variant_id format '1-55051215-G-GA' (GRCh38). Returns AF/popmax/hom/hemi."""
    query = """
    query V($id: String!, $ds: DatasetId!) {
      variant(variantId: $id, dataset: $ds) {
        variant_id
        genome { ac an af homozygote_count hemizygote_count
                 populations { id ac an } }
        exome  { ac an af homozygote_count hemizygote_count
                 populations { id ac an } }
      }
    }"""
    key = "gnomad_" + variant_id + "_" + dataset
    cp = _cache_path("gnomad", key)
    if cp.exists():
        return json.loads(cp.read_text())
    r = requests.post(GNOMAD_API, json={"query": query,
                                        "variables": {"id": variant_id, "ds": dataset}},
                      timeout=40, headers={"User-Agent": "rare-disease-workflow/1.0"})
    r.raise_for_status()
    data = r.json()
    data["_provenance"] = {"url": GNOMAD_API, "variant_id": variant_id, "retrieved_at": _now()}
    cp.write_text(json.dumps(data, indent=2))
    return data


def gnomad_popmax(variant_id: str, dataset: str = "gnomad_r4") -> dict:
    """Convenience: compute global AF and popmax AF from a gnomad_variant payload."""
    d = gnomad_variant(variant_id, dataset).get("data", {}).get("variant")
    if not d:
        return {"af_global": None, "af_popmax": None, "popmax_population": None,
                "hom": None, "hemi": None}
    ac = an = hom = hemi = 0
    pop_ac: dict[str, list[int]] = {}
    for src in ("genome", "exome"):
        s = d.get(src)
        if not s:
            continue
        ac += s.get("ac") or 0
        an += s.get("an") or 0
        hom += s.get("homozygote_count") or 0
        hemi += s.get("hemizygote_count") or 0
        for p in s.get("populations") or []:
            pid = p["id"]
            # skip aggregate / sex strata for popmax
            if any(pid.endswith(sfx) for sfx in ("_XX", "_XY")) or pid in ("ALL",):
                continue
            cur = pop_ac.setdefault(pid, [0, 0])
            cur[0] += p.get("ac") or 0
            cur[1] += p.get("an") or 0
    af_global = (ac / an) if an else None
    popmax = max(((a / n, pid) for pid, (a, n) in pop_ac.items() if n >= 2000 and a > 0),
                 default=(None, None))
    return {"af_global": af_global, "af_popmax": popmax[0], "popmax_population": popmax[1],
            "hom": hom, "hemi": hemi}
