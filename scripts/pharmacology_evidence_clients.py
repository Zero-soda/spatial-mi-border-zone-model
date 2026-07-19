#!/usr/bin/env python3
"""Cache-first clients and parsers for public pharmacology resources."""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from therapeutic_prioritization_utils import (
    normalize_action,
)


OPEN_TARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_ENDPOINT = "https://www.ebi.ac.uk/chembl/api/data"
DGIDB_ENDPOINT = "https://dgidb.org/api/graphql"
USER_AGENT = "spatial-mi-therapeutic-prioritization/0.1 (public-data research)"
SOURCE_RELEASES = {
    "Open Targets": "26.06",
    "ChEMBL target search": "ChEMBL 36",
    "ChEMBL mechanism": "ChEMBL 36",
    "DGIdb": "5.0",
}


OPEN_TARGETS_QUERY = """
query TargetPharmacology($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    isEssential
    tractability { label modality value }
    safetyLiabilities { event datasource }
    associatedDiseases(page: {index: 0, size: 100}) {
      count
      rows { score disease { id name } }
    }
    drugAndClinicalCandidates {
      count
      rows {
        id
        maxClinicalStage
        drug {
          id
          name
          maximumClinicalStage
          drugType
          drugWarnings { warningType description }
        }
        diseases { diseaseFromSource disease { id name } }
      }
    }
  }
}
""".strip()


DGIDB_QUERY_TEMPLATE = """
query {
  genes(names: %s) {
    nodes {
      name
      conceptId
      interactions {
        drug { name conceptId }
        interactionScore
        interactionTypes { type directionality }
        sources { sourceDbName }
      }
    }
  }
}
""".strip()


CARDIOVASCULAR_TERMS = (
    "myocardial infarction",
    "coronary artery",
    "coronary heart",
    "ischemic heart",
    "ischaemic heart",
    "heart failure",
    "cardiac fibrosis",
    "cardiomyopathy",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _cache_key(source: str, request_spec: Mapping[str, Any]) -> str:
    digest = hashlib.sha256(
        json.dumps(request_spec, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:20]
    return f"{source.lower().replace(' ', '_')}_{digest}.json"


def cached_json_request(
    *,
    source: str,
    endpoint: str,
    cache_dir: Path,
    payload: Mapping[str, Any] | None = None,
    query: Mapping[str, Any] | None = None,
    refresh: bool = False,
    timeout: int = 90,
    retries: int = 4,
) -> dict[str, Any]:
    """Perform a JSON HTTP request and preserve an auditable response wrapper."""

    method = "POST" if payload is not None else "GET"
    request_spec: dict[str, Any] = {
        "method": method,
        "endpoint": endpoint,
        "payload": payload,
        "query": query,
    }
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / _cache_key(source, request_spec)
    if cache_path.exists() and not refresh:
        return json.loads(cache_path.read_text())

    url = endpoint
    body: bytes | None = None
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    elif query:
        url = f"{endpoint}?{urlencode(query)}"

    status = 0
    response: Any = None
    error = ""
    for attempt in range(retries + 1):
        try:
            request = Request(url, data=body, headers=headers, method=method)
            with urlopen(request, timeout=timeout) as handle:
                status = int(handle.status)
                response = json.loads(handle.read().decode("utf-8"))
            error = ""
            break
        except HTTPError as exc:
            status = int(exc.code)
            raw = exc.read().decode("utf-8", errors="replace")
            error = f"HTTP {status}: {raw[:1000]}"
            if status != 429 and status < 500:
                break
        except (URLError, TimeoutError, json.JSONDecodeError) as exc:
            error = f"{type(exc).__name__}: {exc}"
        if attempt < retries:
            time.sleep(min(2**attempt, 8))

    wrapper = {
        "source": source,
        "source_release": SOURCE_RELEASES.get(source, "not exposed by endpoint"),
        "retrieved_at_utc": utc_now(),
        "endpoint": endpoint,
        "request": request_spec,
        "http_status": status,
        "error": error,
        "response": response,
    }
    cache_path.write_text(json.dumps(wrapper, indent=2, sort_keys=True) + "\n")
    return wrapper


def open_targets_target(
    ensembl_id: str, cache_dir: Path, refresh: bool = False
) -> dict[str, Any]:
    return cached_json_request(
        source="Open Targets",
        endpoint=OPEN_TARGETS_ENDPOINT,
        cache_dir=cache_dir / "open_targets",
        payload={"query": OPEN_TARGETS_QUERY, "variables": {"ensemblId": ensembl_id}},
        refresh=refresh,
    )


def chembl_target_search(symbol: str, cache_dir: Path, refresh: bool = False) -> dict[str, Any]:
    return cached_json_request(
        source="ChEMBL target search",
        endpoint=f"{CHEMBL_ENDPOINT}/target/search.json",
        cache_dir=cache_dir / "chembl",
        query={"q": symbol, "limit": 20},
        refresh=refresh,
    )


def select_human_chembl_target(payload: Mapping[str, Any] | None, symbol: str) -> dict[str, Any] | None:
    """Select an exact human gene-symbol target from a ChEMBL search response."""

    targets = list((payload or {}).get("targets") or [])
    exact: list[dict[str, Any]] = []
    for target in targets:
        if str(target.get("organism", "")).lower() != "homo sapiens":
            continue
        synonyms: set[str] = set()
        for component in target.get("target_components") or []:
            synonym_rows = component.get("target_component_synonyms") or component.get(
                "component_synonyms"
            ) or []
            for item in synonym_rows:
                if str(item.get("syn_type", "")).upper() == "GENE_SYMBOL":
                    synonyms.add(str(item.get("component_synonym", "")).upper())
        if symbol.upper() in synonyms:
            exact.append(dict(target))
    if not exact:
        return None
    exact.sort(
        key=lambda row: (
            str(row.get("target_type", "")).upper() != "SINGLE PROTEIN",
            str(row.get("target_chembl_id", "")),
        )
    )
    return exact[0]


def chembl_mechanisms(
    target_chembl_id: str, cache_dir: Path, refresh: bool = False
) -> dict[str, Any]:
    return cached_json_request(
        source="ChEMBL mechanism",
        endpoint=f"{CHEMBL_ENDPOINT}/mechanism.json",
        cache_dir=cache_dir / "chembl",
        query={"target_chembl_id": target_chembl_id, "limit": 1000},
        refresh=refresh,
    )


def dgidb_genes(
    symbols: Sequence[str], cache_dir: Path, refresh: bool = False
) -> dict[str, Any]:
    query = DGIDB_QUERY_TEMPLATE % json.dumps(list(symbols))
    return cached_json_request(
        source="DGIdb",
        endpoint=DGIDB_ENDPOINT,
        cache_dir=cache_dir / "dgidb",
        payload={"query": query},
        refresh=refresh,
    )


def unwrap_response(wrapper: Mapping[str, Any]) -> Mapping[str, Any]:
    response = wrapper.get("response")
    return response if isinstance(response, Mapping) else {}


def open_targets_target_record(wrapper: Mapping[str, Any]) -> Mapping[str, Any]:
    response = unwrap_response(wrapper)
    data = response.get("data") if isinstance(response.get("data"), Mapping) else {}
    target = data.get("target") if isinstance(data, Mapping) else None
    return target if isinstance(target, Mapping) else {}


def dgidb_gene_records(wrapper: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    response = unwrap_response(wrapper)
    data = response.get("data") if isinstance(response.get("data"), Mapping) else {}
    genes = data.get("genes") if isinstance(data, Mapping) else {}
    nodes = genes.get("nodes") if isinstance(genes, Mapping) else []
    return [node for node in (nodes or []) if isinstance(node, Mapping)]


def cardiovascular_disease_summary(target: Mapping[str, Any]) -> tuple[float, str]:
    rows = ((target.get("associatedDiseases") or {}).get("rows") or [])
    matches: list[tuple[float, str]] = []
    for row in rows:
        disease = row.get("disease") or {}
        name = str(disease.get("name", ""))
        if any(term in name.lower() for term in CARDIOVASCULAR_TERMS):
            matches.append((float(row.get("score") or 0), name))
    matches.sort(reverse=True)
    return (
        max((score for score, _ in matches), default=0.0),
        "; ".join(name for _, name in matches[:8]),
    )


def action_from_dgidb_types(types: Iterable[Mapping[str, Any]]) -> tuple[str, str]:
    raw_parts: list[str] = []
    normalized: set[str] = set()
    for item in types or []:
        action_type = str(item.get("type", "")).strip()
        directionality = str(item.get("directionality", "")).strip()
        raw_parts.append(": ".join(part for part in (action_type, directionality) if part))
        direction_text = directionality.lower()
        if "inhib" in direction_text or "negative" in direction_text:
            normalized.add("inhibit")
        elif "activ" in direction_text or "positive" in direction_text:
            normalized.add("activate")
        else:
            normalized.add(normalize_action(action_type))
    normalized.discard("unknown")
    call = next(iter(normalized)) if len(normalized) == 1 else "unknown"
    return call, "; ".join(part for part in raw_parts if part)


def cache_manifest(cache_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(cache_dir.rglob("*.json")):
        wrapper = json.loads(path.read_text())
        rows.append(
            {
                "source": wrapper.get("source", ""),
                "source_release": wrapper.get(
                    "source_release", SOURCE_RELEASES.get(str(wrapper.get("source", "")), "not exposed")
                ),
                "retrieved_at_utc": wrapper.get("retrieved_at_utc", ""),
                "endpoint": wrapper.get("endpoint", ""),
                "http_status": wrapper.get("http_status", ""),
                "error": wrapper.get("error", ""),
                "cache_file": str(Path("cache/pharmacology") / path.relative_to(cache_dir)),
            }
        )
    return rows
