#!/usr/bin/env python3
"""Download and cache Ensembl mouse-human one-to-one orthologues."""

from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
CACHE_DIR = Path(__file__).resolve().parents[1] / "cache"
CACHE_JSON = CACHE_DIR / "ensembl_mouse_human_orthologues.json"
OUT_TSV = ROOT / "results" / "tables" / "gse214611_mouse_human_one_to_one_orthologues.tsv"
ENDPOINT = "https://www.ensembl.org/biomart/martservice"

QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="mmusculus_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="hsapiens_homolog_ensembl_gene"/>
    <Attribute name="hsapiens_homolog_associated_gene_name"/>
    <Attribute name="hsapiens_homolog_orthology_type"/>
    <Attribute name="hsapiens_homolog_orthology_confidence"/>
  </Dataset>
</Query>"""


def fetch_biomart() -> dict[str, object]:
    payload = urllib.parse.urlencode({"query": QUERY}).encode()
    request = urllib.request.Request(
        ENDPOINT,
        data=payload,
        headers={"User-Agent": "spatial-mi-border-zone-model/0.2"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        text = response.read().decode("utf-8")
        headers = dict(response.headers.items())
    lines = text.splitlines()
    reader = csv.DictReader(lines, delimiter="\t")
    records = []
    for row in reader:
        if row.get("Human homology type") != "ortholog_one2one":
            continue
        if row.get("Human orthology confidence [0 low, 1 high]") != "1":
            continue
        if not row.get("Gene stable ID") or not row.get("Human gene stable ID"):
            continue
        records.append(
            {
                "mouse_gene_id": row["Gene stable ID"],
                "mouse_gene_name": row.get("Gene name", ""),
                "human_gene_id": row["Human gene stable ID"],
                "human_gene_name": row.get("Human gene name", ""),
                "orthology_type": row["Human homology type"],
                "orthology_confidence": 1,
            }
        )
    if len(records) < 10000:
        raise RuntimeError(f"Unexpectedly small Ensembl one-to-one mapping: {len(records)} records")
    return {
        "source": "Ensembl BioMart",
        "endpoint": ENDPOINT,
        "retrieval_date": date.today().isoformat(),
        "response_headers": headers,
        "query": QUERY,
        "records": records,
    }


def load_or_fetch(refresh: bool = False) -> dict[str, object]:
    if CACHE_JSON.exists() and not refresh:
        return json.loads(CACHE_JSON.read_text())
    payload = fetch_biomart()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_JSON.write_text(json.dumps(payload, ensure_ascii=True, indent=2) + "\n")
    return payload


def main() -> None:
    payload = load_or_fetch()
    records = list(payload["records"])
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", newline="") as handle:
        fields = [
            "mouse_gene_id",
            "mouse_gene_name",
            "human_gene_id",
            "human_gene_name",
            "orthology_type",
            "orthology_confidence",
            "mapping_source",
            "retrieval_date",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in records:
            writer.writerow(
                {
                    **row,
                    "mapping_source": payload["source"],
                    "retrieval_date": payload["retrieval_date"],
                }
            )
    print(f"Wrote {OUT_TSV} ({len(records)} one-to-one orthologues)")


if __name__ == "__main__":
    main()
