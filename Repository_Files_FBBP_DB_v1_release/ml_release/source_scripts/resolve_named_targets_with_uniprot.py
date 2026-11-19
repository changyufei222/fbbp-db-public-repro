#!/usr/bin/env python
"""Resolve name-only FBBP targets through conservative UniProtKB search queries."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_named_target_uniprot_resolution_v2_20260716"
SPECIES_TAXON = {
    "Homo sapiens": "9606",
    "Mus musculus": "10090",
    "Rattus norvegicus": "10116",
    "Oryctolagus cuniculus": "9986",
    "Plasmodium falciparum": "5833",
    "Macaca fascicularis": "9541",
    "Clostridium thermocellum": "1515",
    "Mycobacterium tuberculosis": "1773",
    "HIV-1": "11676",
}

# Only aliases with an explicit organism and an unambiguous canonical gene (or
# explicitly named multi-protein complex) are eligible for deterministic rescue.
CURATED_GENE_TARGETS = {
    "TVAR_0026": [("ELANE", "9606")],
    "TVAR_0022": [("KLKB1", "9606")],
    "TVAR_0024": [("MSTN", "9606")],
    "TVAR_0034": [("ITGAV", "9606"), ("ITGB3", "9606")],
    "TVAR_0037": [("ITGAV", "9606"), ("ITGB3", "9606"), ("ITGB5", "9606"), ("ITGA5", "9606"), ("ITGB1", "9606")],
    "TVAR_0038": [("ITGAV", "9606"), ("ITGB3", "9606"), ("ITGB5", "9606"), ("ITGA5", "9606"), ("ITGB1", "9606")],
    "TVAR_0047": [("IL2RA", "9606")],
    "TVAR_0049": [("GARA", "1773")],
    "TVAR_0052": [("ALB", "9606")],
    "TVAR_0056": [("IGHE", "9606")],
    "TVAR_0062": [("FOLH1", "9541")],
    "TVAR_0067": [("FCGR2A", "9606")],
    "TVAR_0070": [("MET", "9606")],
    "TVAR_0071": [("TFRC", "9606")],
    "TVAR_0072": [("AMA1", "5833")],
    "TVAR_0076": [("CALR", "9986")],
    "TVAR_0082": [("PLG", "9606")],
    "TVAR_0086": [("NFKBIA", "9606")],
    "TVAR_0088": [("TFRC", "9606"), ("EPCAM", "9606"), ("ALB", "9606")],
    "TVAR_0089": [("MET", "9606")],
    "TVAR_0090": [("MET", "9606")],
    "TVAR_0099": [("KCNA1", "10116")],
    "TVAR_0104": [("F7", "9606")],
    "TVAR_0105": [("F7", "9606")],
    "TVAR_0107": [("F3", "9606"), ("F7", "9606")],
}

# Gene-field searches include synonyms and can return multiple reviewed entries
# for otherwise unambiguous labels. These accessions are pinned and then
# verified by exact accession query before use.
CURATED_ACCESSION_TARGETS = {
    "TVAR_0052": ["P02768"],
    "TVAR_0070": ["P08581"],
    "TVAR_0088": ["P02786", "P16422", "P02768"],
    "TVAR_0089": ["P08581"],
    "TVAR_0090": ["P08581"],
}


def clean_gene(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(value)).upper()


def search(query: str) -> tuple[int, list[dict]]:
    url = "https://rest.uniprot.org/uniprotkb/search?" + urlencode({
        "query": query,
        "format": "tsv",
        "fields": "accession,id,protein_name,organism_name,gene_names,reviewed",
        "size": "5",
    })
    request = Request(url, headers={"User-Agent": "FBBP-DB-ML/20260714"})
    with urlopen(request, timeout=60) as response:
        total = int(response.headers.get("X-Total-Results", "0"))
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    return total, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--refresh-curated", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.resume:
        shutil.rmtree(args.output)
    args.output.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(FOUNDATION / "data/entities/target_annotations.csv", dtype=str).fillna("")
    named = targets[(targets["target_split_eligible"].str.lower() == "true") & (targets["uniprot_id"].str.strip() == "")].copy()
    variants = named.sort_values("resource_record_id").drop_duplicates("target_variant_id")
    progress_path = args.output / "named_target_uniprot_resolution.csv"
    if args.resume and progress_path.exists():
        results = pd.read_csv(progress_path, dtype=str).fillna("").to_dict("records")
        if args.refresh_curated:
            results = [
                item
                for item in results
                if not (
                    item["target_variant_id"] in (CURATED_GENE_TARGETS.keys() | CURATED_ACCESSION_TARGETS.keys())
                    and item["resolution_status"] == "manual_review_required"
                )
            ]
    else:
        results = []
    completed = {item["target_variant_id"] for item in results}
    for row in variants.itertuples(index=False):
        if row.target_variant_id in completed:
            continue
        species_taxon = SPECIES_TAXON.get(row.species_name, "")
        gene = clean_gene(row.gene_name or row.gene_name_species)
        protein = str(row.protein_name).strip()
        attempts = []
        curated_targets = CURATED_GENE_TARGETS.get(row.target_variant_id, [])
        curated_accessions = CURATED_ACCESSION_TARGETS.get(row.target_variant_id, [])
        if not curated_targets and not curated_accessions:
            if gene:
                attempts.append(("gene", f"gene:{gene}"))
            if protein and len(protein) >= 4:
                attempts.append(("protein_name", f'protein_name:"{protein}"'))
        selected = None
        query_history = []
        for method, base_query in attempts:
            species_clause = f" AND organism_id:{species_taxon}" if species_taxon else ""
            # Prefer one reviewed, species-matched record before considering all records.
            for label, query in [
                ("reviewed", base_query + species_clause + " AND reviewed:true"),
                ("all", base_query + species_clause),
            ]:
                try:
                    total, rows = search(query)
                    query_history.append(f"{method}_{label}:{query}:{total}")
                except Exception as error:
                    query_history.append(f"{method}_{label}:{query}:ERROR:{error}")
                    continue
                if total == 1 and len(rows) == 1:
                    selected = rows[0]
                    break
            if selected:
                break
        selected_rows = [selected] if selected else []
        if not selected_rows and curated_accessions:
            accession_rows = []
            for accession in curated_accessions:
                query = f"accession:{accession}"
                try:
                    total, rows = search(query)
                    query_history.append(f"curated_accession:{query}:{total}")
                except Exception as error:
                    query_history.append(f"curated_accession:{query}:ERROR:{error}")
                    accession_rows = []
                    break
                if total != 1 or len(rows) != 1:
                    accession_rows = []
                    break
                accession_rows.append(rows[0])
            selected_rows = accession_rows
        if not selected_rows and curated_targets:
            curated_rows = []
            for curated_gene, curated_taxon in curated_targets:
                component = None
                for label, query in [
                    ("reviewed", f"gene_exact:{curated_gene} AND organism_id:{curated_taxon} AND reviewed:true"),
                    ("all", f"gene_exact:{curated_gene} AND organism_id:{curated_taxon}"),
                ]:
                    try:
                        total, rows = search(query)
                        query_history.append(f"curated_gene_{label}:{query}:{total}")
                    except Exception as error:
                        query_history.append(f"curated_gene_{label}:{query}:ERROR:{error}")
                        continue
                    if total == 1 and len(rows) == 1:
                        component = rows[0]
                        break
                if component is None:
                    curated_rows = []
                    break
                curated_rows.append(component)
            selected_rows = curated_rows
        if selected_rows:
            status = "auto_mapped_unique" if len(selected_rows) == 1 and selected is not None else (
                "curated_mapped_unique" if len(selected_rows) == 1 else "curated_mapped_multi_component"
            )
        else:
            status = "manual_review_required"
        results.append({
            "target_variant_id": row.target_variant_id,
            "target_concept_id": row.target_concept_id,
            "gene_name_species": row.gene_name_species,
            "gene_name": row.gene_name,
            "protein_name": row.protein_name,
            "species_name": row.species_name,
            "species_taxon": species_taxon,
            "resolution_status": status,
            "selected_uniprot_accession": ";".join(item.get("Entry", "") for item in selected_rows),
            "selected_entry_name": ";".join(item.get("Entry Name", "") for item in selected_rows),
            "selected_protein_name": ";".join(item.get("Protein names", "") for item in selected_rows),
            "selected_organism": ";".join(item.get("Organism", "") for item in selected_rows),
            "query_history": " | ".join(query_history),
            "resolution_date": str(date.today()),
        })
        pd.DataFrame(results).sort_values("target_variant_id").to_csv(progress_path, index=False)
    output = pd.DataFrame(results)
    output = output.sort_values("target_variant_id")
    output.to_csv(progress_path, index=False)
    summary = output["resolution_status"].value_counts().rename_axis("resolution_status").reset_index(name="target_variants")
    summary.to_csv(args.output / "named_target_uniprot_resolution_summary.csv", index=False)
    print(args.output)


if __name__ == "__main__":
    main()
