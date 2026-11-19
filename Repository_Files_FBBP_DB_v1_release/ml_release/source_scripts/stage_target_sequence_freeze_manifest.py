#!/usr/bin/env python
"""Stage auditable target-sequence retrieval inputs without creating a target-aware split."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_target_sequence_freeze_staging_20260714"
ACCESSION = re.compile(r"^[A-Z0-9]{6,10}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_accessions(value: str) -> list[str]:
    return [part.strip().upper() for part in str(value).split(";") if part.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        shutil.rmtree(args.output)
    (args.output / "data").mkdir(parents=True)
    (args.output / "reports").mkdir(parents=True)
    (args.output / "metadata").mkdir(parents=True)

    targets = pd.read_csv(
        FOUNDATION / "data/entities/target_annotations.csv", dtype=str
    ).fillna("")
    strict = targets[targets["target_split_eligible"].str.lower() == "true"].copy()
    accession_rows = strict[strict["uniprot_id"].str.strip() != ""].copy()
    accession_rows["uniprot_accession"] = accession_rows["uniprot_id"].map(split_accessions)
    accession_rows = accession_rows.explode("uniprot_accession")
    accession_rows["accession_format_valid"] = accession_rows["uniprot_accession"].map(
        lambda value: bool(ACCESSION.fullmatch(value))
    )
    manifest = (
        accession_rows.groupby("uniprot_accession", as_index=False)
        .agg(
            target_variant_count=("target_variant_id", "nunique"),
            strict_target_record_count=("resource_record_id", "nunique"),
            target_variant_ids=("target_variant_id", lambda values: ";".join(sorted(set(values)))),
            resource_record_ids=("resource_record_id", lambda values: ";".join(sorted(set(values)))),
            accession_format_valid=("accession_format_valid", "all"),
        )
        .sort_values("uniprot_accession")
    )
    manifest["retrieval_source"] = "UniProtKB REST /uniprotkb/{accession}.fasta"
    manifest["retrieval_status"] = manifest["accession_format_valid"].map(
        lambda valid: "pending_remote_retrieval" if valid else "manual_accession_normalization_required"
    )
    manifest["frozen_sequence_sha256"] = ""
    manifest["source_record_version"] = ""
    manifest["retrieval_date"] = ""
    manifest.to_csv(args.output / "data/target_accession_freeze_manifest.csv", index=False)

    named = strict[(strict["uniprot_id"].str.strip() == "")].copy()
    named = named[[
        "resource_record_id", "target_variant_id", "target_concept_id", "gene_name_species",
        "gene_name", "protein_name", "species_name", "strict_target_key",
    ]]
    named["status"] = "named_target_without_frozen_sequence"
    named.to_csv(args.output / "data/named_targets_pending_sequence_resolution.csv", index=False)

    audit = {
        "status": "STAGED_NOT_FROZEN",
        "release_date": str(date.today()),
        "strict_target_records": int(len(strict)),
        "records_with_uniprot_accession": int(len(accession_rows["resource_record_id"].unique())),
        "unique_accessions_after_semicolon_split": int(manifest["uniprot_accession"].nunique()),
        "valid_accession_format": int(manifest["accession_format_valid"].sum()),
        "manual_accession_normalization_required": int((~manifest["accession_format_valid"]).sum()),
        "named_target_records_pending_sequence_resolution": int(len(named)),
        "not_created": [
            "target_sequence_fasta",
            "frozen_target_sequence_registry",
            "target_sequence_clusters_40",
            "target_sequence_aware_10fold_v1",
        ],
    }
    (args.output / "reports/target_sequence_freeze_staging_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    report = [
        "# FBBP target-sequence freeze staging audit",
        "",
        "Status: **STAGED - target sequences are not frozen in this directory**",
        "",
        f"- Strict target records: {audit['strict_target_records']}",
        f"- Records with UniProt accessions: {audit['records_with_uniprot_accession']}",
        f"- Unique accessions after splitting multi-accession fields: {audit['unique_accessions_after_semicolon_split']}",
        f"- Accessions requiring manual normalization: {audit['manual_accession_normalization_required']}",
        f"- Name-only strict target records: {audit['named_target_records_pending_sequence_resolution']}",
        "- Fetch each valid accession from the declared source, retain the raw FASTA header, sequence hash, retrieval date and source record version, then rerun target clustering.",
        "- Name-only records cannot enter the sequence-controlled target-aware graph until sequence identity is resolved.",
    ]
    (args.output / "reports/target_sequence_freeze_staging_audit.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    readme = [
        "# FBBP target-sequence freeze staging",
        "",
        "This directory stages accession and name-resolution inputs for target sequence freezing.",
        "It is not a frozen sequence registry and does not contain target clusters or target-aware folds.",
    ]
    (args.output / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    manifest_files = sorted(path for path in args.output.rglob("*") if path.is_file())
    checksums = [
        {"path": path.relative_to(args.output).as_posix(), "sha256": sha256(path)}
        for path in manifest_files
    ]
    (args.output / "metadata/checksums_sha256.json").write_text(
        json.dumps(checksums, indent=2), encoding="utf-8"
    )
    print(args.output)


if __name__ == "__main__":
    main()
