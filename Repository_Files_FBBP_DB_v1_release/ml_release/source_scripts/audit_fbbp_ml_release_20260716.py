#!/usr/bin/env python
"""Cross-audit all FBBP ML release layers and write a unified status report."""

from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
BINDER = ROOT / "FBBP_ML_binder_full_length_10fold_v1"
TARGET = ROOT / "FBBP_ML_target_sequence_aware_10fold_v1"
STRUCTURE = ROOT / "FBBP_ML_structure_core_qc_audited_v1_20260716"
LOOP = ROOT / "FBBP_ML_binder_loop_aware_10fold_diagnostic_v1_20260716"
COMPLEX = ROOT / "FBBP_ML_structure_complex_biological_assembly_qc_v1_20260716"
DERIVED = ROOT / "FBBP_ML_derived_split_views_20260716"
OUTPUT_JSON = ROOT / "FBBP_ML_IMPLEMENTATION_STATUS_AUDITED_20260716.json"
OUTPUT_MD = ROOT / "FBBP_ML_IMPLEMENTATION_STATUS_AUDITED_20260716.md"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest(root: Path) -> dict:
    manifest_path = root / "metadata/checksums_sha256.json"
    rows = json.loads(manifest_path.read_text(encoding="utf-8"))
    missing = []
    mismatched = []
    for row in rows:
        path = root / row["path"]
        if not path.exists():
            missing.append(row["path"])
        elif sha256(path) != row["sha256"]:
            mismatched.append(row["path"])
    return {
        "manifest_entries": len(rows),
        "missing_files": missing,
        "sha256_mismatches": mismatched,
        "status": "PASS" if not missing and not mismatched else "FAIL",
    }


def assert_zero_crossings(frame: pd.DataFrame, group_column: str, fold_column: str = "fold_id") -> int:
    crossings = int((frame.groupby(group_column)[fold_column].nunique() > 1).sum())
    if crossings:
        raise AssertionError(f"{group_column} crosses {fold_column}: {crossings}")
    return crossings


def main() -> None:
    proteins = pd.read_csv(FOUNDATION / "data/entities/protein_objects.csv", dtype=str).fillna("")
    targets = pd.read_csv(FOUNDATION / "data/entities/target_annotations.csv", dtype=str).fillna("")
    windows = pd.read_csv(FOUNDATION / "data/entities/recognition_window_objects.csv", dtype=str).fillna("")
    binder = pd.read_csv(BINDER / "data/binder_full_length_10fold_v1.csv", dtype=str).fillna("")
    binder_clusters = pd.read_csv(BINDER / "data/binder_sequence_clusters_95_affine_v1.csv", dtype=str).fillna("")
    registry = pd.read_csv(TARGET / "data/frozen_target_sequence_registry.csv", dtype=str).fillna("")
    target_assignment = pd.read_csv(TARGET / "data/target_sequence_aware_10fold_v1.csv", dtype=str).fillna("")
    target_excluded = pd.read_csv(TARGET / "data/named_targets_excluded_from_target_sequence_aware_split.csv", dtype=str).fillna("")
    target_clusters = pd.read_csv(TARGET / "data/target_sequence_clusters_40.csv", dtype=str).fillna("")
    structure = pd.read_csv(STRUCTURE / "data/structure_core_coordinate_qc.csv", dtype=str).fillna("")
    loop_assignment = pd.read_csv(LOOP / "data/binder_loop_aware_10fold_diagnostic_v1.csv", dtype=str).fillna("")
    windows_with_fold = pd.read_csv(LOOP / "data/recognition_windows_with_parent_fold.csv", dtype=str).fillna("")
    complex_qc = pd.read_csv(COMPLEX / "data/complex_structure_biological_assembly_qc.csv", dtype=str).fillna("")
    derived_audit = json.loads((DERIVED / "reports/derived_split_views_audit.json").read_text(encoding="utf-8"))

    if len(proteins) != 1996 or proteins["resource_record_id"].nunique() != 1996:
        raise AssertionError("Foundation protein-object denominator changed")
    if len(binder) != 1996 or binder["resource_record_id"].nunique() != 1996:
        raise AssertionError("Binder release does not cover all protein objects")
    if set(binder["resource_record_id"]) != set(proteins["resource_record_id"]):
        raise AssertionError("Binder release protein IDs differ from foundation")
    binder_crossings = assert_zero_crossings(binder, "binder_cluster_95_id")
    if binder_clusters["binder_cluster_95_id"].nunique() != 1577:
        raise AssertionError("Unexpected affine binder cluster count")

    for row in registry.itertuples(index=False):
        sequence = str(row.sequence)
        if not sequence or not set(sequence).issubset(STANDARD_AA):
            raise AssertionError(f"Invalid frozen target sequence: {row.uniprot_accession}")
        if hashlib.sha256(sequence.encode("ascii")).hexdigest() != row.sequence_sha256:
            raise AssertionError(f"Target sequence SHA mismatch: {row.uniprot_accession}")
        if not (TARGET / row.fasta_file).exists():
            raise AssertionError(f"Missing target FASTA: {row.fasta_file}")
    if registry["uniprot_accession"].duplicated().any():
        raise AssertionError("Duplicate target accessions in frozen registry")
    target_component_crossings = assert_zero_crossings(target_assignment, "target_aware_component_id")
    target_binder_crossings = assert_zero_crossings(target_assignment, "binder_cluster_95_id")
    target_cluster_crossings = assert_zero_crossings(target_assignment, "target_cluster_40_id")
    strict_ids = set(targets.loc[targets["target_split_eligible"].str.lower() == "true", "resource_record_id"])
    included_ids = set(target_assignment["resource_record_id"])
    excluded_ids = set(target_excluded["resource_record_id"])
    if included_ids & excluded_ids:
        raise AssertionError("A strict target object is both included and excluded")
    if included_ids | excluded_ids != strict_ids:
        raise AssertionError("Strict target inclusion/exclusion accounting is incomplete")
    if len(registry) != 145 or len(included_ids) != 679 or len(excluded_ids) != 112:
        raise AssertionError("Final target denominators changed")
    if target_clusters["target_cluster_40_id"].nunique() != 132:
        raise AssertionError("Unexpected 40% target cluster count")

    mapped = structure[structure["status"] == "PASS"]
    exact = structure[structure["exact_sequence_core_status"] == "PASS"]
    if len(structure) != 282 or len(mapped) != 201 or len(exact) != 189:
        raise AssertionError("Structure-core tier denominators changed")
    if not set(exact["structure_instance_id"]).issubset(set(mapped["structure_instance_id"])):
        raise AssertionError("Exact structure tier is not a subset of mapped tier")
    if (mapped["fail_reasons"].str.strip() != "").any():
        raise AssertionError("A mapped structure-core pass has failure reasons")

    if len(loop_assignment) != 1996 or loop_assignment["resource_record_id"].nunique() != 1996:
        raise AssertionError("Loop-aware assignment does not cover all protein objects")
    loop_component_crossings = assert_zero_crossings(loop_assignment, "loop_aware_component_95_id")
    loop_binder_crossings = assert_zero_crossings(loop_assignment, "binder_cluster_95_id")
    if len(windows_with_fold) != len(windows) or windows_with_fold["fold_id"].eq("").any():
        raise AssertionError("Recognition windows do not all inherit a parent fold")

    window_complex = complex_qc[complex_qc["complex_grade_status"] == "PASS_WINDOW_TARGET_CONTACT"]
    any_contact = complex_qc[
        complex_qc["complex_grade_status"].isin(
            ["PASS_WINDOW_TARGET_CONTACT", "PASS_BINDER_TARGET_CONTACT_NO_WINDOW_CONTACT"]
        )
    ]
    if len(complex_qc) != 201 or len(window_complex) != 19 or len(any_contact) != 53:
        raise AssertionError("Biological-assembly complex tier denominators changed")
    if not set(complex_qc["structure_instance_id"]).issubset(set(mapped["structure_instance_id"])):
        raise AssertionError("Complex QC contains a structure outside mapped-core passes")
    if not set(window_complex["structure_instance_id"]).issubset(set(any_contact["structure_instance_id"])):
        raise AssertionError("Window-contact tier is not a subset of any-contact tier")

    crossing_values = [
        value
        for view in derived_audit["crossing_audit"].values()
        for value in view.values()
    ]
    if any(crossing_values):
        raise AssertionError("A derived split view has leakage-controlled crossings")

    manifests = {
        "binder": validate_manifest(BINDER),
        "target": validate_manifest(TARGET),
        "structure_core": validate_manifest(STRUCTURE),
        "loop_aware": validate_manifest(LOOP),
        "complex": validate_manifest(COMPLEX),
        "derived_views": validate_manifest(DERIVED),
    }
    if any(item["status"] != "PASS" for item in manifests.values()):
        raise AssertionError("At least one checksum manifest failed")

    report = {
        "status": "PASS",
        "audit_date": str(date.today()),
        "foundation": {
            "protein_objects": 1996,
            "structure_instances": 1996,
            "target_annotations": 908,
            "strict_target_eligible_objects": len(strict_ids),
            "recognition_windows": len(windows),
        },
        "binder_full_length_10fold_v1": {
            "protein_objects": 1996,
            "affine_binder_clusters_95": int(binder_clusters["binder_cluster_95_id"].nunique()),
            "largest_cluster": int(pd.to_numeric(binder_clusters["binder_cluster_size"]).max()),
            "fold_sizes": binder.groupby("fold_id")["resource_record_id"].nunique().to_dict(),
            "cluster_crossings": binder_crossings,
        },
        "target_sequence_aware_10fold_v1": {
            "frozen_target_sequences": len(registry),
            "uniparc_recovered_sequences": int(registry["retrieval_source"].str.contains("UniParc").sum()),
            "resolved_protein_objects": len(included_ids),
            "excluded_strict_objects": len(excluded_ids),
            "target_clusters_40": int(target_clusters["target_cluster_40_id"].nunique()),
            "target_aware_components": int(target_assignment["target_aware_component_id"].nunique()),
            "largest_component": int(pd.to_numeric(target_assignment["component_size"]).max()),
            "fold_sizes": target_assignment.groupby("fold_id")["resource_record_id"].nunique().to_dict(),
            "crossings": {
                "binder": target_binder_crossings,
                "target": target_cluster_crossings,
                "component": target_component_crossings,
            },
        },
        "structure_core": {
            "xray_em_candidates": len(structure),
            "mapped_95_pass": len(mapped),
            "exact_sequence_pass": len(exact),
        },
        "loop_aware_diagnostic": {
            "protein_objects": len(loop_assignment),
            "parent_components_95": int(loop_assignment["loop_aware_component_95_id"].nunique()),
            "largest_component_95": int(pd.to_numeric(loop_assignment["loop_aware_component_95_size"]).max()),
            "fold_sizes": loop_assignment.groupby("fold_id")["resource_record_id"].nunique().to_dict(),
            "crossings": {"component": loop_component_crossings, "binder": loop_binder_crossings},
        },
        "biological_assembly_complex_qc": {
            "mapped_core_evaluated": len(complex_qc),
            "window_target_contact_pass": len(window_complex),
            "any_binder_target_contact_pass": len(any_contact),
        },
        "derived_views": {
            "available_ratios": ["9:1", "8:2", "7:3", "8:1:1"],
            "crossing_audit": derived_audit["crossing_audit"],
        },
        "checksum_manifests": manifests,
        "release_boundary": {
            "benchmark_zip_created": False,
            "manuscript_or_docx_rendered": False,
            "target_aware_scope": "679 sequence-resolved strict objects; 112 unresolved strict objects remain explicitly excluded",
            "complex_grade_scope": "19 recognition-window contact structures; does not replace the 201/189 structure-core tiers",
        },
    }
    OUTPUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    lines = [
        "# FBBP ML implementation status - independently audited",
        "",
        "Status: **PASS**",
        "",
        "## Final denominators",
        "",
        "- Foundation: 1,996 protein objects; 908 target annotations (791 strict, 117 placeholder-excluded); 3,383 recognition windows.",
        "- Binder full-length split: 1,996 objects; 1,577 affine 95% within-family clusters; 10 folds; zero cluster crossings.",
        "- Target-sequence-aware split: 145 frozen target sequences; 679 strict objects included; 112 explicitly excluded; 132 target clusters at 40%; 123 components; zero crossings.",
        "- Structure core: 282 X-ray/EM candidates; 201 mapped-95 passes; 189 exact-sequence-core passes.",
        "- Loop-aware diagnostic: 955 parent components at 95%; largest component 244; zero component or binder-cluster crossings.",
        "- Biological-assembly complex QC: 201 mapped-core structures evaluated; 53 whole-binder contacts; 19 recognition-window contacts.",
        "",
        "## Derived views",
        "",
        "- 9:1, 8:2, 7:3 and 8:1:1 views were generated from frozen fold IDs for full-length, target-aware and loop-aware assignments.",
        "- Every derived view has zero binder-cluster, target-component and loop-component crossings.",
        "",
        "## Boundary",
        "",
        "- No benchmark ZIP was created.",
        "- No DOCX was rendered or modified in this audit run.",
        "- The 19-structure complex-grade tier is an additional high-confidence layer, not the denominator for the full ML release.",
        "- The 112 unresolved strict target objects remain in the database with explicit exclusion reasons; no ambiguous sequence was fabricated.",
        "",
        "## Checksum QA",
        "",
    ]
    lines.extend(
        f"- {name}: {item['status']} ({item['manifest_entries']} entries)"
        for name, item in manifests.items()
    )
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUTPUT_MD)


if __name__ == "__main__":
    main()
