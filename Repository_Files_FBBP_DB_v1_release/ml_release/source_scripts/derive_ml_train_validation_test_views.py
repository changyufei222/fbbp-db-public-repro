#!/usr/bin/env python
"""Derive fixed ratio views from frozen FBBP fold IDs without moving any cluster."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BINDER = ROOT / "FBBP_ML_binder_full_length_10fold_v1/data/binder_full_length_10fold_v1.csv"
TARGET = ROOT / "FBBP_ML_target_sequence_aware_10fold_v1/data/target_sequence_aware_10fold_v1.csv"
STRUCTURE_QC = ROOT / "FBBP_ML_structure_core_qc_audited_v1_20260716/data/structure_core_coordinate_qc.csv"
LOOP_AWARE = ROOT / "FBBP_ML_binder_loop_aware_10fold_diagnostic_v1_20260716/data/binder_loop_aware_10fold_diagnostic_v1.csv"
COMPLEX_QC = ROOT / "FBBP_ML_structure_complex_biological_assembly_qc_v1_20260716/data/complex_structure_biological_assembly_qc.csv"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_derived_split_views_20260716"
VIEW_MAPPINGS = {
    "9_1": {fold_id: ("train" if fold_id <= 8 else "test") for fold_id in range(10)},
    "8_2": {fold_id: ("train" if fold_id <= 7 else "test") for fold_id in range(10)},
    "7_3": {fold_id: ("train" if fold_id <= 6 else "test") for fold_id in range(10)},
    "8_1_1": {
        fold_id: ("train" if fold_id <= 7 else "validation" if fold_id == 8 else "test")
        for fold_id in range(10)
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_views(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for label, mapping in VIEW_MAPPINGS.items():
        result[f"split_role_{label}"] = result["fold_id"].map(mapping)
    return result


def split_counts(frame: pd.DataFrame, id_column: str) -> dict:
    return {
        label: frame.groupby(f"split_role_{label}")[id_column].nunique().to_dict()
        for label in VIEW_MAPPINGS
    }


def main() -> None:
    output = DEFAULT_OUTPUT
    if output.exists():
        shutil.rmtree(output)
    (output / "data").mkdir(parents=True)
    (output / "reports").mkdir(parents=True)
    (output / "metadata").mkdir(parents=True)
    binder = pd.read_csv(BINDER)
    target = pd.read_csv(TARGET)
    qc = pd.read_csv(STRUCTURE_QC)
    complex_qc = pd.read_csv(COMPLEX_QC)
    loop_aware = pd.read_csv(LOOP_AWARE)
    binder = add_views(binder)
    target = add_views(target)
    loop_aware = add_views(loop_aware)
    for label in VIEW_MAPPINGS:
        binder.to_csv(output / f"data/binder_full_length_{label}.csv", index=False)
        target.to_csv(output / f"data/target_sequence_aware_{label}.csv", index=False)
        loop_aware.to_csv(output / f"data/binder_loop_aware_{label}.csv", index=False)
    structure_core = qc[qc["status"] == "PASS"].merge(
        binder[["resource_record_id", "fold_id"] + [f"split_role_{label}" for label in VIEW_MAPPINGS]],
        on="resource_record_id", how="left", validate="one_to_one"
    )
    if structure_core["fold_id"].isna().any():
        raise RuntimeError("A structure-core record lacks a parent binder fold")
    structure_core.to_csv(
        output / "data/experimental_structure_core_mapped_95_parent_binder_fold.csv",
        index=False,
    )
    exact_structure_core = structure_core[
        structure_core["exact_sequence_core_status"] == "PASS"
    ].copy()
    exact_structure_core.to_csv(
        output / "data/experimental_structure_core_exact_sequence_parent_binder_fold.csv",
        index=False,
    )
    complex_with_fold = complex_qc.merge(
        binder[["resource_record_id", "fold_id"] + [f"split_role_{label}" for label in VIEW_MAPPINGS]],
        on="resource_record_id",
        how="left",
        validate="one_to_one",
    )
    if complex_with_fold["fold_id"].isna().any():
        raise RuntimeError("A biological-assembly complex record lacks a parent binder fold")
    complex_window = complex_with_fold[
        complex_with_fold["complex_grade_status"] == "PASS_WINDOW_TARGET_CONTACT"
    ].copy()
    complex_contact = complex_with_fold[
        complex_with_fold["complex_grade_status"].isin(
            ["PASS_WINDOW_TARGET_CONTACT", "PASS_BINDER_TARGET_CONTACT_NO_WINDOW_CONTACT"]
        )
    ].copy()
    complex_window.to_csv(
        output / "data/complex_grade_window_contact_parent_binder_fold.csv", index=False
    )
    complex_contact.to_csv(
        output / "data/complex_grade_any_binder_contact_parent_binder_fold.csv", index=False
    )
    crossing_audit = {}
    for label in VIEW_MAPPINGS:
        role_column = f"split_role_{label}"
        binder_crossings = binder.groupby("binder_cluster_95_id")[role_column].nunique()
        target_crossings = target.groupby("target_aware_component_id")[role_column].nunique()
        loop_crossings = loop_aware.groupby("loop_aware_component_95_id")[role_column].nunique()
        crossing_audit[label] = {
            "binder_cluster_crossings": int((binder_crossings > 1).sum()),
            "target_component_crossings": int((target_crossings > 1).sum()),
            "loop_aware_component_crossings": int((loop_crossings > 1).sum()),
        }
        if any(crossing_audit[label].values()):
            raise RuntimeError(f"A leakage-controlled group crosses derived split roles in {label}")
    audit = {
        "release_date": str(date.today()),
        "ratio_mappings": VIEW_MAPPINGS,
        "binder_protein_objects": int(binder["resource_record_id"].nunique()),
        "binder_split_counts": split_counts(binder, "resource_record_id"),
        "target_sequence_resolved_protein_objects": int(target["resource_record_id"].nunique()),
        "target_pair_rows": int(len(target)),
        "target_split_counts": split_counts(target, "resource_record_id"),
        "loop_aware_protein_objects": int(loop_aware["resource_record_id"].nunique()),
        "loop_aware_split_counts": split_counts(loop_aware, "resource_record_id"),
        "mapped_95_structure_core_records": int(len(structure_core)),
        "mapped_95_structure_core_split_counts": split_counts(structure_core, "resource_record_id"),
        "exact_sequence_structure_core_records": int(len(exact_structure_core)),
        "exact_sequence_structure_core_split_counts": split_counts(exact_structure_core, "resource_record_id"),
        "complex_grade_window_contact_records": int(len(complex_window)),
        "complex_grade_window_contact_split_counts": split_counts(complex_window, "resource_record_id"),
        "complex_grade_any_binder_contact_records": int(len(complex_contact)),
        "complex_grade_any_binder_contact_split_counts": split_counts(complex_contact, "resource_record_id"),
        "crossing_audit": crossing_audit,
        "not_created": ["benchmark_zip", "complex_grade_structure_benchmark"],
    }
    (output / "reports/derived_split_views_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (output / "README.md").write_text(
        "# FBBP ML derived split views\n\n"
        "The 9:1, 8:2, 7:3 and 8:1:1 roles in this directory are derived only from frozen 10-fold IDs. "
        "No binder cluster, target cluster, or target-aware component was moved. "
        "Structure-core mapped-95 and exact-sequence tiers inherit their parent binder fold.\n",
        encoding="utf-8",
    )
    files = sorted(path for path in output.rglob("*") if path.is_file())
    (output / "metadata/checksums_sha256.json").write_text(
        json.dumps([{"path": path.relative_to(output).as_posix(), "sha256": sha256(path)} for path in files], indent=2),
        encoding="utf-8",
    )
    print(output)


if __name__ == "__main__":
    main()
