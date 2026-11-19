"""Build self-contained ML-ready sequence and structure layers for FBBP ML-grade v1."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser


PACKAGE = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(r"<local_path_removed>")
SOURCE_COORDINATES = SOURCE_ROOT / "结构序列整理/final_sequence_structure_web_package"
ASSEMBLY_CACHE = SOURCE_ROOT / "FBBP_RCSB_BIOASSEMBLY_CACHE_20260716"

FOUNDATION = PACKAGE / "01_foundation/FBBP_ML_foundation_20260714"
BINDER = PACKAGE / "02_binder_full_length_10fold_v1/FBBP_ML_binder_full_length_10fold_v1"
TARGET = PACKAGE / "03_target_sequence_aware_10fold_v1/FBBP_ML_target_sequence_aware_10fold_v1"
STRUCTURE = PACKAGE / "04_structure_core_qc_audited_v1/FBBP_ML_structure_core_qc_audited_v1_20260716"
COMPLEX = PACKAGE / "05_biological_assembly_complex_qc_v1/FBBP_ML_structure_complex_biological_assembly_qc_v1_20260716"
DERIVED = PACKAGE / "07_derived_split_views/FBBP_ML_derived_split_views_20260716"
OUTPUT = PACKAGE / "12_ml_ready_data"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def assert_no_missing(frame: pd.DataFrame, columns: list[str], context: str) -> None:
    missing = {
        column: int((frame[column].astype(str).str.strip() == "").sum())
        for column in columns
        if int((frame[column].astype(str).str.strip() == "").sum())
    }
    if missing:
        raise ValueError(f"{context}: required values missing: {missing}")


def write_table(frame: pd.DataFrame, stem: str) -> dict[str, str | int | bool]:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT / f"{stem}.csv"
    frame.to_csv(csv_path, index=False)
    result: dict[str, str | int | bool] = {
        "rows": int(len(frame)),
        "csv": csv_path.relative_to(PACKAGE).as_posix(),
        "parquet_written": False,
    }
    try:
        parquet_path = OUTPUT / f"{stem}.parquet"
        frame.to_parquet(parquet_path, index=False)
        result["parquet"] = parquet_path.relative_to(PACKAGE).as_posix()
        result["parquet_written"] = True
    except (ImportError, ModuleNotFoundError, ValueError):
        pass
    return result


def parse_pdb(path: Path, compressed: bool = False) -> tuple[int, str]:
    parser = PDBParser(QUIET=True)
    if compressed:
        with gzip.open(path, "rt") as handle:
            structure = parser.get_structure(path.stem, handle)
    else:
        structure = parser.get_structure(path.stem, path)
    model = next(iter(structure))
    atoms = sum(1 for _ in structure.get_atoms())
    chains = ";".join(sorted(chain.id for chain in model))
    if atoms == 0 or not chains:
        raise ValueError(f"Unparseable or empty coordinate file: {path}")
    return atoms, chains


def copy_coordinate(source: Path, destination: Path) -> tuple[str, int, int, str]:
    if not source.exists():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists() or sha256(source) != sha256(destination):
        shutil.copy2(source, destination)
    compressed = destination.suffix.lower() == ".gz"
    atoms, chains = parse_pdb(destination, compressed=compressed)
    return sha256(destination), destination.stat().st_size, atoms, chains


def main() -> None:
    proteins = read_csv(FOUNDATION / "data/entities/protein_objects.csv")
    assert_no_missing(proteins, ["resource_record_id", "sequence_object_id", "final_sequence"], "protein objects")
    if not proteins["final_sequence"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", na=False).all():
        raise ValueError("Protein sequences contain missing or noncanonical characters")

    binder_folds = read_csv(BINDER / "data/binder_full_length_10fold_v1.csv")
    binder_view = read_csv(DERIVED / "data/binder_full_length_8_1_1.csv")
    binder_ready = proteins.merge(
        binder_folds[["resource_record_id", "binder_cluster_95_id", "binder_cluster_size", "fold_id"]],
        on="resource_record_id",
        how="inner",
        validate="one_to_one",
    ).merge(
        binder_view[["resource_record_id", "split_role_8_1_1"]],
        on="resource_record_id",
        how="inner",
        validate="one_to_one",
    )
    if len(binder_ready) != 1996:
        raise ValueError(f"Expected 1,996 binder rows, found {len(binder_ready)}")
    binder_ready = binder_ready.rename(
        columns={
            "final_sequence": "binder_sequence",
            "sequence_length": "binder_sequence_length",
            "split_role_8_1_1": "split_role",
        }
    )
    binder_ready["task_id"] = "fbbp_scaffold_family_classification_v1"
    binder_ready["supervised_label"] = binder_ready["scaffold_category"]
    binder_ready["label_definition"] = "normalized_scaffold_category"
    binder_ready["source_release"] = "binder_full_length_10fold_v1"
    assert_no_missing(binder_ready, ["binder_sequence", "fold_id", "split_role", "supervised_label"], "binder ML-ready")
    binder_summary = write_table(binder_ready, "binder_full_length_ml_ready_8_1_1_v1")

    target_folds = read_csv(TARGET / "data/target_sequence_aware_10fold_v1.csv")
    target_view = read_csv(DERIVED / "data/target_sequence_aware_8_1_1.csv")
    registry = read_csv(TARGET / "data/frozen_target_sequence_registry.csv")
    target_ready = target_folds.merge(
        proteins[["resource_record_id", "sequence_object_id", "canonical_name", "scaffold_category", "organism", "sequence_length", "final_sequence"]],
        on="resource_record_id",
        how="inner",
        validate="many_to_one",
    ).merge(
        registry[["target_sequence_id", "uniprot_accession", "sequence_length", "sequence", "sequence_sha256", "fasta_file"]],
        on=["target_sequence_id", "uniprot_accession"],
        how="inner",
        suffixes=("_binder", "_target"),
        validate="many_to_one",
    ).merge(
        target_view[["resource_record_id", "target_sequence_id", "split_role_8_1_1"]],
        on=["resource_record_id", "target_sequence_id"],
        how="inner",
        validate="one_to_one",
    )
    if len(target_ready) != 710 or target_ready["resource_record_id"].nunique() != 679:
        raise ValueError(
            "Expected 710 target-aware protein-target pairs from 679 unique protein objects, "
            f"found {len(target_ready)} pairs from {target_ready['resource_record_id'].nunique()} objects"
        )
    target_ready = target_ready.rename(
        columns={
            "final_sequence": "binder_sequence",
            "sequence_length_binder": "binder_sequence_length",
            "sequence": "target_sequence",
            "sequence_length_target": "target_sequence_length",
            "split_role_8_1_1": "split_role",
        }
    )
    target_ready["task_id"] = "fbbp_target_association_retrieval_v1"
    target_ready["curated_target_association_label"] = 1
    target_ready["label_caveat"] = "positive curated association only; no experimentally verified negative pairs are supplied"
    target_ready["source_release"] = "target_sequence_aware_10fold_v1"
    assert_no_missing(target_ready, ["binder_sequence", "target_sequence", "fold_id", "split_role"], "target ML-ready")
    target_summary = write_table(target_ready, "target_sequence_aware_pair_ml_ready_8_1_1_v1")

    qc = read_csv(STRUCTURE / "data/structure_core_coordinate_qc.csv")
    mapped = qc[qc["status"] == "PASS"].copy()
    exact = mapped[mapped["exact_sequence_core_status"] == "PASS"].copy()
    if len(mapped) != 201 or len(exact) != 189:
        raise ValueError(f"Unexpected structure tiers: mapped={len(mapped)}, exact={len(exact)}")
    split_roles = binder_ready[["resource_record_id", "fold_id", "split_role", "binder_sequence"]]
    coordinate_root = OUTPUT / "coordinates/core_mapped_95"
    coordinate_rows = []
    for row in mapped.to_dict("records"):
        source = SOURCE_COORDINATES / row["structure_path"]
        destination = coordinate_root / f"{row['structure_instance_id']}.pdb"
        digest, size, atoms, chains = copy_coordinate(source, destination)
        coordinate_rows.append(
            {
                "structure_instance_id": row["structure_instance_id"],
                "resource_record_id": row["resource_record_id"],
                "pdb_code": row["pdb_code"],
                "resolved_chain_id": row["resolved_chain_id"],
                "coordinate_relpath": destination.relative_to(PACKAGE).as_posix(),
                "coordinate_sha256": digest,
                "coordinate_bytes": size,
                "parsed_atom_count": atoms,
                "parsed_chain_ids": chains,
                "structure_core_tier": "mapped_95",
            }
        )
    coordinate_manifest = pd.DataFrame(coordinate_rows).merge(
        split_roles, on="resource_record_id", how="left", validate="one_to_one"
    )
    assert_no_missing(coordinate_manifest, ["coordinate_relpath", "coordinate_sha256", "fold_id", "split_role"], "mapped coordinates")
    coordinate_summary = write_table(coordinate_manifest, "structure_core_mapped_95_coordinate_manifest_v1")
    mapped_ready = mapped.merge(coordinate_manifest, on=["structure_instance_id", "resource_record_id", "pdb_code", "resolved_chain_id"], how="inner", validate="one_to_one")
    mapped_summary = write_table(mapped_ready, "structure_core_mapped_95_ml_ready_8_1_1_v1")
    exact_ready = mapped_ready[mapped_ready["exact_sequence_core_status"] == "PASS"].copy()
    exact_ready["structure_core_tier"] = "exact_sequence"
    exact_summary = write_table(exact_ready, "structure_core_exact_sequence_ml_ready_8_1_1_v1")

    complex_binder = read_csv(COMPLEX / "data/complex_grade_binder_contact_candidates.csv")
    complex_window = read_csv(COMPLEX / "data/complex_grade_window_contact_candidates.csv")
    registry_by_accession = registry[["uniprot_accession", "target_sequence_id", "sequence", "sequence_sha256"]]

    def build_complex(frame: pd.DataFrame, tier: str) -> tuple[pd.DataFrame, dict[str, str | int | bool]]:
        ready = frame.merge(
            proteins[["resource_record_id", "sequence_object_id", "canonical_name", "scaffold_category", "final_sequence"]],
            on="resource_record_id",
            how="inner",
            validate="one_to_one",
        ).merge(
            registry_by_accession,
            left_on="best_target_accession",
            right_on="uniprot_accession",
            how="inner",
            validate="many_to_one",
        ).merge(
            split_roles[["resource_record_id", "fold_id", "split_role"]],
            on="resource_record_id",
            how="inner",
            validate="one_to_one",
        )
        if len(ready) != len(frame):
            raise ValueError(f"{tier}: target sequence or split join lost rows")
        paths = []
        for row in ready.to_dict("records"):
            source = ASSEMBLY_CACHE / row["pdb_code"] / f"{row['pdb_code']}.pdb{row['best_assembly_id']}.gz"
            destination = OUTPUT / "coordinates/complex_biological_assemblies" / f"{row['pdb_code']}.pdb{row['best_assembly_id']}.pdb.gz"
            digest, size, atoms, chains = copy_coordinate(source, destination)
            paths.append((destination.relative_to(PACKAGE).as_posix(), digest, size, atoms, chains))
        ready[["assembly_coordinate_relpath", "assembly_coordinate_sha256", "assembly_coordinate_bytes", "assembly_parsed_atom_count", "assembly_parsed_chain_ids"]] = pd.DataFrame(paths, index=ready.index)
        ready = ready.rename(columns={"final_sequence": "binder_sequence", "sequence": "target_sequence"})
        ready["complex_tier"] = tier
        ready["source_release"] = "structure_complex_biological_assembly_qc_v1"
        assert_no_missing(ready, ["binder_sequence", "target_sequence", "assembly_coordinate_relpath", "fold_id", "split_role"], tier)
        return ready, write_table(ready, f"complex_{tier}_ml_ready_8_1_1_v1")

    complex_binder_ready, complex_binder_summary = build_complex(complex_binder, "binder_target_contact")
    complex_window_ready, complex_window_summary = build_complex(complex_window, "recognition_window_target_contact")

    summary = {
        "release_name": "FBBP_ML_ready_release_v1",
        "binder": binder_summary,
        "target_pair": target_summary,
        "structure_mapped_95": mapped_summary,
        "structure_exact_sequence": exact_summary,
        "structure_coordinate_manifest": coordinate_summary,
        "complex_binder_contact": complex_binder_summary,
        "complex_window_contact": complex_window_summary,
        "coordinate_counts": {
            "mapped_95_pdb": int(len(coordinate_manifest)),
            "exact_sequence_subset": int(len(exact_ready)),
            "binder_contact_rows": int(len(complex_binder_ready)),
            "window_contact_rows": int(len(complex_window_ready)),
            "unique_complex_assembly_files": int(complex_binder_ready["assembly_coordinate_relpath"].nunique()),
        },
        "task_boundary": {
            "binder": "supervised scaffold-family classification is available; binding/non-binding labels are not supplied",
            "target_pair": "positive curated target associations are available; no validated negative pairs are supplied",
            "structure": "coordinates are self-contained for the mapped-95 and exact-sequence tiers",
            "complex": "biological assemblies are self-contained for the 53/19 contact tiers; the 19-record tier is validation-only",
        },
    }
    (OUTPUT / "ML_READY_BUILD_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
