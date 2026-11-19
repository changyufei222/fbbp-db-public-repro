"""Minimal loader and integrity checks for FBBP ML-ready release v1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parent


def _read(name: str) -> pd.DataFrame:
    csv_path = ROOT / f"{name}.csv"
    parquet_path = ROOT / f"{name}.parquet"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    return pd.read_csv(csv_path, dtype=str, keep_default_na=False)


def load_binder() -> pd.DataFrame:
    return _read("binder_full_length_ml_ready_8_1_1_v1")


def load_target_pairs() -> pd.DataFrame:
    return _read("target_sequence_aware_pair_ml_ready_8_1_1_v1")


def load_structure_core(tier: str = "mapped_95") -> pd.DataFrame:
    names = {
        "mapped_95": "structure_core_mapped_95_ml_ready_8_1_1_v1",
        "exact_sequence": "structure_core_exact_sequence_ml_ready_8_1_1_v1",
    }
    if tier not in names:
        raise ValueError(f"Unknown structure tier: {tier}")
    return _read(names[tier])


def load_complex(tier: str = "recognition_window_target_contact") -> pd.DataFrame:
    names = {
        "binder_target_contact": "complex_binder_target_contact_ml_ready_8_1_1_v1",
        "recognition_window_target_contact": "complex_recognition_window_target_contact_ml_ready_8_1_1_v1",
    }
    if tier not in names:
        raise ValueError(f"Unknown complex tier: {tier}")
    return _read(names[tier])


def package_path(relative_path: str) -> Path:
    path = PACKAGE / relative_path
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def coordinate_path(relative_path: str) -> Path:
    return package_path(relative_path)


def target_fasta_path(relative_path: str) -> Path:
    return package_path(relative_path)


def smoke_test() -> dict[str, int | str]:
    binder = load_binder()
    pairs = load_target_pairs()
    mapped = load_structure_core("mapped_95")
    exact = load_structure_core("exact_sequence")
    complex_binder = load_complex("binder_target_contact")
    complex_window = load_complex("recognition_window_target_contact")
    checks = {
        "binder_rows": len(binder),
        "target_pair_rows": len(pairs),
        "target_pair_unique_binders": pairs["resource_record_id"].nunique(),
        "mapped_structure_rows": len(mapped),
        "exact_structure_rows": len(exact),
        "binder_contact_rows": len(complex_binder),
        "window_contact_rows": len(complex_window),
        "missing_core_coordinates": sum(not coordinate_path(path).exists() for path in mapped["coordinate_relpath"]),
        "missing_complex_coordinates": sum(
            not coordinate_path(path).exists() for path in complex_binder["assembly_coordinate_relpath"]
        ),
        "missing_target_fastas": sum(
            not target_fasta_path(path).exists() for path in pairs["fasta_file"]
        ),
    }
    expected = {
        "binder_rows": 1996,
        "target_pair_rows": 710,
        "target_pair_unique_binders": 679,
        "mapped_structure_rows": 201,
        "exact_structure_rows": 189,
        "binder_contact_rows": 53,
        "window_contact_rows": 19,
        "missing_core_coordinates": 0,
        "missing_complex_coordinates": 0,
        "missing_target_fastas": 0,
    }
    for key, value in expected.items():
        if checks[key] != value:
            raise AssertionError(f"{key}: expected {value}, observed {checks[key]}")
    checks["status"] = "PASS"
    return checks


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test:
        print(json.dumps(smoke_test(), indent=2))
    else:
        print("Import this module, then call load_binder(), load_target_pairs(), load_structure_core(), or load_complex().")
