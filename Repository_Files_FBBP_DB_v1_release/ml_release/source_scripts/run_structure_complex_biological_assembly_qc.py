#!/usr/bin/env python
"""Validate FBBP complex candidates using RCSB biological assemblies and loop contacts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
from datetime import date
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
STRUCTURE_QC = ROOT / "FBBP_ML_structure_core_qc_audited_v1_20260716/data/structure_core_coordinate_qc.csv"
TARGET_RELEASE = ROOT / "FBBP_ML_target_sequence_aware_10fold_v1"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_structure_complex_biological_assembly_qc_v1_20260716"
DEFAULT_CACHE = ROOT / "FBBP_RCSB_BIOASSEMBLY_CACHE_20260716"
CONTACT_DISTANCE = 5.0
STANDARD_RESIDUES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
USER_AGENT = "FBBP-DB-ML/20260716"

ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.match_score = 2.0
ALIGNER.mismatch_score = -1.0
ALIGNER.open_gap_score = -8.0
ALIGNER.extend_gap_score = -0.5
ALIGNER.open_left_insertion_score = 0.0
ALIGNER.extend_left_insertion_score = 0.0
ALIGNER.open_right_insertion_score = 0.0
ALIGNER.extend_right_insertion_score = 0.0
ALIGNER.open_left_deletion_score = 0.0
ALIGNER.extend_left_deletion_score = 0.0
ALIGNER.open_right_deletion_score = 0.0
ALIGNER.extend_right_deletion_score = 0.0


def fetch_bytes(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=120) as response:
        return response.read()


def entry_assembly_count(pdb_code: str, cache_dir: Path) -> int:
    metadata_path = cache_dir / pdb_code / "entry.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    if not metadata_path.exists():
        raw = fetch_bytes(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_code}")
        metadata_path.write_bytes(raw)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return max(1, int(metadata.get("rcsb_entry_info", {}).get("assembly_count", 1) or 1))


def assembly_path(pdb_code: str, assembly_id: int, cache_dir: Path) -> Path | None:
    path = cache_dir / pdb_code / f"{pdb_code}.pdb{assembly_id}.gz"
    if path.exists():
        return path
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb{assembly_id}.gz"
    try:
        path.write_bytes(fetch_bytes(url))
    except HTTPError as error:
        if error.code == 404:
            return None
        raise
    return path


def alignment_metrics(expected: str, observed: str) -> dict:
    if not expected or not observed:
        return {
            "identity": 0.0,
            "expected_coverage": 0.0,
            "observed_coverage": 0.0,
            "aligned_count": 0,
            "expected_to_observed": {},
        }
    alignment = ALIGNER.align(expected, observed)[0]
    expected_blocks, observed_blocks = alignment.aligned
    pairs = []
    for (expected_start, expected_end), (observed_start, observed_end) in zip(expected_blocks, observed_blocks):
        length = min(int(expected_end - expected_start), int(observed_end - observed_start))
        pairs.extend(
            (int(expected_start) + offset, int(observed_start) + offset)
            for offset in range(length)
        )
    if not pairs:
        return alignment_metrics("", "")
    exact = sum(expected[left] == observed[right] for left, right in pairs)
    return {
        "identity": exact / len(pairs),
        "expected_coverage": len({left for left, _ in pairs}) / len(expected),
        "observed_coverage": len({right for _, right in pairs}) / len(observed),
        "aligned_count": len(pairs),
        "expected_to_observed": {left: right for left, right in pairs},
    }


def chain_record(chain) -> dict:
    residues = []
    for residue in chain.get_residues():
        residue_name = residue.get_resname().strip().upper()
        if residue.id[0] == " ":
            pass
        elif not is_aa(residue, standard=False):
            continue
        atoms = [
            atom
            for atom in residue.get_atoms()
            if str(atom.element).upper() != "H" and atom.get_occupancy() != 0
        ]
        if not atoms:
            continue
        residues.append(
            {
                "residue_name": residue_name,
                "one_letter": STANDARD_RESIDUES.get(residue_name, "X"),
                "atoms": atoms,
            }
        )
    return {
        "sequence": "".join(residue["one_letter"] for residue in residues),
        "residues": residues,
    }


def coordinates(atoms: list) -> np.ndarray:
    if not atoms:
        return np.empty((0, 3), dtype=float)
    return np.asarray([atom.coord for atom in atoms], dtype=float)


def minimum_distance(left: np.ndarray, right: np.ndarray) -> float:
    if not len(left) or not len(right):
        return float("inf")
    best = float("inf")
    for start in range(0, len(left), 128):
        block = left[start:start + 128]
        distances = np.linalg.norm(block[:, None, :] - right[None, :, :], axis=2)
        best = min(best, float(distances.min()))
        if best <= CONTACT_DISTANCE:
            return best
    return best


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def evaluate_structure(
    row: pd.Series,
    target_sequences: list[dict],
    recognition_windows: list[dict],
    cache_dir: Path,
) -> dict:
    result = {
        "structure_instance_id": row["structure_instance_id"],
        "resource_record_id": row["resource_record_id"],
        "pdb_code": row["pdb_code"],
        "file_core_status": row["status"],
        "exact_sequence_core_status": row["exact_sequence_core_status"],
        "resolved_target_sequence_count": len(target_sequences),
        "recognition_window_count": len(recognition_windows),
        "assembly_count_reported": 0,
        "assemblies_downloaded": 0,
        "best_assembly_id": "",
        "best_model_id": "",
        "best_binder_chain_id": "",
        "best_target_chain_id": "",
        "best_target_accession": "",
        "binder_chain_identity": 0.0,
        "binder_expected_coverage": 0.0,
        "target_chain_identity": 0.0,
        "target_observed_chain_coverage": 0.0,
        "mapped_recognition_window_count": 0,
        "minimum_window_target_distance_angstrom": "",
        "minimum_binder_target_distance_angstrom": "",
        "complex_grade_status": "FAIL",
        "failure_reason": "",
    }
    if not target_sequences:
        result["failure_reason"] = "no_frozen_target_sequence"
        return result
    pdb_code = str(row["pdb_code"]).upper()
    try:
        assembly_count = entry_assembly_count(pdb_code, cache_dir)
    except Exception as error:
        result["failure_reason"] = f"entry_metadata_failed:{type(error).__name__}"
        return result
    result["assembly_count_reported"] = assembly_count
    best_evidence = None
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    for assembly_id in range(1, assembly_count + 1):
        try:
            path = assembly_path(pdb_code, assembly_id, cache_dir)
        except Exception:
            continue
        if path is None:
            continue
        result["assemblies_downloaded"] += 1
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                structure = parser.get_structure(f"{pdb_code}_{assembly_id}", handle)
        except Exception:
            continue
        for model in structure:
            chains = []
            for chain in model:
                record = chain_record(chain)
                if record["sequence"]:
                    chains.append({"chain_id": str(chain.id), **record})
            binder_candidates = []
            for chain in chains:
                metrics = alignment_metrics(row["final_sequence"], chain["sequence"])
                if metrics["identity"] >= 0.95 and metrics["expected_coverage"] >= 0.80 and metrics["observed_coverage"] >= 0.80:
                    binder_candidates.append({**chain, "metrics": metrics})
            if not binder_candidates:
                continue
            target_candidates = []
            for chain in chains:
                for target in target_sequences:
                    metrics = alignment_metrics(target["sequence"], chain["sequence"])
                    minimum_aligned = min(20, len(chain["sequence"]), len(target["sequence"]))
                    if metrics["identity"] >= 0.90 and metrics["observed_coverage"] >= 0.80 and metrics["aligned_count"] >= minimum_aligned:
                        target_candidates.append({**chain, "metrics": metrics, **target})
            for binder in binder_candidates:
                binder_atoms = [
                    atom
                    for observed_index in binder["metrics"]["expected_to_observed"].values()
                    for atom in binder["residues"][observed_index]["atoms"]
                ]
                binder_coords = coordinates(binder_atoms)
                window_atoms = []
                mapped_window_count = 0
                mapping = binder["metrics"]["expected_to_observed"]
                for window in recognition_windows:
                    start = int(window["canonical_loop_start"]) - 1
                    end = int(window["canonical_loop_end"])
                    expected_positions = list(range(start, end))
                    observed_positions = [mapping[position] for position in expected_positions if position in mapping]
                    if expected_positions and len(observed_positions) / len(expected_positions) >= 0.80:
                        mapped_window_count += 1
                        window_atoms.extend(
                            atom
                            for observed_index in observed_positions
                            for atom in binder["residues"][observed_index]["atoms"]
                        )
                window_coords = coordinates(window_atoms)
                for target in target_candidates:
                    if target["chain_id"] == binder["chain_id"]:
                        continue
                    target_coords = coordinates(
                        [atom for residue in target["residues"] for atom in residue["atoms"]]
                    )
                    binder_distance = minimum_distance(binder_coords, target_coords)
                    window_distance = minimum_distance(window_coords, target_coords)
                    evidence = {
                        "assembly_id": assembly_id,
                        "model_id": str(model.id),
                        "binder_chain_id": binder["chain_id"],
                        "target_chain_id": target["chain_id"],
                        "target_accession": target["uniprot_accession"],
                        "binder_identity": binder["metrics"]["identity"],
                        "binder_coverage": binder["metrics"]["expected_coverage"],
                        "target_identity": target["metrics"]["identity"],
                        "target_coverage": target["metrics"]["observed_coverage"],
                        "mapped_window_count": mapped_window_count,
                        "window_distance": window_distance,
                        "binder_distance": binder_distance,
                    }
                    rank = (
                        window_distance > CONTACT_DISTANCE,
                        binder_distance > CONTACT_DISTANCE,
                        window_distance,
                        binder_distance,
                        -target["metrics"]["identity"],
                    )
                    if best_evidence is None or rank < best_evidence[0]:
                        best_evidence = (rank, evidence)
    if best_evidence is None:
        result["failure_reason"] = "no_sequence_resolved_binder_target_chain_pair_in_deposited_assembly"
        return result
    evidence = best_evidence[1]
    result.update(
        {
            "best_assembly_id": evidence["assembly_id"],
            "best_model_id": evidence["model_id"],
            "best_binder_chain_id": evidence["binder_chain_id"],
            "best_target_chain_id": evidence["target_chain_id"],
            "best_target_accession": evidence["target_accession"],
            "binder_chain_identity": round(evidence["binder_identity"], 6),
            "binder_expected_coverage": round(evidence["binder_coverage"], 6),
            "target_chain_identity": round(evidence["target_identity"], 6),
            "target_observed_chain_coverage": round(evidence["target_coverage"], 6),
            "mapped_recognition_window_count": evidence["mapped_window_count"],
            "minimum_window_target_distance_angstrom": round(evidence["window_distance"], 3) if np.isfinite(evidence["window_distance"]) else "",
            "minimum_binder_target_distance_angstrom": round(evidence["binder_distance"], 3) if np.isfinite(evidence["binder_distance"]) else "",
        }
    )
    if evidence["mapped_window_count"] and evidence["window_distance"] <= CONTACT_DISTANCE:
        result["complex_grade_status"] = "PASS_WINDOW_TARGET_CONTACT"
    elif evidence["binder_distance"] <= CONTACT_DISTANCE:
        result["complex_grade_status"] = "PASS_BINDER_TARGET_CONTACT_NO_WINDOW_CONTACT"
        result["failure_reason"] = "no_mapped_recognition_window_within_5A"
    else:
        result["failure_reason"] = "sequence_resolved_chains_without_5A_contact"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    for directory in (args.output / "data", args.output / "reports", args.output / "metadata", args.cache):
        directory.mkdir(parents=True, exist_ok=True)
    progress_path = args.output / "data/complex_structure_biological_assembly_qc.csv"
    if args.resume and progress_path.exists():
        results = pd.read_csv(progress_path, dtype=str).fillna("").to_dict("records")
    else:
        results = []
    completed = {row["structure_instance_id"] for row in results}
    qc = pd.read_csv(STRUCTURE_QC, dtype=str).fillna("")
    qc = qc[qc["status"] == "PASS"].copy()
    proteins = pd.read_csv(FOUNDATION / "data/entities/protein_objects.csv", dtype=str).fillna("")
    qc = qc.merge(
        proteins[["resource_record_id", "final_sequence"]],
        on="resource_record_id",
        how="left",
        validate="one_to_one",
    )
    links = pd.read_csv(TARGET_RELEASE / "data/protein_target_sequence_links.csv", dtype=str).fillna("")
    registry = pd.read_csv(TARGET_RELEASE / "data/frozen_target_sequence_registry.csv", dtype=str).fillna("")
    links = links.merge(
        registry[["target_sequence_id", "uniprot_accession", "sequence"]],
        on=["target_sequence_id", "uniprot_accession"],
        how="left",
        validate="many_to_one",
    )
    target_for_parent = {
        resource_record_id: group[["uniprot_accession", "target_sequence_id", "sequence"]].drop_duplicates().to_dict("records")
        for resource_record_id, group in links.groupby("resource_record_id")
    }
    windows = pd.read_csv(FOUNDATION / "data/entities/recognition_window_objects.csv", dtype=str).fillna("")
    windows = windows[
        (windows["canonical_loop_sequence"].str.strip() != "")
        & (windows["canonical_loop_start"].str.strip() != "")
        & (windows["canonical_loop_end"].str.strip() != "")
    ]
    windows_for_parent = {
        resource_record_id: group[["loop_id", "canonical_loop_start", "canonical_loop_end", "canonical_loop_sequence"]].to_dict("records")
        for resource_record_id, group in windows.groupby("resource_record_id")
    }
    for row in qc.sort_values("structure_instance_id").to_dict("records"):
        if row["structure_instance_id"] in completed:
            continue
        result = evaluate_structure(
            pd.Series(row),
            target_for_parent.get(row["resource_record_id"], []),
            windows_for_parent.get(row["resource_record_id"], []),
            args.cache,
        )
        results.append(result)
        pd.DataFrame(results).sort_values("structure_instance_id").to_csv(progress_path, index=False)
        time.sleep(0.05)
    frame = pd.DataFrame(results).sort_values("structure_instance_id")
    frame.to_csv(progress_path, index=False)
    window_pass = frame[frame["complex_grade_status"] == "PASS_WINDOW_TARGET_CONTACT"].copy()
    binder_pass = frame[
        frame["complex_grade_status"].isin(
            ["PASS_WINDOW_TARGET_CONTACT", "PASS_BINDER_TARGET_CONTACT_NO_WINDOW_CONTACT"]
        )
    ].copy()
    window_pass.to_csv(args.output / "data/complex_grade_window_contact_candidates.csv", index=False)
    binder_pass.to_csv(args.output / "data/complex_grade_binder_contact_candidates.csv", index=False)
    audit = {
        "status": "PASS_WITH_COMPLEX_CANDIDATE_CLASSIFICATION",
        "release_date": str(date.today()),
        "mapped_core_structures_evaluated": int(len(frame)),
        "structures_with_frozen_target_sequence": int((frame["resolved_target_sequence_count"].astype(int) > 0).sum()),
        "window_target_contact_pass": int(len(window_pass)),
        "whole_binder_target_contact_pass_including_window_pass": int(len(binder_pass)),
        "status_counts": frame["complex_grade_status"].value_counts().to_dict(),
        "failure_reason_counts": frame["failure_reason"].value_counts().to_dict(),
        "rules": {
            "assembly_source": "all deposited RCSB biological assemblies",
            "binder_chain": ">=95% aligned identity, >=80% expected and observed chain coverage",
            "target_chain": ">=90% aligned identity, >=80% observed chain coverage, >=20 aligned residues when available",
            "recognition_window_mapping": ">=80% window residues mapped to the assembly binder chain",
            "contact": "minimum non-hydrogen atom distance <=5.0 Angstrom",
            "primary_complex_grade": "resolved target sequence plus recognition-window-to-target contact",
        },
        "not_created": ["benchmark_zip"],
    }
    (args.output / "reports/structure_complex_biological_assembly_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    (args.output / "README.md").write_text(
        "# FBBP biological-assembly complex QC\n\n"
        "This layer evaluates mapped-core experimental structures against all deposited RCSB biological "
        "assemblies. Binder and target chains are reidentified by sequence; the primary complex-grade tier "
        "requires a frozen target sequence and a recognition-window heavy-atom contact within 5 A.\n",
        encoding="utf-8",
    )
    checksum_path = args.output / "metadata/checksums_sha256.json"
    files = sorted(path for path in args.output.rglob("*") if path.is_file() and path != checksum_path)
    checksum_path.write_text(
        json.dumps(
            [{"path": path.relative_to(args.output).as_posix(), "sha256": sha256(path)} for path in files],
            indent=2,
        ),
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
