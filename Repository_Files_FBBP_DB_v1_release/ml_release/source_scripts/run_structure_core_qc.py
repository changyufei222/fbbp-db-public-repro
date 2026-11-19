#!/usr/bin/env python
"""Run conservative coordinate QC for FBBP X-ray/EM structure-core candidates."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections import OrderedDict
from datetime import date
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB.Polypeptide import is_aa


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
PACKAGE_ROOT = ROOT / "结构序列整理/final_sequence_structure_web_package"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_structure_core_qc_audited_v1_20260716"
STANDARD_RESIDUES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
REQUIRED_BACKBONE = {"N", "CA", "C"}

ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.match_score = 2.0
ALIGNER.mismatch_score = -1.0
ALIGNER.open_gap_score = -8.0
ALIGNER.extend_gap_score = -0.5
# Terminal overhangs are eligible for deterministic trimming; internal gaps are penalized.
ALIGNER.open_left_insertion_score = 0.0
ALIGNER.extend_left_insertion_score = 0.0
ALIGNER.open_right_insertion_score = 0.0
ALIGNER.extend_right_insertion_score = 0.0
ALIGNER.open_left_deletion_score = 0.0
ALIGNER.extend_left_deletion_score = 0.0
ALIGNER.open_right_deletion_score = 0.0
ALIGNER.extend_right_deletion_score = 0.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_chain(path: Path, chain_id: str) -> list[dict]:
    residues: OrderedDict[tuple[str, str, str], dict] = OrderedDict()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[21].strip() != str(chain_id).strip():
                continue
            atom = line[12:16].strip()
            residue = line[17:20].strip().upper()
            record_type = line[:6].strip()
            # Water, ions and ligands can share a chain ID. Retain only polymer
            # ATOM records and amino-acid-like HETATM records such as MSE.
            if record_type == "HETATM" and residue not in STANDARD_RESIDUES and not is_aa(residue, standard=False):
                continue
            key = (line[21], line[22:26].strip(), line[26].strip())
            if key not in residues:
                residues[key] = {"residue_names": set(), "atoms": set(), "record_types": set()}
            residues[key]["residue_names"].add(residue)
            residues[key]["atoms"].add(atom)
            residues[key]["record_types"].add(record_type)
    result = []
    for key, value in residues.items():
        residue_names = sorted(value["residue_names"])
        standard_names = [name for name in residue_names if name in STANDARD_RESIDUES]
        residue = standard_names[0] if len(standard_names) == 1 else residue_names[0]
        result.append(
            {
                "chain": key[0],
                "resseq": key[1],
                "icode": key[2],
                "residue": residue,
                "residue_names": ";".join(residue_names),
                "ambiguous_residue_name": len(residue_names) != 1,
                "atoms": value["atoms"],
                "record_types": value["record_types"],
            }
        )
    return result


def available_chains(path: Path) -> list[str]:
    chains = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                chains.add(line[21].strip())
    return sorted(chains)


def empty_sequence_metrics(expected: str, observed: str) -> dict:
    return {
        "aligned_pair_count": 0,
        "exact_match_count": 0,
        "aligned_sequence_identity": 0.0,
        "retained_core_coordinate_coverage": 0.0,
        "full_sequence_coordinate_coverage": 0.0,
        "internal_missing_residue_count": len(expected),
        "internal_observed_insertion_count": len(observed),
        "n_terminal_trimmed_residues": len(expected),
        "c_terminal_trimmed_residues": 0,
        "_observed_core_indices": set(),
    }


def sequence_metrics(expected: str, observed: str) -> dict:
    if not expected or not observed:
        return empty_sequence_metrics(expected, observed)
    alignment = ALIGNER.align(expected, observed)[0]
    expected_blocks, observed_blocks = alignment.aligned
    if len(expected_blocks) == 0:
        return empty_sequence_metrics(expected, observed)
    aligned_pairs: list[tuple[int, int]] = []
    for (expected_start, expected_end), (observed_start, observed_end) in zip(expected_blocks, observed_blocks):
        block_length = min(int(expected_end - expected_start), int(observed_end - observed_start))
        aligned_pairs.extend(
            (int(expected_start) + offset, int(observed_start) + offset)
            for offset in range(block_length)
        )
    if not aligned_pairs:
        return empty_sequence_metrics(expected, observed)
    expected_indices = {expected_index for expected_index, _ in aligned_pairs}
    observed_indices = {observed_index for _, observed_index in aligned_pairs}
    exact_matches = sum(
        expected[expected_index] == observed[observed_index]
        for expected_index, observed_index in aligned_pairs
    )
    expected_start = min(expected_indices)
    expected_end = max(expected_indices) + 1
    observed_start = min(observed_indices)
    observed_end = max(observed_indices) + 1
    expected_core_span = expected_end - expected_start
    observed_core_span = observed_end - observed_start
    aligned_count = len(aligned_pairs)
    return {
        "aligned_pair_count": aligned_count,
        "exact_match_count": exact_matches,
        "aligned_sequence_identity": round(exact_matches / aligned_count, 6),
        "retained_core_coordinate_coverage": round(aligned_count / expected_core_span, 6),
        "full_sequence_coordinate_coverage": round(len(expected_indices) / len(expected), 6),
        "internal_missing_residue_count": max(0, expected_core_span - len(expected_indices)),
        "internal_observed_insertion_count": max(0, observed_core_span - len(observed_indices)),
        "n_terminal_trimmed_residues": expected_start,
        "c_terminal_trimmed_residues": len(expected) - expected_end,
        "_observed_core_indices": observed_indices,
    }


def assess_chain(path: Path, expected_sequence: str, chain_id: str) -> dict:
    residues = parse_chain(path, chain_id)
    observed = "".join(STANDARD_RESIDUES.get(record["residue"], "X") for record in residues)
    metrics = sequence_metrics(expected_sequence, observed)
    core_indices = metrics.pop("_observed_core_indices")
    core_residues = [record for index, record in enumerate(residues) if index in core_indices]
    noncanonical = sum(
        record["residue"] not in STANDARD_RESIDUES or record["ambiguous_residue_name"]
        for record in core_residues
    )
    missing_backbone = sum(
        not REQUIRED_BACKBONE.issubset(record["atoms"])
        for record in core_residues
    )
    reasons = []
    if not residues:
        reasons.append("chain_has_no_polymer_coordinates")
    if noncanonical:
        reasons.append("noncanonical_residue_in_retained_core")
    if missing_backbone:
        reasons.append("missing_N_CA_C_backbone_atom")
    if metrics["aligned_sequence_identity"] < 0.95:
        reasons.append("aligned_sequence_identity_lt_0_95")
    if metrics["retained_core_coordinate_coverage"] < 0.95:
        reasons.append("retained_core_coordinate_coverage_lt_0_95")
    if metrics["full_sequence_coordinate_coverage"] < 0.80:
        reasons.append("full_sequence_coordinate_coverage_lt_0_80")
    if metrics["internal_missing_residue_count"]:
        reasons.append("internal_coordinate_gap")
    if metrics["internal_observed_insertion_count"]:
        reasons.append("internal_observed_insertion")
    mapped_status = "PASS" if not reasons else "FAIL"
    exact_status = "PASS" if mapped_status == "PASS" and metrics["aligned_sequence_identity"] == 1.0 else "FAIL"
    return {
        "resolved_chain_id": chain_id,
        "polymer_coordinate_residue_count": len(residues),
        "retained_core_residue_count": len(core_residues),
        "noncanonical_core_residue_count": noncanonical,
        "missing_backbone_core_residue_count": missing_backbone,
        **metrics,
        "status": mapped_status,
        "exact_sequence_core_status": exact_status,
        "fail_reasons": ";".join(reasons),
    }


def qc_row(row: pd.Series) -> dict:
    path = PACKAGE_ROOT / row["package_best_structure_path"]
    result = {
        "structure_instance_id": row["structure_instance_id"],
        "resource_record_id": row["resource_record_id"],
        "scaffold_category": row["scaffold_category"],
        "pdb_code": row["pdb_code"],
        "chain_id": row["chain_id"],
        "resolved_chain_id": "",
        "chain_assignment_status": "",
        "structure_method": row["structure_method"],
        "resolution_angstrom": row["resolution_angstrom"],
        "expected_sequence_length": len(row["final_sequence"]),
        "source_structure_match_identity": row["structure_match_identity"],
        "source_structure_match_coverage": row["structure_match_coverage"],
        "source_structure_match_chain_coverage": row["structure_match_chain_coverage"],
        "structure_path": row["package_best_structure_path"],
        "file_present": path.exists(),
        "polymer_coordinate_residue_count": 0,
        "retained_core_residue_count": 0,
        "noncanonical_core_residue_count": 0,
        "missing_backbone_core_residue_count": 0,
        "aligned_pair_count": 0,
        "exact_match_count": 0,
        "aligned_sequence_identity": 0.0,
        "retained_core_coordinate_coverage": 0.0,
        "full_sequence_coordinate_coverage": 0.0,
        "internal_missing_residue_count": 0,
        "internal_observed_insertion_count": 0,
        "n_terminal_trimmed_residues": 0,
        "c_terminal_trimmed_residues": 0,
        "status": "FAIL",
        "exact_sequence_core_status": "FAIL",
        "fail_reasons": "",
    }
    if not path.exists():
        result["fail_reasons"] = "missing_coordinate_file"
        return result
    declared_chain = str(row["chain_id"]).strip()
    candidate_chains = available_chains(path)
    if declared_chain and declared_chain not in candidate_chains:
        candidate_chains.append(declared_chain)
    assessments = [assess_chain(path, row["final_sequence"], chain) for chain in candidate_chains]
    if not assessments:
        result["fail_reasons"] = "no_chain_in_coordinate_file"
        return result
    assessments.sort(
        key=lambda item: (
            item["exact_sequence_core_status"] != "PASS",
            item["status"] != "PASS",
            -item["aligned_sequence_identity"],
            -item["retained_core_coordinate_coverage"],
            -item["full_sequence_coordinate_coverage"],
            -item["aligned_pair_count"],
            item["missing_backbone_core_residue_count"],
            item["resolved_chain_id"],
        )
    )
    chosen = assessments[0]
    result.update(chosen)
    if declared_chain and chosen["resolved_chain_id"] == declared_chain:
        result["chain_assignment_status"] = "declared_chain"
    elif declared_chain and chosen["status"] == "PASS":
        result["chain_assignment_status"] = "declared_chain_reassigned_by_sequence"
    elif declared_chain:
        result["chain_assignment_status"] = "declared_chain_unresolved_best_available"
    else:
        pass_count = sum(item["status"] == "PASS" for item in assessments)
        if pass_count == 0:
            result["chain_assignment_status"] = "sequence_recovery_failed_best_available"
        elif pass_count == 1:
            result["chain_assignment_status"] = "sequence_recovered_unique"
        else:
            result["chain_assignment_status"] = "sequence_recovered_multiple_equivalent"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        shutil.rmtree(args.output)
    (args.output / "data").mkdir(parents=True)
    (args.output / "reports").mkdir(parents=True)
    (args.output / "metadata").mkdir(parents=True)

    structures = pd.read_csv(FOUNDATION / "data/entities/structure_instances.csv").fillna("")
    proteins = pd.read_csv(FOUNDATION / "data/entities/protein_objects.csv").fillna("")
    candidates = structures[
        structures["ml_structure_grade"] == "xray_em_core_candidate_pending_coordinate_qc"
    ].merge(
        proteins[["resource_record_id", "scaffold_category", "final_sequence"]],
        on="resource_record_id",
        how="left",
        validate="one_to_one",
    )
    rows = [qc_row(row) for _, row in candidates.iterrows()]
    qc = pd.DataFrame(rows).sort_values("resource_record_id")
    qc.to_csv(args.output / "data/structure_core_coordinate_qc.csv", index=False)
    qc[qc["status"] == "PASS"].to_csv(
        args.output / "data/structure_core_mapped_95_candidates.csv", index=False
    )
    qc[qc["exact_sequence_core_status"] == "PASS"].to_csv(
        args.output / "data/structure_core_exact_sequence_candidates.csv", index=False
    )
    qc[qc["status"] == "FAIL"].to_csv(
        args.output / "data/structure_core_qc_exclusions.csv", index=False
    )
    method_summary = (
        qc.groupby(["structure_method", "status"], as_index=False)
        .agg(records=("resource_record_id", "size"))
    )
    method_summary.to_csv(args.output / "reports/structure_core_qc_method_summary.csv", index=False)
    reason_summary = (
        qc.assign(fail_reason=qc["fail_reasons"].str.split(";"))
        .explode("fail_reason")
    )
    reason_summary = reason_summary[reason_summary["fail_reason"] != ""]
    reason_counts = reason_summary.groupby("fail_reason", as_index=False).agg(
        records=("resource_record_id", "size")
    )
    reason_counts.to_csv(args.output / "reports/structure_core_qc_failure_reasons.csv", index=False)
    pass_count = int((qc["status"] == "PASS").sum())
    exact_pass_count = int((qc["exact_sequence_core_status"] == "PASS").sum())
    chain_assignment_counts = {
        str(key): int(value)
        for key, value in qc["chain_assignment_status"].value_counts().sort_index().items()
    }
    coverage_sensitivity = {
        str(threshold): int(
            (
                (qc["aligned_sequence_identity"] >= 0.95)
                & (qc["retained_core_coordinate_coverage"] >= 0.95)
                & (qc["full_sequence_coordinate_coverage"] >= threshold)
                & (qc["internal_missing_residue_count"] == 0)
                & (qc["internal_observed_insertion_count"] == 0)
                & (qc["noncanonical_core_residue_count"] == 0)
                & (qc["missing_backbone_core_residue_count"] == 0)
            ).sum()
        )
        for threshold in (0.50, 0.80, 0.95)
    }
    audit = {
        "status": "PASS_WITH_CANDIDATE_CLASSIFICATION",
        "release_date": str(date.today()),
        "candidate_records": int(len(qc)),
        "mapped_95_structure_core_pass_records": pass_count,
        "exact_sequence_structure_core_pass_records": exact_pass_count,
        "structure_core_fail_records": int(len(qc) - pass_count),
        "chain_assignment_counts": chain_assignment_counts,
        "full_sequence_coverage_sensitivity_counts": coverage_sensitivity,
        "rules": {
            "coordinate_source": "declared local package_best_structure_path",
            "chain_resolution": "declared chain retained when valid; otherwise all coordinate chains are sequence-audited and reassignment is explicit",
            "polymer_filter": "ATOM residues plus amino-acid-like HETATM records; water, ions and non-polymer ligands excluded",
            "residue_identity": "exact residue matches counted position by position; nonstandard residues fail only inside the retained core",
            "backbone_atoms": "N, CA and C required for every retained core residue",
            "aligned_sequence_identity_minimum": 0.95,
            "retained_core_coordinate_coverage_minimum": 0.95,
            "full_sequence_coordinate_coverage_minimum": 0.80,
            "internal_coordinate_gaps": "not permitted",
            "internal_observed_insertions": "not permitted",
            "terminal_overhangs": "free in alignment and reported as deterministic trimming",
            "exact_sequence_tier": "all mapped-core rules plus 100% aligned sequence identity",
        },
    }
    (args.output / "reports/structure_core_qc_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    report = [
        "# FBBP structure-core coordinate QC",
        "",
        "Status: **candidate classification only; a complex-grade benchmark still requires chain-role and biological-assembly QC**",
        "",
        f"- X-ray/EM candidates inspected: {len(qc)}",
        f"- Mapped-core QC pass (>=95% identity; >=80% full-sequence coverage): {pass_count}",
        f"- Exact-sequence core pass: {exact_pass_count}",
        f"- Coordinate QC fail: {len(qc) - pass_count}",
        f"- Declared chains reassigned by sequence evidence: {chain_assignment_counts.get('declared_chain_reassigned_by_sequence', 0)}",
        "- Passing this file-level QC does not validate binder/target roles or biological assembly.",
    ]
    (args.output / "reports/structure_core_qc_audit.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    manifest = []
    for path in sorted(args.output.rglob("*")):
        if path.is_file():
            manifest.append({"path": path.relative_to(args.output).as_posix(), "sha256": sha256(path)})
    (args.output / "metadata/checksums_sha256.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(args.output)


if __name__ == "__main__":
    main()
