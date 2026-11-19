#!/usr/bin/env python
"""Build the FBBP full-length binder-only 10-fold assignment from audited clusters.

This script intentionally does not create a train/validation/test export or a ZIP.
It assigns whole 95% within-family full-length binder clusters to frozen fold IDs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections import defaultdict
from datetime import date
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner


RELEASE_NAME = "binder_full_length_10fold_v1"
DEFAULT_SEED = "20260714"
FOLD_COUNT = 10
ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_binder_full_length_10fold_v1"


class UnionFind:
    def __init__(self, items: list[str]):
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in items}

    def find(self, item: str) -> str:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: str, right: str) -> None:
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if self.size[left] < self.size[right]:
            left, right = right, left
        self.parent[right] = left
        self.size[left] += self.size[right]

    def groups(self) -> dict[str, list[str]]:
        grouped: dict[str, list[str]] = defaultdict(list)
        for item in self.parent:
            grouped[self.find(item)].append(item)
        return grouped


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def sequence_identity(left: str, right: str, aligner: PairwiseAligner) -> float:
    longer = max(len(left), len(right))
    shorter = min(len(left), len(right))
    if not longer or shorter / longer < 0.95:
        return 0.0
    if left == right:
        return 1.0
    alignment = aligner.align(left, right)[0]
    left_blocks, right_blocks = alignment.aligned
    matches = 0
    for (left_start, left_end), (right_start, right_end) in zip(left_blocks, right_blocks):
        matches += sum(
            a == b
            for a, b in zip(
                left[int(left_start):int(left_end)],
                right[int(right_start):int(right_end)],
            )
        )
    return matches / longer


def stable_rank(seed: str, cluster_id: str, fold_id: int) -> int:
    value = f"{seed}|{cluster_id}|{fold_id}".encode("utf-8")
    return int(hashlib.sha256(value).hexdigest()[:16], 16)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_inputs() -> pd.DataFrame:
    proteins = pd.read_csv(FOUNDATION / "data/entities/protein_objects.csv").fillna("")
    if proteins["resource_record_id"].duplicated().any():
        raise ValueError("Protein IDs are not unique")
    required = [
        "resource_record_id",
        "protein_row_id",
        "sequence_object_id",
        "scaffold_category",
        "final_sequence",
    ]
    if proteins[required].eq("").any().any():
        raise ValueError("Binder clustering inputs contain blank required values")
    return proteins[required].sort_values("resource_record_id").reset_index(drop=True)


def build_clusters(proteins: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    union_find = UnionFind(proteins["resource_record_id"].tolist())
    aligner = make_aligner()
    comparison_count = 0
    length_eligible_comparison_count = 0
    linked_pair_count = 0
    sequence_for_id = proteins.set_index("resource_record_id")["final_sequence"].to_dict()
    for _, family in proteins.groupby("scaffold_category", sort=True):
        members = family.sort_values("resource_record_id")["resource_record_id"].tolist()
        for index, left_id in enumerate(members):
            left = sequence_for_id[left_id]
            for right_id in members[index + 1:]:
                comparison_count += 1
                right = sequence_for_id[right_id]
                if min(len(left), len(right)) / max(len(left), len(right)) < 0.95:
                    continue
                length_eligible_comparison_count += 1
                if sequence_identity(left, right, aligner) >= 0.95:
                    union_find.union(left_id, right_id)
                    linked_pair_count += 1
    groups = sorted(
        (sorted(members) for members in union_find.groups().values()),
        key=lambda members: members[0],
    )
    cluster_for_id = {}
    cluster_size = {}
    for index, members in enumerate(groups, start=1):
        cluster_id = f"FBBPSEQ95_{index:05d}"
        for resource_record_id in members:
            cluster_for_id[resource_record_id] = cluster_id
            cluster_size[resource_record_id] = len(members)
    result = proteins.drop(columns=["final_sequence"]).copy()
    result["binder_cluster_95_id"] = result["resource_record_id"].map(cluster_for_id)
    result["binder_cluster_size"] = result["resource_record_id"].map(cluster_size).astype(int)
    result["identity_threshold"] = 0.95
    result["clustering_scope"] = "within_scaffold_family_full_length_sequence"
    result["clustering_method"] = "global_affine_alignment_exact_matches_divided_by_longer_sequence_length"
    cluster_family_count = result.groupby("binder_cluster_95_id")["scaffold_category"].nunique()
    if (cluster_family_count != 1).any():
        raise ValueError("A binder cluster spans multiple scaffold families")
    diagnostics = {
        "comparison_count": comparison_count,
        "length_eligible_comparison_count": length_eligible_comparison_count,
        "linked_pair_count": linked_pair_count,
        "cluster_count": len(groups),
        "largest_cluster": max(len(members) for members in groups),
    }
    return result.sort_values("resource_record_id").reset_index(drop=True), diagnostics


def cluster_table(records: pd.DataFrame) -> pd.DataFrame:
    summary = (
        records.groupby(["binder_cluster_95_id", "scaffold_category"], as_index=False)
        .agg(cluster_member_count=("resource_record_id", "size"))
        .sort_values(
            ["cluster_member_count", "scaffold_category", "binder_cluster_95_id"],
            ascending=[False, True, True],
        )
        .reset_index(drop=True)
    )
    return summary


def assign_folds(clusters: pd.DataFrame, seed: str) -> pd.DataFrame:
    total_records = int(clusters["cluster_member_count"].sum())
    family_totals = clusters.groupby("scaffold_category")["cluster_member_count"].sum()
    fold_totals = [0] * FOLD_COUNT
    fold_family_totals: dict[str, list[int]] = {
        family: [0] * FOLD_COUNT for family in family_totals.index
    }
    assignments: dict[str, int] = {}

    for row in clusters.itertuples(index=False):
        cluster_id = row.binder_cluster_95_id
        family = row.scaffold_category
        size = int(row.cluster_member_count)
        candidates = []
        for fold_id in range(FOLD_COUNT):
            # Allocate every family across folds before placing a second cluster
            # into a fold for that family. Global total then breaks equal family
            # states. This hardens the assignment against empty folds.
            candidates.append(
                (
                    fold_family_totals[family][fold_id],
                    fold_totals[fold_id],
                    stable_rank(seed, cluster_id, fold_id),
                    fold_id,
                )
            )
        fold_id = min(candidates)[-1]
        assignments[cluster_id] = fold_id
        fold_totals[fold_id] += size
        fold_family_totals[family][fold_id] += size

    assigned = clusters.copy()
    assigned["fold_id"] = assigned["binder_cluster_95_id"].map(assignments).astype(int)
    return assigned


def write_outputs(
    records: pd.DataFrame,
    clusters: pd.DataFrame,
    output: Path,
    seed: str,
    clustering_diagnostics: dict,
) -> None:
    if output.exists():
        shutil.rmtree(output)
    (output / "data").mkdir(parents=True)
    (output / "reports").mkdir(parents=True)
    (output / "metadata").mkdir(parents=True)

    fold_map = clusters.set_index("binder_cluster_95_id")["fold_id"]
    assignments = records.copy()
    assignments["fold_id"] = assignments["binder_cluster_95_id"].map(fold_map).astype(int)
    assignments["release_name"] = RELEASE_NAME
    assignments["split_role"] = "frozen_10fold"
    assignments = assignments[
        [
            "resource_record_id",
            "protein_row_id",
            "sequence_object_id",
            "scaffold_category",
            "binder_cluster_95_id",
            "binder_cluster_size",
            "clustering_method",
            "fold_id",
            "release_name",
            "split_role",
        ]
    ].sort_values("resource_record_id")
    assignments_path = output / "data/binder_full_length_10fold_v1.csv"
    assignments.to_csv(assignments_path, index=False)
    records.to_csv(output / "data/binder_sequence_clusters_95_affine_v1.csv", index=False)

    clusters_path = output / "data/binder_cluster_to_fold_v1.csv"
    clusters.sort_values("binder_cluster_95_id").to_csv(clusters_path, index=False)

    global_balance = (
        assignments.groupby("fold_id", as_index=False)
        .agg(
            protein_objects=("resource_record_id", "size"),
            binder_clusters=("binder_cluster_95_id", "nunique"),
            largest_cluster_in_fold=("binder_cluster_size", "max"),
        )
        .set_index("fold_id")
        .reindex(range(FOLD_COUNT), fill_value=0)
        .reset_index()
    )
    family_balance = (
        assignments.groupby(["scaffold_category", "fold_id"], as_index=False)
        .agg(
            protein_objects=("resource_record_id", "size"),
            binder_clusters=("binder_cluster_95_id", "nunique"),
        )
    )
    families = sorted(assignments["scaffold_category"].unique())
    full_index = pd.MultiIndex.from_product(
        [families, range(FOLD_COUNT)], names=["scaffold_category", "fold_id"]
    )
    family_balance = (
        family_balance.set_index(["scaffold_category", "fold_id"])
        .reindex(full_index, fill_value=0)
        .reset_index()
    )
    balance = family_balance.merge(global_balance, on="fold_id", how="left")
    balance_path = output / "reports/binder_full_length_10fold_v1_family_balance.csv"
    balance.to_csv(balance_path, index=False)
    global_balance.to_csv(
        output / "reports/binder_full_length_10fold_v1_global_balance.csv", index=False
    )

    cluster_crossings = assignments.groupby("binder_cluster_95_id")["fold_id"].nunique()
    family_count_range = (
        family_balance.groupby("scaffold_category")["protein_objects"].agg(["min", "max"])
    )
    audit = {
        "release_name": RELEASE_NAME,
        "status": "PASS",
        "release_date": str(date.today()),
        "seed": seed,
        "fold_count": FOLD_COUNT,
        "node_definition": "all protein objects with audited normalized full-length sequences",
        "cluster_definition": "within-scaffold-family full-length single-linkage clusters at global affine alignment identity >=0.95",
        "identity_definition": "exact aligned residue matches divided by the longer full-length sequence",
        "protein_objects": int(len(assignments)),
        "binder_clusters": int(assignments["binder_cluster_95_id"].nunique()),
        "largest_binder_cluster": int(assignments["binder_cluster_size"].max()),
        "clustering_diagnostics": clustering_diagnostics,
        "fold_ids": sorted(assignments["fold_id"].unique().tolist()),
        "binder_cluster_crossings": int((cluster_crossings > 1).sum()),
        "unassigned_protein_objects": int(assignments["fold_id"].isna().sum()),
        "duplicate_protein_ids": int(assignments["resource_record_id"].duplicated().sum()),
        "global_fold_sizes": {
            str(row.fold_id): int(row.protein_objects)
            for row in global_balance.itertuples(index=False)
        },
        "family_fold_size_ranges": {
            family: {"min": int(row["min"]), "max": int(row["max"])}
            for family, row in family_count_range.iterrows()
        },
        "not_created": [
            "target_sequence_aware_10fold_v1",
            "benchmark_zip",
            "experimental_structure_core_task",
        ],
        "derived_train_validation_test_view": "FBBP_ML_derived_split_views_20260716",
    }
    if audit["binder_cluster_crossings"] != 0:
        raise RuntimeError("Binder clusters cross folds")
    if audit["unassigned_protein_objects"] != 0:
        raise RuntimeError("Some protein objects are unassigned")
    if audit["duplicate_protein_ids"] != 0:
        raise RuntimeError("Duplicate protein objects in assignment")
    if audit["fold_ids"] != list(range(FOLD_COUNT)):
        raise RuntimeError("Not all 10 folds are populated")

    audit_path = output / "reports/binder_full_length_10fold_v1_leakage_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = [
        "# FBBP binder-only full-length 10-fold audit",
        "",
        "Status: **PASS - frozen fold IDs; no train/validation/test export or ZIP**",
        "",
        f"- Protein objects: {audit['protein_objects']}",
        f"- 95% within-family binder clusters: {audit['binder_clusters']}",
        f"- Largest cluster: {audit['largest_binder_cluster']}",
        f"- Seed: {seed}",
        f"- Fold sizes: {audit['global_fold_sizes']}",
        f"- Binder-cluster crossings: {audit['binder_cluster_crossings']}",
        "- Recognition windows are not clustered and inherit their parent protein fold.",
        "- Target-aware assignment is not included: target sequences are not yet frozen or clustered.",
        "- Structure-core assignment is not included: X-ray/EM candidates still require coordinate QC.",
    ]
    (output / "reports/binder_full_length_10fold_v1_audit.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    readme = [
        "# FBBP binder-only full-length 10-fold release",
        "",
        "This directory contains a deterministic fold assignment for the full-length FBBP binder-only task.",
        "It does not contain a target-sequence-aware split, a train/validation/test view, or a benchmark ZIP.",
        "All recognition-window objects inherit the fold of their parent protein.",
    ]
    (output / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    manifest_files = sorted(path for path in output.rglob("*") if path.is_file())
    manifest = [
        {"path": path.relative_to(output).as_posix(), "sha256": sha256(path)}
        for path in manifest_files
    ]
    (output / "metadata/checksums_sha256.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", default=DEFAULT_SEED)
    args = parser.parse_args()
    proteins = load_inputs()
    records, clustering_diagnostics = build_clusters(proteins)
    clusters = assign_folds(cluster_table(records), args.seed)
    write_outputs(records, clusters, args.output, args.seed, clustering_diagnostics)
    print(args.output)


if __name__ == "__main__":
    main()
