#!/usr/bin/env python
"""Build a parent-level 10-fold diagnostic constrained by full-length and loop similarity."""

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


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
BINDER_RELEASE = ROOT / "FBBP_ML_binder_full_length_10fold_v1"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_binder_loop_aware_10fold_diagnostic_v1_20260716"
FOLD_COUNT = 10
SEED = "20260716-loop-aware"
THRESHOLDS = {"80": 0.80, "95": 0.95, "exact": 1.0}


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
        result: dict[str, list[str]] = defaultdict(list)
        for item in self.parent:
            result[self.find(item)].append(item)
        return result


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def sequence_identity(left: str, right: str, aligner: PairwiseAligner, threshold: float) -> float:
    longer = max(len(left), len(right))
    shorter = min(len(left), len(right))
    if not longer or shorter / longer < threshold:
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


def cluster_windows(windows: pd.DataFrame, threshold: float, label: str) -> tuple[pd.DataFrame, dict]:
    work = windows[windows["canonical_loop_sequence"].str.strip() != ""].copy()
    union_find = UnionFind(work["loop_id"].tolist())
    aligner = make_aligner()
    comparisons = 0
    eligible = 0
    linked = 0
    sequence_for_id = work.set_index("loop_id")["canonical_loop_sequence"].to_dict()
    for _, group in work.groupby(["family_name", "loop_slot"], sort=True):
        members = group.sort_values("loop_id")["loop_id"].tolist()
        for index, left_id in enumerate(members):
            left = sequence_for_id[left_id]
            for right_id in members[index + 1:]:
                comparisons += 1
                right = sequence_for_id[right_id]
                if min(len(left), len(right)) / max(len(left), len(right)) < threshold:
                    continue
                eligible += 1
                if sequence_identity(left, right, aligner, threshold) >= threshold:
                    union_find.union(left_id, right_id)
                    linked += 1
    groups = sorted(
        (sorted(members) for members in union_find.groups().values()),
        key=lambda members: members[0],
    )
    cluster_for_id = {}
    size_for_id = {}
    for index, members in enumerate(groups, start=1):
        cluster_id = f"FBBPLOOP{label}_{index:05d}"
        for loop_id in members:
            cluster_for_id[loop_id] = cluster_id
            size_for_id[loop_id] = len(members)
    result = work.copy()
    result[f"loop_cluster_{label}_id"] = result["loop_id"].map(cluster_for_id)
    result[f"loop_cluster_{label}_size"] = result["loop_id"].map(size_for_id).astype(int)
    result["identity_threshold"] = threshold
    result["clustering_method"] = "within_family_and_loop_slot_global_affine_exact_matches_divided_by_longer_window"
    stats = {
        "threshold": threshold,
        "cluster_count": len(groups),
        "largest_window_cluster": max(len(members) for members in groups),
        "comparison_count": comparisons,
        "length_eligible_comparison_count": eligible,
        "linked_pair_count": linked,
    }
    return result, stats


def parent_components(
    proteins: pd.DataFrame,
    binder_clusters: pd.DataFrame,
    window_clusters: pd.DataFrame,
    cluster_column: str,
    label: str,
) -> pd.DataFrame:
    union_find = UnionFind(proteins["resource_record_id"].tolist())
    for _, group in binder_clusters.groupby("binder_cluster_95_id"):
        members = sorted(group["resource_record_id"].unique())
        for member in members[1:]:
            union_find.union(members[0], member)
    for _, group in window_clusters.groupby(cluster_column):
        members = sorted(group["resource_record_id"].unique())
        for member in members[1:]:
            union_find.union(members[0], member)
    groups = sorted(
        (sorted(members) for members in union_find.groups().values()),
        key=lambda members: members[0],
    )
    component_for_id = {}
    size_for_id = {}
    for index, members in enumerate(groups, start=1):
        component_id = f"FBBPLOOPCOMP{label}_{index:05d}"
        for resource_record_id in members:
            component_for_id[resource_record_id] = component_id
            size_for_id[resource_record_id] = len(members)
    result = proteins[["resource_record_id", "scaffold_category"]].copy()
    result[f"loop_aware_component_{label}_id"] = result["resource_record_id"].map(component_for_id)
    result[f"loop_aware_component_{label}_size"] = result["resource_record_id"].map(size_for_id).astype(int)
    return result


def stable_rank(component_id: str, fold_id: int) -> int:
    return int(hashlib.sha256(f"{SEED}|{component_id}|{fold_id}".encode()).hexdigest()[:16], 16)


def assign_folds(components: pd.DataFrame, component_column: str, size_column: str) -> pd.DataFrame:
    summary = components.groupby([component_column, "scaffold_category"], as_index=False).agg(
        component_size=("resource_record_id", "size")
    )
    if summary.groupby(component_column)["scaffold_category"].nunique().gt(1).any():
        raise RuntimeError("A loop-aware component spans scaffold families")
    summary = summary.sort_values(
        ["component_size", "scaffold_category", component_column],
        ascending=[False, True, True],
    )
    fold_totals = [0] * FOLD_COUNT
    family_fold_totals: dict[str, list[int]] = defaultdict(lambda: [0] * FOLD_COUNT)
    mapping = {}
    for row in summary.itertuples(index=False):
        component_id = getattr(row, component_column)
        candidates = [
            (
                fold_totals[fold_id],
                family_fold_totals[row.scaffold_category][fold_id],
                stable_rank(component_id, fold_id),
                fold_id,
            )
            for fold_id in range(FOLD_COUNT)
        ]
        fold_id = min(candidates)[-1]
        mapping[component_id] = fold_id
        fold_totals[fold_id] += row.component_size
        family_fold_totals[row.scaffold_category][fold_id] += row.component_size
    result = components.copy()
    result["fold_id"] = result[component_column].map(mapping).astype(int)
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        shutil.rmtree(args.output)
    for directory in ("data", "reports", "metadata"):
        (args.output / directory).mkdir(parents=True, exist_ok=True)
    proteins = pd.read_csv(FOUNDATION / "data/entities/protein_objects.csv", dtype=str).fillna("")
    windows = pd.read_csv(FOUNDATION / "data/entities/recognition_window_objects.csv", dtype=str).fillna("")
    binder = pd.read_csv(BINDER_RELEASE / "data/binder_sequence_clusters_95_affine_v1.csv", dtype=str).fillna("")
    cluster_frames = {}
    component_frames = {}
    stats = {}
    for label, threshold in THRESHOLDS.items():
        cluster_frame, cluster_stats = cluster_windows(windows, threshold, label)
        component_frame = parent_components(
            proteins,
            binder,
            cluster_frame,
            f"loop_cluster_{label}_id",
            label,
        )
        cluster_frames[label] = cluster_frame
        component_frames[label] = component_frame
        stats[label] = {
            **cluster_stats,
            "parent_component_count": int(component_frame[f"loop_aware_component_{label}_id"].nunique()),
            "largest_parent_component": int(component_frame[f"loop_aware_component_{label}_size"].max()),
        }
        cluster_frame.to_csv(args.output / f"data/recognition_window_clusters_{label}.csv", index=False)
        component_frame.to_csv(args.output / f"data/parent_components_{label}.csv", index=False)
    primary = component_frames["95"]
    primary = primary.merge(
        binder[["resource_record_id", "binder_cluster_95_id", "binder_cluster_size"]],
        on="resource_record_id",
        how="left",
        validate="one_to_one",
    )
    primary = assign_folds(primary, "loop_aware_component_95_id", "loop_aware_component_95_size")
    primary["release_name"] = "binder_loop_aware_10fold_diagnostic_v1"
    primary.to_csv(args.output / "data/binder_loop_aware_10fold_diagnostic_v1.csv", index=False)
    windows_with_fold = windows.merge(
        primary[["resource_record_id", "loop_aware_component_95_id", "fold_id"]],
        on="resource_record_id",
        how="left",
        validate="many_to_one",
    )
    windows_with_fold.to_csv(args.output / "data/recognition_windows_with_parent_fold.csv", index=False)
    fold_sizes = primary.groupby("fold_id")["resource_record_id"].nunique().reindex(range(FOLD_COUNT), fill_value=0)
    component_crossings = primary.groupby("loop_aware_component_95_id")["fold_id"].nunique()
    binder_crossings = primary.groupby("binder_cluster_95_id")["fold_id"].nunique()
    window_crossings = windows_with_fold.merge(
        cluster_frames["95"][["loop_id", "loop_cluster_95_id"]],
        on="loop_id",
        how="inner",
        validate="one_to_one",
    ).groupby("loop_cluster_95_id")["fold_id"].nunique()
    audit = {
        "status": "PASS_DIAGNOSTIC",
        "release_date": str(date.today()),
        "protein_objects": int(len(primary)),
        "recognition_windows": int(len(windows)),
        "recognition_windows_with_sequence": int((windows["canonical_loop_sequence"].str.strip() != "").sum()),
        "recognition_windows_without_sequence_excluded_from_clustering": int((windows["canonical_loop_sequence"].str.strip() == "").sum()),
        "primary_threshold": 0.95,
        "sensitivity": stats,
        "fold_sizes": {str(index): int(value) for index, value in fold_sizes.items()},
        "component_crossings": int((component_crossings > 1).sum()),
        "binder_cluster_crossings": int((binder_crossings > 1).sum()),
        "loop_cluster_crossings": int((window_crossings > 1).sum()),
        "not_created": ["benchmark_zip"],
    }
    if any(audit[key] for key in ("component_crossings", "binder_cluster_crossings", "loop_cluster_crossings")):
        raise RuntimeError("A leakage-controlled group crosses loop-aware folds")
    (args.output / "reports/binder_loop_aware_10fold_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    (args.output / "README.md").write_text(
        "# FBBP binder-loop-aware 10-fold diagnostic\n\n"
        "The primary 95% graph connects parent proteins through either the audited full-length binder cluster "
        "or a recognition-window cluster within the same family and loop slot. Exact and 80% results are "
        "sensitivity diagnostics. Windows always inherit the final parent fold. No benchmark ZIP is created.\n",
        encoding="utf-8",
    )
    files = sorted(path for path in args.output.rglob("*") if path.is_file())
    (args.output / "metadata/checksums_sha256.json").write_text(
        json.dumps(
            [{"path": path.relative_to(args.output).as_posix(), "sha256": sha256(path)} for path in files],
            indent=2,
        ),
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
