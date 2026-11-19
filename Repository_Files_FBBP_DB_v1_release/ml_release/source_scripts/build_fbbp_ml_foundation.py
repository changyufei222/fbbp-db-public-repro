#!/usr/bin/env python
"""Build full-length FBBP ML object registries without freezing a final split."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio.Align import PairwiseAligner


BINDER_CLUSTER_THRESHOLD = 0.95


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def clean_sequence(value: object) -> str:
    return re.sub(r"[^A-Za-z]", "", str(value)).upper()


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def normalized_key(value: object) -> str:
    return clean_text(value).upper()


def sequence_object_id(sequence: str) -> str:
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()[:12].upper()
    return f"FBBPSEQ_{digest}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def structure_grade(row: pd.Series) -> str:
    layer = row["structure_source_layer"]
    method = row["structure_method"].upper()
    resolution = row["resolution_angstrom"]
    if layer == "predicted_local_AF2":
        return "predicted_context_only"
    if "NMR" in method:
        return "nmr_flexibility_candidate_pending_model_qc"
    if "X-RAY" in method or method == "EM":
        if pd.isna(resolution):
            return "experimental_excluded_missing_resolution"
        if float(resolution) <= 3.5:
            return "xray_em_core_candidate_pending_coordinate_qc"
        return "experimental_excluded_resolution_gt_3_5"
    return "experimental_excluded_unsupported_method"


def target_fields(row: pd.Series) -> pd.Series:
    variant_id = clean_text(row.get("target_variant_id", ""))
    uniprot = normalized_key(row.get("uniprot_id", ""))
    species = normalized_key(row.get("species_name", ""))
    name = (
        normalized_key(row.get("gene_name_species", ""))
        or normalized_key(row.get("gene_name", ""))
        or normalized_key(row.get("protein_name", ""))
    )
    if uniprot:
        return pd.Series(
            [f"UP:{uniprot}", "accession_resolved_sequence_pending", True]
        )
    if name:
        return pd.Series([f"NAME:{name}|SPECIES:{species}", "named_no_sequence", True])
    if variant_id:
        return pd.Series(["", "placeholder_excluded", False])
    return pd.Series(["", "unannotated", False])


class UnionFind:
    def __init__(self, items: Iterable[str]):
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in items}

    def find(self, item: str) -> str:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]

    def groups(self) -> dict[str, list[str]]:
        groups: dict[str, list[str]] = defaultdict(list)
        for item in self.parent:
            groups[self.find(item)].append(item)
        return groups


def sequence_identity(seq_a: str, seq_b: str, aligner: PairwiseAligner) -> float:
    if not seq_a or not seq_b:
        return 0.0
    longer = max(len(seq_a), len(seq_b))
    shorter = min(len(seq_a), len(seq_b))
    if shorter / longer < BINDER_CLUSTER_THRESHOLD:
        return 0.0
    if len(seq_a) == len(seq_b):
        return sum(a == b for a, b in zip(seq_a, seq_b)) / longer
    return float(aligner.score(seq_a, seq_b)) / longer


def build_binder_sequence_clusters(
    protein_objects: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    work = protein_objects[
        [
            "resource_record_id",
            "protein_row_id",
            "sequence_object_id",
            "scaffold_category",
            "final_sequence",
        ]
    ].copy()
    record_ids = work["resource_record_id"].tolist()
    union_find = UnionFind(record_ids)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    comparison_count = 0
    linked_pair_count = 0
    for _, family in work.groupby("scaffold_category", dropna=False):
        records = family.sort_values("resource_record_id").to_dict("records")
        for left_index, left in enumerate(records):
            for right in records[left_index + 1 :]:
                longer = max(len(left["final_sequence"]), len(right["final_sequence"]))
                shorter = min(len(left["final_sequence"]), len(right["final_sequence"]))
                if longer == 0 or shorter / longer < BINDER_CLUSTER_THRESHOLD:
                    continue
                comparison_count += 1
                if (
                    sequence_identity(
                        left["final_sequence"], right["final_sequence"], aligner
                    )
                    >= BINDER_CLUSTER_THRESHOLD
                ):
                    union_find.union(
                        left["resource_record_id"], right["resource_record_id"]
                    )
                    linked_pair_count += 1

    groups = union_find.groups()
    ordered_groups = sorted(
        groups.values(), key=lambda members: (min(members), len(members))
    )
    cluster_for_record: dict[str, str] = {}
    cluster_size: dict[str, int] = {}
    for index, members in enumerate(ordered_groups, start=1):
        cluster_id = f"FBBPSEQ95_{index:05d}"
        for record_id in members:
            cluster_for_record[record_id] = cluster_id
            cluster_size[record_id] = len(members)

    work["binder_cluster_95_id"] = work["resource_record_id"].map(
        cluster_for_record
    )
    work["binder_cluster_size"] = work["resource_record_id"].map(cluster_size)
    work["identity_threshold"] = BINDER_CLUSTER_THRESHOLD
    work["clustering_scope"] = "within_scaffold_family_full_length_sequence"
    stats = {
        "identity_threshold": BINDER_CLUSTER_THRESHOLD,
        "comparison_count": comparison_count,
        "linked_pair_count": linked_pair_count,
        "cluster_count": int(work["binder_cluster_95_id"].nunique()),
        "largest_cluster_size": int(work["binder_cluster_size"].max()),
        "implementation": (
            "Within-family single-linkage graph; equal-length pairs use positional "
            "identity, unequal-length pairs use a global match score divided by the "
            "longer sequence length, and pairs below the length-ratio threshold are skipped."
        ),
    }
    return work.drop(columns="final_sequence"), stats


def build_strict_target_key_components(
    binder_clusters: pd.DataFrame, target_annotations: pd.DataFrame
) -> tuple[pd.DataFrame, dict]:
    eligible = target_annotations[target_annotations["target_split_eligible"]].merge(
        binder_clusters[
            [
                "resource_record_id",
                "protein_row_id",
                "scaffold_category",
                "binder_cluster_95_id",
            ]
        ],
        on=["resource_record_id", "protein_row_id"],
        how="left",
        validate="one_to_one",
    )
    record_ids = eligible["resource_record_id"].tolist()
    union_find = UnionFind(record_ids)
    for _, group in eligible.groupby("binder_cluster_95_id"):
        members = group["resource_record_id"].tolist()
        for member in members[1:]:
            union_find.union(members[0], member)
    for _, group in eligible.groupby("strict_target_key"):
        members = group["resource_record_id"].tolist()
        for member in members[1:]:
            union_find.union(members[0], member)

    groups = union_find.groups()
    ordered_groups = sorted(
        groups.values(), key=lambda members: (min(members), len(members))
    )
    component_for_record: dict[str, str] = {}
    component_size: dict[str, int] = {}
    for index, members in enumerate(ordered_groups, start=1):
        component_id = f"FBBPTKEY95_{index:05d}"
        for record_id in members:
            component_for_record[record_id] = component_id
            component_size[record_id] = len(members)

    eligible["strict_target_key_component_95_id"] = eligible[
        "resource_record_id"
    ].map(component_for_record)
    eligible["component_size"] = eligible["resource_record_id"].map(component_size)
    eligible["diagnostic_status"] = (
        "provisional_exact_target_key_graph_not_target_sequence_clustered"
    )
    largest_component_size = int(eligible["component_size"].max())
    stats = {
        "status": "provisional_diagnostic_not_a_released_split",
        "eligible_records": int(len(eligible)),
        "component_count": int(
            eligible["strict_target_key_component_95_id"].nunique()
        ),
        "largest_component_size": largest_component_size,
        "largest_component_fraction": round(largest_component_size / len(eligible), 6),
        "binder_identity_threshold": BINDER_CLUSTER_THRESHOLD,
        "target_edge": "exact strict target key; target sequences are not yet frozen or clustered",
    }
    columns = [
        "resource_record_id",
        "protein_row_id",
        "sequence_object_id",
        "structure_unique_sequence_id",
        "scaffold_category",
        "binder_cluster_95_id",
        "strict_target_key",
        "target_annotation_grade",
        "strict_target_key_component_95_id",
        "component_size",
        "diagnostic_status",
    ]
    return eligible[columns], stats


def resolve_loop_coordinates(row: pd.Series, sequences: dict[str, str]) -> pd.Series:
    sequence = sequences.get(str(row["protein_row_id"]), "")
    start = int(row["loop_start"])
    end = int(row["loop_end"])
    length = int(row["loop_length"])
    domain_start = int(row["domain_start"])
    if end - start + 1 == length and 1 <= start <= end <= len(sequence):
        return pd.Series(["direct", start, end, sequence[start - 1 : end]])
    canonical_start = start - domain_start + 1
    canonical_end = end - domain_start + 1
    if (
        end - start + 1 == length
        and 1 <= canonical_start <= canonical_end <= len(sequence)
    ):
        return pd.Series(
            [
                "domain_offset",
                canonical_start,
                canonical_end,
                sequence[canonical_start - 1 : canonical_end],
            ]
        )
    return pd.Series(["unresolved", "", "", ""])


def build(source_root: Path, output_root: Path) -> dict:
    source_root = source_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    db = source_root / "数据库表"
    manifest_path = (
        source_root
        / "结构序列整理"
        / "final_sequence_structure_web_package"
        / "00_manifest"
        / "web_ready_final_manifest.csv"
    )
    authority_path = (
        source_root
        / "数据整理"
        / "主表核心"
        / "source_core_main_1996_preprediction.csv"
    )

    proteins = read_csv(db / "proteins.csv").rename(columns={"protein_id": "protein_row_id"})
    manifest = read_csv(manifest_path)
    authority = read_csv(authority_path)
    active = read_csv(db / "active_structure_context.csv")
    interactions = read_csv(db / "interactions.csv")
    variants = read_csv(db / "target_species_variants.csv")
    concepts = read_csv(db / "targets_conceptual.csv")
    loops = read_csv(db / "loop_flexibility_public_summary.csv")

    source_paths = [
        db / "proteins.csv",
        manifest_path,
        authority_path,
        db / "active_structure_context.csv",
        db / "interactions.csv",
        db / "target_species_variants.csv",
        db / "targets_conceptual.csv",
        db / "loop_flexibility_public_summary.csv",
    ]

    manifest["resource_record_id"] = (
        manifest["protein_label"].str.extract(r"^(R\d+)_", expand=False).fillna("")
    )
    proteins["final_sequence"] = proteins["final_tested_sequence"].map(clean_sequence)
    proteins["sequence_object_id"] = proteins["final_sequence"].map(sequence_object_id)

    protein_objects = proteins.merge(
        manifest[
            [
                "resource_record_id",
                "structure_unique_sequence_id",
                "package_sequence_path",
                "package_best_structure_path",
                "selected_structure_type",
                "structure_source_layer",
                "structure_method",
                "submitted_structure_available",
            ]
        ],
        on="structure_unique_sequence_id",
        how="left",
        validate="one_to_one",
    )
    protein_objects["snapshot_identity_status"] = "unique_sequence_in_20260714_snapshot"
    if not protein_objects["resource_record_id"].str.match(r"^R\d+$").all():
        raise ValueError("Every protein object must have a public R identifier")
    if protein_objects["resource_record_id"].duplicated().any():
        raise ValueError("resource_record_id must be unique")
    if protein_objects["sequence_object_id"].duplicated().any():
        raise ValueError("Current snapshot unexpectedly contains duplicate full-length sequences")
    protein_columns = [
        "resource_record_id",
        "protein_row_id",
        "structure_unique_sequence_id",
        "sequence_object_id",
        "canonical_name",
        "scaffold_category",
        "organism",
        "sequence_length",
        "final_sequence",
        "snapshot_identity_status",
    ]
    write_csv(protein_objects[protein_columns], output_root / "data/entities/protein_objects.csv")

    authority_structure = authority[
        [
            "structure_unique_sequence_id",
            "Structural_Info_resolution",
            "Sequence_PDB_ID",
            "pdb_code",
            "chain_id",
            "structure_match_identity",
            "structure_match_coverage",
            "structure_match_chain_coverage",
        ]
    ].rename(
        columns={
            "pdb_code": "authority_pdb_code",
            "chain_id": "authority_chain_id",
            "structure_match_identity": "authority_structure_match_identity",
            "structure_match_coverage": "authority_structure_match_coverage",
            "structure_match_chain_coverage": "authority_structure_match_chain_coverage",
        }
    )
    active_structure = active[
        [
            "protein_row_id",
            "structure_unique_sequence_id",
            "pdb_code",
            "chain_id",
            "structure_match_identity",
            "structure_match_coverage",
            "structure_match_chain_coverage",
        ]
    ]
    structure_instances = (
        protein_objects.merge(
            authority_structure,
            on="structure_unique_sequence_id",
            how="left",
            validate="one_to_one",
        )
        .merge(
            active_structure,
            on=["protein_row_id", "structure_unique_sequence_id"],
            how="left",
            validate="one_to_one",
        )
        .fillna("")
    )
    structure_instances["pdb_code"] = structure_instances["pdb_code"].where(
        structure_instances["pdb_code"].ne(""),
        structure_instances["authority_pdb_code"].where(
            structure_instances["authority_pdb_code"].ne(""),
            structure_instances["Sequence_PDB_ID"],
        ),
    )
    structure_instances["chain_id"] = structure_instances["chain_id"].where(
        structure_instances["chain_id"].ne(""), structure_instances["authority_chain_id"]
    )
    for field in [
        "structure_match_identity",
        "structure_match_coverage",
        "structure_match_chain_coverage",
    ]:
        structure_instances[field] = structure_instances[field].where(
            structure_instances[field].ne(""),
            structure_instances[f"authority_{field}"],
        )
    structure_instances["resolution_angstrom"] = pd.to_numeric(
        structure_instances["Structural_Info_resolution"], errors="coerce"
    )
    structure_instances["structure_instance_id"] = (
        "FBBPSTR_" + structure_instances["resource_record_id"] + "_01"
    )
    structure_instances["ml_structure_grade"] = structure_instances.apply(
        structure_grade, axis=1
    )
    structure_columns = [
        "structure_instance_id",
        "resource_record_id",
        "protein_row_id",
        "sequence_object_id",
        "structure_unique_sequence_id",
        "structure_source_layer",
        "selected_structure_type",
        "structure_method",
        "resolution_angstrom",
        "pdb_code",
        "chain_id",
        "structure_match_identity",
        "structure_match_coverage",
        "structure_match_chain_coverage",
        "package_best_structure_path",
        "ml_structure_grade",
    ]
    write_csv(
        structure_instances[structure_columns],
        output_root / "data/entities/structure_instances.csv",
    )

    target_annotations = (
        interactions.merge(variants, on="target_variant_id", how="left")
        .merge(concepts, on="target_concept_id", how="left")
        .fillna("")
        .merge(
            protein_objects[
                [
                    "resource_record_id",
                    "protein_row_id",
                    "sequence_object_id",
                    "structure_unique_sequence_id",
                ]
            ],
            on=["protein_row_id", "structure_unique_sequence_id"],
            how="left",
            validate="one_to_one",
        )
    )
    target_annotations[
        ["strict_target_key", "target_annotation_grade", "target_split_eligible"]
    ] = target_annotations.apply(target_fields, axis=1)
    target_columns = [
        "resource_record_id",
        "protein_row_id",
        "sequence_object_id",
        "structure_unique_sequence_id",
        "target_variant_id",
        "target_concept_id",
        "uniprot_id",
        "gene_name_species",
        "gene_name",
        "protein_name",
        "species_name",
        "strict_target_key",
        "target_annotation_grade",
        "target_split_eligible",
    ]
    write_csv(
        target_annotations[target_columns],
        output_root / "data/entities/target_annotations.csv",
    )

    binder_clusters, binder_cluster_stats = build_binder_sequence_clusters(
        protein_objects
    )
    write_csv(
        binder_clusters,
        output_root / "data/diagnostics/binder_sequence_clusters_95.csv",
    )
    target_components, target_component_stats = build_strict_target_key_components(
        binder_clusters, target_annotations
    )
    write_csv(
        target_components,
        output_root / "data/diagnostics/strict_target_key_components_95.csv",
    )

    loop_objects = loops.merge(
        protein_objects[
            [
                "resource_record_id",
                "protein_row_id",
                "sequence_object_id",
                "structure_unique_sequence_id",
                "structure_source_layer",
            ]
        ],
        on=["protein_row_id", "structure_unique_sequence_id"],
        how="left",
        validate="many_to_one",
    )
    for column in ["loop_start", "loop_end", "loop_length", "domain_start"]:
        loop_objects[column] = pd.to_numeric(loop_objects[column], errors="raise").astype(int)
    sequences = protein_objects.set_index("protein_row_id")["final_sequence"].to_dict()
    loop_objects[
        [
            "coordinate_mapping_status",
            "canonical_loop_start",
            "canonical_loop_end",
            "canonical_loop_sequence",
        ]
    ] = loop_objects.apply(resolve_loop_coordinates, axis=1, sequences=sequences)
    loop_objects["split_participation"] = "inherits_parent_protein_split_only"
    loop_columns = [
        "loop_id",
        "resource_record_id",
        "protein_row_id",
        "sequence_object_id",
        "structure_unique_sequence_id",
        "family_name",
        "loop_slot",
        "loop_chain_id",
        "loop_start",
        "loop_end",
        "loop_length",
        "domain_start",
        "domain_end",
        "coordinate_mapping_status",
        "canonical_loop_start",
        "canonical_loop_end",
        "canonical_loop_sequence",
        "main_figure_use_class",
        "structure_source_layer",
        "split_participation",
    ]
    write_csv(
        loop_objects[loop_columns],
        output_root / "data/entities/recognition_window_objects.csv",
    )

    counts = {
        "protein_records": int(len(protein_objects)),
        "scaffold_families": int(protein_objects["scaffold_category"].nunique()),
        "unique_full_length_sequence_objects": int(
            protein_objects["sequence_object_id"].nunique()
        ),
        "structure_ready_records": int(
            protein_objects["submitted_structure_available"]
            .str.lower()
            .eq("true")
            .sum()
        ),
        "structure_instances": int(len(structure_instances)),
        "experimental_structure_instances": int(
            (structure_instances["structure_source_layer"] == "experimental_PDB").sum()
        ),
        "xray_em_core_candidates_pending_coordinate_qc": int(
            (
                structure_instances["ml_structure_grade"]
                == "xray_em_core_candidate_pending_coordinate_qc"
            ).sum()
        ),
        "nmr_flexibility_candidates": int(
            (
                structure_instances["ml_structure_grade"]
                == "nmr_flexibility_candidate_pending_model_qc"
            ).sum()
        ),
        "predicted_context_only": int(
            (structure_instances["ml_structure_grade"] == "predicted_context_only").sum()
        ),
        "target_annotated_records": int(target_annotations["target_variant_id"].ne("").sum()),
        "strict_target_split_eligible_records": int(
            target_annotations["target_split_eligible"].sum()
        ),
        "placeholder_target_records_excluded": int(
            (target_annotations["target_annotation_grade"] == "placeholder_excluded").sum()
        ),
        "recognition_window_objects": int(len(loop_objects)),
        "coordinate_resolved_window_objects": int(
            (loop_objects["coordinate_mapping_status"] != "unresolved").sum()
        ),
        "coordinate_unresolved_window_objects": int(
            (loop_objects["coordinate_mapping_status"] == "unresolved").sum()
        ),
        "binder_sequence_clusters_95": binder_cluster_stats["cluster_count"],
        "largest_binder_sequence_cluster_95": binder_cluster_stats[
            "largest_cluster_size"
        ],
        "provisional_strict_target_key_components_95": target_component_stats[
            "component_count"
        ],
        "largest_provisional_strict_target_key_component_95": target_component_stats[
            "largest_component_size"
        ],
    }
    scaffold_counts = {
        str(key): int(value)
        for key, value in protein_objects["scaffold_category"].value_counts().items()
    }
    summary = {
        "release_status": "foundation_audit_only_not_a_final_benchmark",
        "benchmark_scope": "full_length_sequences_across_12_scaffold_families",
        "window_policy": "recognition windows remain linked objects and inherit the parent protein split; they do not participate in clustering or partition assignment",
        "split_release_status": "no_fold_or_train_validation_test_assignment_released",
        "binder_cluster_diagnostics": binder_cluster_stats,
        "strict_target_key_component_diagnostics": target_component_stats,
        "counts": counts,
        "normalized_scaffold_label_counts": scaffold_counts,
        "blocking_conditions": [
            "Run atom-level coordinate completeness and residue-identity QC before freezing structure-core records.",
            "Retrieve and version target sequences, then cluster them before creating the combined binder-target split.",
            "Revalidate binder-chain and target-chain roles before releasing any complex-modeling benchmark.",
            "Do not treat predicted AF2 structures as experimental benchmark labels.",
        ],
    }
    metadata_dir = output_root / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    source_inventory = pd.DataFrame(
        [
            {
                "source_path": str(path.relative_to(source_root)).replace("\\", "/"),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in source_paths
        ]
    )
    write_csv(source_inventory, metadata_dir / "source_inventory.csv")
    (metadata_dir / "foundation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = [
        "# FBBP full-length ML foundation audit",
        "",
        "Status: **FOUNDATION ONLY - NOT A FINAL BENCHMARK**",
        "",
        "The benchmark scope is full-length sequence data across 12 scaffold families. Recognition-window objects are retained and inherit the parent protein split, but they are not used for clustering or partition assignment. The 95% binder clusters and strict-target-key components below are diagnostics only; no fold or train/validation/test assignment is released.",
        "",
        "## Counts",
        "",
    ]
    for key, value in counts.items():
        lines.append(f"- {key}: {value}")
    lines += ["", "## Normalized scaffold labels", ""]
    for key, value in scaffold_counts.items():
        lines.append(f"- {key}: {value}")
    lines += ["", "## Provisional graph diagnostics", ""]
    for key, value in binder_cluster_stats.items():
        lines.append(f"- binder_{key}: {value}")
    for key, value in target_component_stats.items():
        lines.append(f"- strict_target_key_{key}: {value}")
    lines += ["", "## Blocking conditions", ""]
    for item in summary["blocking_conditions"]:
        lines.append(f"- {item}")
    lines.append("")
    reports_dir = output_root / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "foundation_audit.md").write_text("\n".join(lines), encoding="utf-8")

    checksum_lines = []
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "checksums_sha256.txt":
            relative = str(path.relative_to(output_root)).replace("\\", "/")
            checksum_lines.append(f"{sha256_file(path)}  {relative}")
    (metadata_dir / "checksums_sha256.txt").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.source_root, args.output_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
