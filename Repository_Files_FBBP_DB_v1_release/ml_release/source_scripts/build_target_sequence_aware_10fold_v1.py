#!/usr/bin/env python
"""Freeze UniProt target sequences and build the primary 40% target-aware 10-fold split."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import time
from collections import defaultdict
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd
from Bio.Align import PairwiseAligner


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = ROOT / "FBBP_ML_foundation_20260714"
BINDER_RELEASE = ROOT / "FBBP_ML_binder_full_length_10fold_v1"
STAGING = ROOT / "FBBP_ML_target_sequence_freeze_staging_20260714"
NAME_RESOLUTION = ROOT / "FBBP_ML_named_target_uniprot_resolution_v2_20260716/named_target_uniprot_resolution.csv"
UNIPARC_RECOVERY = ROOT / "FBBP_ML_uniprot_uniparc_recovery_20260716/unavailable_accession_uniparc_recovery.csv"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_target_sequence_aware_10fold_v1"
RELEASE_NAME = "target_sequence_aware_10fold_v1"
SEED = "20260714"
FOLD_COUNT = 10
THRESHOLDS = {"30": 0.30, "40": 0.40, "50": 0.50}


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


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_accessions(value: str) -> list[str]:
    return [item.strip().upper() for item in str(value).split(";") if item.strip()]


def fetch_fasta(accession: str) -> tuple[str, str, str]:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    request = Request(url, headers={"User-Agent": "FBBP-DB-ML/20260714"})
    with urlopen(request, timeout=60) as response:
        raw = response.read()
    text = raw.decode("utf-8").strip() + "\n"
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Invalid FASTA response for {accession}")
    sequence = "".join(lines[1:]).strip().upper()
    if not sequence or re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
        raise ValueError(f"Invalid protein sequence for {accession}")
    return text, lines[0][1:], sequence


def freeze_sequences(manifest: pd.DataFrame, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    fasta_dir = output / "data/fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    registry = []
    failures = []
    recovery = {}
    if UNIPARC_RECOVERY.exists():
        recovery_frame = pd.read_csv(UNIPARC_RECOVERY, dtype=str).fillna("")
        recovery = recovery_frame.set_index("uniprot_accession").to_dict("index")
    for index, row in manifest.sort_values("uniprot_accession").reset_index(drop=True).iterrows():
        accession = row["uniprot_accession"]
        last_error = None
        sequence = None
        retrieval_source = "UniProtKB REST /uniprotkb/{accession}.fasta"
        source_record_version = ""
        for attempt in range(3):
            try:
                raw, header, sequence = fetch_fasta(accession)
                break
            except Exception as error:  # Network and server failures are retried deterministically.
                last_error = error
                if attempt == 2:
                    raw = header = sequence = None
                    break
                time.sleep(2 * (attempt + 1))
        if sequence is None:
            archived = recovery.get(accession, {})
            if archived.get("recovery_status") == "recovered_unique_uniparc_sequence":
                sequence = archived["sequence"]
                uniparc_id = archived["selected_uniparc_id"]
                header = f"archived|{accession}|{uniparc_id} recovered by UniProt ID Mapping"
                raw = f">{header}\n{sequence}\n"
                retrieval_source = "UniProt ID Mapping UniProtKB_AC-ID -> UniParc; UniParc FASTA"
                source_record_version = uniparc_id
            else:
                failures.append({
                    "uniprot_accession": accession,
                    "retrieval_status": "sequence_unavailable_or_retrieval_failed",
                    "failure_detail": str(last_error),
                    "retrieval_date": str(date.today()),
                })
                continue
        path = fasta_dir / f"{accession}.fasta"
        path.write_text(raw, encoding="utf-8")
        version_match = re.search(r"\bSV=(\d+)\b", header)
        if version_match:
            source_record_version = version_match.group(1)
        target_sequence_id = "TSEQ_" + sha256_bytes(sequence.encode("ascii"))[:12].upper()
        registry.append(
            {
                "target_sequence_id": target_sequence_id,
                "uniprot_accession": accession,
                "fasta_header": header,
                "sequence_length": len(sequence),
                "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
                "source_record_version": source_record_version,
                "retrieval_date": str(date.today()),
                "retrieval_source": retrieval_source,
                "fasta_file": f"data/fasta/{accession}.fasta",
                "sequence": sequence,
            }
        )
        if index < len(manifest) - 1:
            time.sleep(0.1)
    result = pd.DataFrame(registry).sort_values("uniprot_accession").reset_index(drop=True)
    result.to_csv(output / "data/frozen_target_sequence_registry.csv", index=False)
    with (output / "data/target_sequences.fasta").open("w", encoding="utf-8") as handle:
        for row in result.itertuples(index=False):
            handle.write(f">{row.target_sequence_id}|{row.uniprot_accession}|SV={row.source_record_version}\n")
            handle.write(f"{row.sequence}\n")
    failures_frame = pd.DataFrame(failures)
    failures_frame.to_csv(output / "data/target_accessions_without_frozen_sequence.csv", index=False)
    return result, failures_frame


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
    alignment = aligner.align(left, right)[0]
    left_blocks, right_blocks = alignment.aligned
    matches = 0
    for (left_start, left_end), (right_start, right_end) in zip(left_blocks, right_blocks):
        matches += sum(
            a == b
            for a, b in zip(left[int(left_start):int(left_end)], right[int(right_start):int(right_end)])
        )
    return matches / longer


def target_clusters(registry: pd.DataFrame, threshold: float, label: str) -> tuple[pd.DataFrame, dict]:
    rows = registry.sort_values("uniprot_accession").to_dict("records")
    union_find = UnionFind([row["uniprot_accession"] for row in rows])
    aligner = make_aligner()
    compared = 0
    linked = 0
    for index, left in enumerate(rows):
        for right in rows[index + 1 :]:
            if min(len(left["sequence"]), len(right["sequence"])) / max(len(left["sequence"]), len(right["sequence"])) < threshold:
                continue
            compared += 1
            if sequence_identity(left["sequence"], right["sequence"], aligner, threshold) >= threshold:
                union_find.union(left["uniprot_accession"], right["uniprot_accession"])
                linked += 1
    groups = sorted(union_find.groups().values(), key=lambda members: (min(members), len(members)))
    cluster_for_accession = {}
    for index, members in enumerate(groups, start=1):
        cluster_id = f"TGSEQ{label}_{index:05d}"
        for accession in members:
            cluster_for_accession[accession] = cluster_id
    result = registry[["target_sequence_id", "uniprot_accession", "sequence_length", "sequence_sha256"]].copy()
    result[f"target_cluster_{label}_id"] = result["uniprot_accession"].map(cluster_for_accession)
    result[f"target_cluster_{label}_size"] = result.groupby(f"target_cluster_{label}_id")["uniprot_accession"].transform("size")
    result["identity_threshold"] = threshold
    result["clustering_method"] = "global_affine_alignment_exact_matches_divided_by_longer_sequence_length"
    return result, {"threshold": threshold, "cluster_count": len(groups), "comparison_count": compared, "linked_pair_count": linked}


def exact_accession_clusters(registry: pd.DataFrame) -> pd.DataFrame:
    result = registry[["target_sequence_id", "uniprot_accession", "sequence_length", "sequence_sha256"]].copy()
    result["target_cluster_exact_accession_id"] = "TGEXACT_" + result["uniprot_accession"]
    result["target_cluster_exact_accession_size"] = 1
    result["identity_threshold"] = "exact_accession"
    result["clustering_method"] = "exact_query_accession"
    return result


def build_links(registry: pd.DataFrame, failures: pd.DataFrame, name_resolution: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    targets = pd.read_csv(FOUNDATION / "data/entities/target_annotations.csv", dtype=str).fillna("")
    clusters = pd.read_csv(
        BINDER_RELEASE / "data/binder_sequence_clusters_95_affine_v1.csv",
        dtype=str,
    ).fillna("")
    proteins = pd.read_csv(FOUNDATION / "data/entities/protein_objects.csv", dtype=str).fillna("")
    strict = targets[targets["target_split_eligible"].str.lower() == "true"].copy()
    accession_rows = strict[strict["uniprot_id"].str.strip() != ""].copy()
    accession_rows["uniprot_accession"] = accession_rows["uniprot_id"].map(split_accessions)
    accession_rows = accession_rows.explode("uniprot_accession")
    resolved_statuses = {
        "auto_mapped_unique",
        "curated_mapped_unique",
        "curated_mapped_multi_component",
    }
    mapped_variants = name_resolution[name_resolution["resolution_status"].isin(resolved_statuses)][
        ["target_variant_id", "selected_uniprot_accession"]
    ].rename(columns={"selected_uniprot_accession": "uniprot_accession"})
    mapped_variants["uniprot_accession"] = mapped_variants["uniprot_accession"].map(split_accessions)
    mapped_variants = mapped_variants.explode("uniprot_accession")
    named_rows = strict[strict["uniprot_id"].str.strip() == ""].merge(
        mapped_variants, on="target_variant_id", how="inner", validate="many_to_many"
    )
    candidate_rows = pd.concat([accession_rows, named_rows], ignore_index=True)
    links = candidate_rows.merge(
        registry[["target_sequence_id", "uniprot_accession"]], on="uniprot_accession", how="inner", validate="many_to_one"
    ).merge(
        clusters[["resource_record_id", "binder_cluster_95_id", "binder_cluster_size"]], on="resource_record_id", how="left", validate="many_to_one"
    ).merge(
        proteins[["resource_record_id", "scaffold_category"]], on="resource_record_id", how="left", validate="many_to_one"
    )
    named_excluded = strict[strict["uniprot_id"].str.strip() == ""].merge(
        name_resolution[["target_variant_id", "resolution_status"]], on="target_variant_id", how="left", validate="many_to_one"
    )
    named_excluded = named_excluded[~named_excluded["resolution_status"].isin(resolved_statuses)].copy()
    failed_records = candidate_rows.merge(
        failures[["uniprot_accession", "retrieval_status", "failure_detail"]],
        on="uniprot_accession", how="inner"
    )
    named_excluded["exclusion_reason"] = "name_target_unresolved_or_ambiguous"
    failed_records["exclusion_reason"] = "accession_resolved_but_sequence_unavailable"
    excluded_columns = ["resource_record_id", "target_variant_id", "target_concept_id", "uniprot_id", "gene_name_species", "gene_name", "protein_name", "species_name", "strict_target_key", "exclusion_reason"]
    for column in excluded_columns:
        if column not in named_excluded:
            named_excluded[column] = ""
        if column not in failed_records:
            failed_records[column] = ""
    excluded = pd.concat([named_excluded[excluded_columns], failed_records[excluded_columns]], ignore_index=True)
    covered = set(links["resource_record_id"]) | set(failed_records["resource_record_id"]) | set(named_excluded["resource_record_id"])
    if covered != set(strict["resource_record_id"]):
        raise ValueError("Some strict target records were neither frozen nor explicitly excluded")
    return links, excluded


def components(links: pd.DataFrame, target_cluster_column: str) -> pd.DataFrame:
    node_ids = sorted(links["resource_record_id"].unique())
    union_find = UnionFind(node_ids)
    for _, group in links.groupby("binder_cluster_95_id"):
        members = sorted(group["resource_record_id"].unique())
        for member in members[1:]:
            union_find.union(members[0], member)
    for _, group in links.groupby(target_cluster_column):
        members = sorted(group["resource_record_id"].unique())
        for member in members[1:]:
            union_find.union(members[0], member)
    groups = sorted(union_find.groups().values(), key=lambda members: (min(members), len(members)))
    component_for_record = {}
    component_size = {}
    for index, members in enumerate(groups, start=1):
        component_id = f"TACOMP_{index:05d}"
        for member in members:
            component_for_record[member] = component_id
            component_size[member] = len(members)
    result = links.copy()
    result["target_aware_component_id"] = result["resource_record_id"].map(component_for_record)
    result["component_size"] = result["resource_record_id"].map(component_size)
    return result


def stable_rank(component_id: str, fold_id: int) -> int:
    return int(hashlib.sha256(f"{SEED}|{component_id}|{fold_id}".encode()).hexdigest()[:16], 16)


def assign_components(primary: pd.DataFrame) -> pd.DataFrame:
    components_table = (
        primary.groupby("target_aware_component_id", as_index=False)
        .agg(component_size=("resource_record_id", "nunique"))
    )
    component_families = primary.pivot_table(
        index="target_aware_component_id", columns="scaffold_category", values="resource_record_id", aggfunc="nunique", fill_value=0
    )
    components_table = components_table.join(component_families, on="target_aware_component_id")
    families = sorted(primary["scaffold_category"].unique())
    components_table = components_table.sort_values(["component_size", "target_aware_component_id"], ascending=[False, True]).reset_index(drop=True)
    family_targets = primary.groupby("scaffold_category")["resource_record_id"].nunique().to_dict()
    family_targets = {family: value / FOLD_COUNT for family, value in family_targets.items()}
    fold_total = [0] * FOLD_COUNT
    fold_family = {family: [0] * FOLD_COUNT for family in families}
    mapping = {}
    for index, row in components_table.iterrows():
        candidates = []
        for fold_id in range(FOLD_COUNT):
            total_score = ((fold_total[fold_id] + row.component_size - len(primary["resource_record_id"].unique()) / FOLD_COUNT) / max(1, len(primary["resource_record_id"].unique()) / FOLD_COUNT)) ** 2
            family_score = sum(
                ((fold_family[family][fold_id] + row.get(family, 0) - family_targets[family]) / max(1, family_targets[family])) ** 2
                for family in families
            )
            empty_penalty = 0 if index >= FOLD_COUNT or fold_total[fold_id] == 0 else 1
            candidates.append((empty_penalty, fold_total[fold_id], family_score, stable_rank(row.target_aware_component_id, fold_id), fold_id))
        fold_id = min(candidates)[-1]
        mapping[row.target_aware_component_id] = fold_id
        fold_total[fold_id] += int(row.component_size)
        for family in families:
            fold_family[family][fold_id] += int(row.get(family, 0))
    result = primary.copy()
    result["fold_id"] = result["target_aware_component_id"].map(mapping).astype(int)
    return result


def write_release(output: Path, registry: pd.DataFrame, links: pd.DataFrame, excluded: pd.DataFrame, diagnostics: dict) -> None:
    primary = assign_components(components(links, "target_cluster_40_id"))
    assignment = primary[[
        "resource_record_id", "scaffold_category", "binder_cluster_95_id", "binder_cluster_size", "uniprot_accession", "target_sequence_id", "target_cluster_40_id", "target_aware_component_id", "component_size", "fold_id"
    ]].drop_duplicates().sort_values(["resource_record_id", "uniprot_accession"])
    assignment["release_name"] = RELEASE_NAME
    assignment.to_csv(output / "data/target_sequence_aware_10fold_v1.csv", index=False)
    links.to_csv(output / "data/protein_target_sequence_links.csv", index=False)
    excluded.to_csv(output / "data/named_targets_excluded_from_target_sequence_aware_split.csv", index=False)
    for label, cluster_frame in diagnostics["cluster_frames"].items():
        cluster_frame.to_csv(output / f"data/target_sequence_clusters_{label}.csv", index=False)
        component_column = f"target_cluster_{label}_id" if label != "exact_accession" else "target_cluster_exact_accession_id"
        if label == "40":
            component_input = links
        else:
            component_input = links.merge(
                cluster_frame[["uniprot_accession", component_column]],
                on="uniprot_accession",
                how="left",
            )
        component_frame = components(component_input, component_column)
        component_frame[["resource_record_id", "uniprot_accession", component_column, "target_aware_component_id", "component_size"]].drop_duplicates().to_csv(output / f"data/target_aware_components_{label}.csv", index=False)
    fold_sizes = assignment.groupby("fold_id")["resource_record_id"].nunique().reindex(range(FOLD_COUNT), fill_value=0)
    component_crossings = assignment.groupby("target_aware_component_id")["fold_id"].nunique()
    binder_crossings = assignment.groupby("binder_cluster_95_id")["fold_id"].nunique()
    target_crossings = assignment.groupby("target_cluster_40_id")["fold_id"].nunique()
    family_balance = assignment.groupby(["scaffold_category", "fold_id"])["resource_record_id"].nunique().reset_index(name="protein_objects")
    family_balance.to_csv(output / "reports/target_sequence_aware_10fold_family_balance.csv", index=False)
    audit = {
        "release_name": RELEASE_NAME,
        "status": "PASS",
        "release_date": str(date.today()),
        "primary_target_identity_threshold": 0.40,
        "sensitivity_versions": [0.30, 0.50, "exact_accession"],
        "frozen_target_sequences": int(len(registry)),
        "target_sequence_resolved_protein_objects": int(assignment["resource_record_id"].nunique()),
        "target_aware_component_count": int(assignment["target_aware_component_id"].nunique()),
        "largest_target_aware_component": int(assignment["component_size"].max()),
        "records_excluded_from_target_sequence_graph": int(excluded["resource_record_id"].nunique()),
        "name_target_records_excluded": int(excluded.loc[
            excluded["exclusion_reason"] == "name_target_unresolved_or_ambiguous",
            "resource_record_id",
        ].nunique()),
        "accession_records_excluded_sequence_unavailable": int(excluded.loc[
            excluded["exclusion_reason"] == "accession_resolved_but_sequence_unavailable",
            "resource_record_id",
        ].nunique()),
        "uniparc_recovered_target_sequences": int(
            registry["retrieval_source"].str.contains("UniParc", na=False).sum()
        ),
        "fold_sizes": {str(index): int(value) for index, value in fold_sizes.items()},
        "binder_cluster_crossings": int((binder_crossings > 1).sum()),
        "target_cluster_40_crossings": int((target_crossings > 1).sum()),
        "component_crossings": int((component_crossings > 1).sum()),
        "cluster_diagnostics": diagnostics["stats"],
        "not_created": ["benchmark_zip", "complex_grade_structure_benchmark"],
        "derived_train_validation_test_view": "FBBP_ML_derived_split_views_20260716",
    }
    if audit["binder_cluster_crossings"] or audit["target_cluster_40_crossings"] or audit["component_crossings"]:
        raise RuntimeError("Target-aware graph crosses a fold boundary")
    if sorted(fold_sizes.index[fold_sizes > 0].tolist()) != list(range(FOLD_COUNT)):
        raise RuntimeError("Not all target-aware folds are populated")
    (output / "reports/target_sequence_aware_10fold_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# FBBP target-sequence-aware 10-fold audit", "", "Status: **PASS - frozen target sequences and 40% target-aware fold IDs**", ""]
    report += [f"- {key}: {value}" for key, value in audit.items() if key not in {"cluster_diagnostics", "not_created"}]
    report += ["- Unresolved or ambiguous name-only target records are excluded rather than inferred.", "- No train/validation/test CSV or benchmark ZIP is created."]
    (output / "reports/target_sequence_aware_10fold_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# FBBP target-sequence-aware 10-fold release\n\n"
        "Primary target clustering uses 40% sequence identity after UniProt sequence freezing. "
        "Sensitivity diagnostics use 30%, 50%, and exact accession. Ambiguous name-only targets are excluded.\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.resume:
        shutil.rmtree(args.output)
    for directory in ("data", "reports", "metadata"):
        (args.output / directory).mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(STAGING / "data/target_accession_freeze_manifest.csv", dtype=str).fillna("")
    name_resolution = pd.read_csv(NAME_RESOLUTION, dtype=str).fillna("")
    resolved_statuses = {
        "auto_mapped_unique",
        "curated_mapped_unique",
        "curated_mapped_multi_component",
    }
    name_accessions = name_resolution[
        name_resolution["resolution_status"].isin(resolved_statuses)
    ]["selected_uniprot_accession"].map(split_accessions).explode()
    manifest = pd.concat([manifest, pd.DataFrame({"uniprot_accession": name_accessions})], ignore_index=True)
    manifest = manifest.drop_duplicates("uniprot_accession").reset_index(drop=True)
    registry_path = args.output / "data/frozen_target_sequence_registry.csv"
    failure_path = args.output / "data/target_accessions_without_frozen_sequence.csv"
    if args.resume and registry_path.exists() and failure_path.exists():
        registry = pd.read_csv(registry_path, dtype=str).fillna("")
        retrieval_failures = pd.read_csv(failure_path, dtype=str).fillna("")
    else:
        # UniProt confirms the accession-qualified isoform P00533-1, so include it verbatim.
        registry, retrieval_failures = freeze_sequences(manifest, args.output)
    links, excluded = build_links(registry, retrieval_failures, name_resolution)
    diagnostics = {"cluster_frames": {}, "stats": {}}
    for label, threshold in THRESHOLDS.items():
        frame, stats = target_clusters(registry, threshold, label)
        diagnostics["cluster_frames"][label] = frame
        diagnostics["stats"][label] = stats
    exact = exact_accession_clusters(registry)
    diagnostics["cluster_frames"]["exact_accession"] = exact
    diagnostics["stats"]["exact_accession"] = {"cluster_count": int(len(exact)), "threshold": "exact_accession"}
    links = links.merge(diagnostics["cluster_frames"]["40"][["uniprot_accession", "target_cluster_40_id"]], on="uniprot_accession", how="left")
    write_release(args.output, registry, links, excluded, diagnostics)
    checksum_path = args.output / "metadata/checksums_sha256.json"
    manifest_files = sorted(
        path for path in args.output.rglob("*") if path.is_file() and path != checksum_path
    )
    checksum_path.write_text(
        json.dumps([{"path": path.relative_to(args.output).as_posix(), "sha256": sha256_file(path)} for path in manifest_files], indent=2),
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
