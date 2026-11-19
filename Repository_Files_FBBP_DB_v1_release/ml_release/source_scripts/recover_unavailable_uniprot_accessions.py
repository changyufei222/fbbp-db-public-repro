#!/usr/bin/env python
"""Recover obsolete UniProtKB accessions through the official UniParc ID mapping API."""

from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FAILURES = ROOT / "FBBP_ML_target_sequence_aware_10fold_v1/data/target_accessions_without_frozen_sequence.csv"
DEFAULT_OUTPUT = ROOT / "FBBP_ML_uniprot_uniparc_recovery_20260716"
USER_AGENT = "FBBP-DB-ML/20260716"


def request_json(url: str, data: bytes | None = None) -> dict:
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    request = Request(url, data=data, headers=headers)
    with urlopen(request, timeout=90) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_text(url: str) -> str:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=90) as response:
        return response.read().decode("utf-8")


def submit_mapping(accessions: list[str]) -> str:
    payload = urlencode(
        {
            "from": "UniProtKB_AC-ID",
            "to": "UniParc",
            "ids": ",".join(accessions),
        }
    ).encode("utf-8")
    result = request_json("https://rest.uniprot.org/idmapping/run", data=payload)
    return result["jobId"]


def wait_for_job(job_id: str) -> None:
    for _ in range(90):
        result = request_json(f"https://rest.uniprot.org/idmapping/status/{job_id}")
        if result.get("jobStatus") == "FINISHED" or "results" in result or "failedIds" in result:
            return
        if result.get("jobStatus") not in {"RUNNING", "NEW"}:
            raise RuntimeError(f"UniProt mapping failed: {result}")
        time.sleep(1)
    raise TimeoutError(f"UniProt mapping job did not finish: {job_id}")


def extract_uniparc_id(value: object) -> str:
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        return ""
    for key in ("uniParcId", "id", "primaryAccession"):
        if value.get(key):
            return str(value[key])
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fasta_dir = args.output / "fasta"
    fasta_dir.mkdir(exist_ok=True)
    failures = pd.read_csv(FAILURES, dtype=str).fillna("")
    accessions = sorted(failures["uniprot_accession"].unique())
    job_id = submit_mapping(accessions)
    wait_for_job(job_id)
    result_url = f"https://rest.uniprot.org/idmapping/uniparc/results/stream/{job_id}?format=json"
    mapping_result = request_json(result_url)
    (args.output / "idmapping_raw_response.json").write_text(
        json.dumps(mapping_result, indent=2), encoding="utf-8"
    )
    mappings: dict[str, set[str]] = {accession: set() for accession in accessions}
    for item in mapping_result.get("results", []):
        source = str(item.get("from", "")).upper()
        uniparc_id = extract_uniparc_id(item.get("to"))
        if source in mappings and uniparc_id:
            mappings[source].add(uniparc_id)
    rows = []
    for accession in accessions:
        candidates = sorted(mappings[accession])
        status = "unresolved"
        sequence = ""
        selected_uniparc_id = ""
        detail = ""
        if len(candidates) == 1:
            selected_uniparc_id = candidates[0]
            fasta = fetch_text(f"https://rest.uniprot.org/uniparc/{selected_uniparc_id}.fasta")
            lines = [line.strip() for line in fasta.splitlines() if line.strip()]
            sequence = "".join(lines[1:]).upper() if lines and lines[0].startswith(">") else ""
            if sequence and set(sequence).issubset(set("ACDEFGHIKLMNPQRSTVWY")):
                status = "recovered_unique_uniparc_sequence"
                (fasta_dir / f"{accession}__{selected_uniparc_id}.fasta").write_text(
                    fasta.strip() + "\n", encoding="utf-8"
                )
            else:
                detail = "UniParc FASTA missing or contains noncanonical residues"
        elif len(candidates) > 1:
            status = "ambiguous_multiple_uniparc_sequences"
            detail = ";".join(candidates)
        else:
            detail = "ID Mapping returned no UniParc record"
        rows.append(
            {
                "uniprot_accession": accession,
                "recovery_status": status,
                "selected_uniparc_id": selected_uniparc_id,
                "candidate_uniparc_ids": ";".join(candidates),
                "sequence_length": len(sequence),
                "sequence": sequence,
                "failure_detail": detail,
                "mapping_job_id": job_id,
                "recovery_date": str(date.today()),
            }
        )
    pd.DataFrame(rows).to_csv(args.output / "unavailable_accession_uniparc_recovery.csv", index=False)
    print(args.output)


if __name__ == "__main__":
    main()
