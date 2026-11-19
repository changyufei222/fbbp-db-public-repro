#!/usr/bin/env python
"""Recompute the PCSK9/1459D05 contacts and buried-SASA source tables.

The input is PDB 4OV6 clipped to chains D, E, and G. PCSK9 chains D/E are the
target and Adnectin chain G is the binder. Contacts are non-hydrogen atom-pair
contacts at the specified cutoff and are distinct from unique residue-pair
counts used for database-wide FRR contact summaries.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley


FRRS = {
    "L2": (22, 30),
    "L6": (60, 66),
    "L7": (76, 83),
}


def standard_residues(chain):
    return [residue for residue in chain if residue.id[0] == " "]


def heavy_atoms(residue):
    for atom in residue:
        element = (atom.element or "").strip().upper()
        if element in {"H", "D"}:
            continue
        if not element and atom.get_name().upper().startswith(("H", "D")):
            continue
        yield atom


def frr_assignment(residue_number):
    for frr, (start, end) in FRRS.items():
        if start <= residue_number <= end:
            return frr
    return "non_FRR_interface"


def pruned_structure(structure, chain_ids):
    copied = deepcopy(structure)
    model = copied[0]
    for chain_id in list(model.child_dict):
        if chain_id not in chain_ids:
            model.detach_child(chain_id)
    return copied


def residue_key(residue):
    return residue.id[1], residue.id[2].strip(), residue.resname


def mean_asa(keys, values):
    usable = [key for key in keys if key in values]
    return sum(values[key] for key in usable) / len(usable) if usable else 0.0


def compute_contacts(model, binder_chain, target_chains, cutoff):
    binder_atoms = [
        atom
        for residue in standard_residues(model[binder_chain])
        for atom in heavy_atoms(residue)
    ]
    target_atoms = [
        atom
        for chain_id in target_chains
        for residue in standard_residues(model[chain_id])
        for atom in heavy_atoms(residue)
    ]

    search = NeighborSearch(target_atoms)
    rows = []
    for binder_atom in binder_atoms:
        binder_residue = binder_atom.get_parent()
        for target_atom in search.search(binder_atom.coord, cutoff, "A"):
            distance = float(binder_atom - target_atom)
            if distance > cutoff:
                continue
            target_residue = target_atom.get_parent()
            rows.append(
                {
                    "pdb_id": "4OV6",
                    "binder_chain": binder_chain,
                    "target_chain": target_residue.get_parent().id,
                    "cutoff_A": f"{cutoff:.1f}",
                    "binder_residue_number": binder_residue.id[1],
                    "binder_insertion_code": binder_residue.id[2].strip(),
                    "binder_residue_name": binder_residue.resname,
                    "binder_atom_name": binder_atom.get_name(),
                    "target_residue_number": target_residue.id[1],
                    "target_insertion_code": target_residue.id[2].strip(),
                    "target_residue_name": target_residue.resname,
                    "target_atom_name": target_atom.get_name(),
                    "distance_A": round(distance, 3),
                    "frr": frr_assignment(binder_residue.id[1]),
                }
            )

    rows.sort(
        key=lambda row: (
            row["binder_residue_number"],
            row["binder_insertion_code"],
            row["binder_atom_name"],
            row["target_chain"],
            row["target_residue_number"],
            row["target_insertion_code"],
            row["target_atom_name"],
        )
    )
    return rows


def compute_delta_asa(structure, binder_chain, target_chains):
    binder_alone = pruned_structure(structure, {binder_chain})
    complex_pair = pruned_structure(structure, {binder_chain, *target_chains})
    calculator = ShrakeRupley(n_points=960, probe_radius=1.4)
    calculator.compute(binder_alone, level="R")
    calculator.compute(complex_pair, level="R")

    alone = {
        residue_key(residue): float(getattr(residue, "sasa", 0.0))
        for residue in standard_residues(binder_alone[0][binder_chain])
    }
    in_complex = {
        residue_key(residue): float(getattr(residue, "sasa", 0.0))
        for residue in standard_residues(complex_pair[0][binder_chain])
    }
    return alone, in_complex


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--binder-chain", default="G")
    parser.add_argument("--target-chains", default="D,E")
    parser.add_argument("--cutoff", default=4.5, type=float)
    args = parser.parse_args()

    target_chains = [value.strip() for value in args.target_chains.split(",") if value.strip()]
    structure = PDBParser(QUIET=True).get_structure("4OV6", str(args.pdb))
    model = structure[0]
    contacts = compute_contacts(model, args.binder_chain, target_chains, args.cutoff)
    alone, in_complex = compute_delta_asa(structure, args.binder_chain, target_chains)

    atom_pairs_by_residue = defaultdict(int)
    target_residues_by_residue = defaultdict(set)
    for row in contacts:
        key = (
            row["binder_residue_number"],
            row["binder_insertion_code"],
            row["binder_residue_name"],
        )
        atom_pairs_by_residue[key] += 1
        target_residues_by_residue[key].add(
            (
                row["target_chain"],
                row["target_residue_number"],
                row["target_insertion_code"],
                row["target_residue_name"],
            )
        )

    residue_rows = []
    for residue in standard_residues(model[args.binder_chain]):
        key = residue_key(residue)
        delta_asa = max(alone[key] - in_complex[key], 0.0)
        atom_pair_count = atom_pairs_by_residue[key]
        if atom_pair_count == 0 and delta_asa <= 0.001:
            continue
        residue_rows.append(
            {
                "pdb_id": "4OV6",
                "binder_chain": args.binder_chain,
                "target_chains": ";".join(target_chains),
                "frr": frr_assignment(key[0]),
                "residue_number": key[0],
                "insertion_code": key[1],
                "residue_name": key[2],
                "heavy_atom_atom_pair_count_4p5A": atom_pair_count,
                "unique_target_residue_count_4p5A": len(target_residues_by_residue[key]),
                "asa_binder_alone_A2_biopython_sr960": round(alone[key], 3),
                "asa_binder_in_complex_A2_biopython_sr960": round(in_complex[key], 3),
                "deltaASA_A2_biopython_sr960": round(delta_asa, 3),
            }
        )

    summary_rows = []
    for frr, (start, end) in FRRS.items():
        frr_contacts = [row for row in contacts if row["frr"] == frr]
        frr_residues = [row for row in residue_rows if row["frr"] == frr]
        loop_keys = [
            residue_key(residue)
            for residue in standard_residues(model[args.binder_chain])
            if start <= residue.id[1] <= end
        ]
        mean_alone = mean_asa(loop_keys, alone)
        mean_complex = mean_asa(loop_keys, in_complex)
        summary_rows.append(
            {
                "entity": frr,
                "residue_range": f"{start}-{end}",
                "binder_contact_residue_count_4p5A": len(
                    {
                        (row["binder_residue_number"], row["binder_insertion_code"])
                        for row in frr_contacts
                    }
                ),
                "heavy_atom_atom_pair_count_4p5A": len(frr_contacts),
                "target_residue_count_4p5A": len(
                    {
                        (row["target_chain"], row["target_residue_number"], row["target_insertion_code"])
                        for row in frr_contacts
                    }
                ),
                "deltaASA_A2_sum_biopython_sr960": round(
                    sum(row["deltaASA_A2_biopython_sr960"] for row in frr_residues), 3
                ),
                "deltaASA_A2_rounded_1dp": round(
                    sum(row["deltaASA_A2_biopython_sr960"] for row in frr_residues), 1
                ),
                "mean_asa_binder_alone_A2_biopython_sr960": round(mean_alone, 3),
                "mean_asa_binder_in_complex_A2_biopython_sr960": round(mean_complex, 3),
                "mean_deltaASA_A2_biopython_sr960": round(mean_alone - mean_complex, 3),
            }
        )

    frr_contacts = [row for row in contacts if row["frr"] in FRRS]
    frr_residues = [row for row in residue_rows if row["frr"] in FRRS]
    all_binder_residues = [row for row in residue_rows]
    all_frr_keys = [
        residue_key(residue)
        for residue in standard_residues(model[args.binder_chain])
        if any(start <= residue.id[1] <= end for start, end in FRRS.values())
    ]
    all_interface_keys = [
        (row["residue_number"], row["insertion_code"], row["residue_name"])
        for row in all_binder_residues
    ]
    all_frr_mean_alone = mean_asa(all_frr_keys, alone)
    all_frr_mean_complex = mean_asa(all_frr_keys, in_complex)
    all_interface_mean_alone = mean_asa(all_interface_keys, alone)
    all_interface_mean_complex = mean_asa(all_interface_keys, in_complex)
    summary_rows.append(
        {
            "entity": "all_FRRs",
            "residue_range": "L2+L6+L7",
            "binder_contact_residue_count_4p5A": len(
                {
                    (row["binder_residue_number"], row["binder_insertion_code"])
                    for row in frr_contacts
                }
            ),
            "heavy_atom_atom_pair_count_4p5A": len(frr_contacts),
            "target_residue_count_4p5A": len(
                {
                    (row["target_chain"], row["target_residue_number"], row["target_insertion_code"])
                    for row in frr_contacts
                }
            ),
            "deltaASA_A2_sum_biopython_sr960": round(
                sum(row["deltaASA_A2_biopython_sr960"] for row in frr_residues), 3
            ),
            "deltaASA_A2_rounded_1dp": round(
                sum(row["deltaASA_A2_biopython_sr960"] for row in frr_residues), 1
            ),
            "mean_asa_binder_alone_A2_biopython_sr960": round(all_frr_mean_alone, 3),
            "mean_asa_binder_in_complex_A2_biopython_sr960": round(all_frr_mean_complex, 3),
            "mean_deltaASA_A2_biopython_sr960": round(
                all_frr_mean_alone - all_frr_mean_complex, 3
            ),
        }
    )
    summary_rows.append(
        {
            "entity": "all_binder_interface",
            "residue_range": "all chain G interface residues",
            "binder_contact_residue_count_4p5A": len(
                {
                    (row["binder_residue_number"], row["binder_insertion_code"])
                    for row in contacts
                }
            ),
            "heavy_atom_atom_pair_count_4p5A": len(contacts),
            "target_residue_count_4p5A": len(
                {
                    (row["target_chain"], row["target_residue_number"], row["target_insertion_code"])
                    for row in contacts
                }
            ),
            "deltaASA_A2_sum_biopython_sr960": round(
                sum(row["deltaASA_A2_biopython_sr960"] for row in all_binder_residues), 3
            ),
            "deltaASA_A2_rounded_1dp": round(
                sum(row["deltaASA_A2_biopython_sr960"] for row in all_binder_residues), 1
            ),
            "mean_asa_binder_alone_A2_biopython_sr960": round(all_interface_mean_alone, 3),
            "mean_asa_binder_in_complex_A2_biopython_sr960": round(
                all_interface_mean_complex, 3
            ),
            "mean_deltaASA_A2_biopython_sr960": round(
                all_interface_mean_alone - all_interface_mean_complex, 3
            ),
        }
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "pcsk9_1459D05_heavy_atom_contacts.tsv", contacts)
    write_tsv(args.outdir / "pcsk9_1459D05_residue_burial_deltaASA.tsv", residue_rows)
    write_tsv(args.outdir / "pcsk9_1459D05_loop_interface_summary.tsv", summary_rows)


if __name__ == "__main__":
    main()
