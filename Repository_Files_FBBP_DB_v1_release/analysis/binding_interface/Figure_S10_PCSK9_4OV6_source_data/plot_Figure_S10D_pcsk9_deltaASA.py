#!/usr/bin/env python
"""Plot the PCSK9/4OV6 Figure S10D loop-level ΔASA bars from source data."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


ORDER = ["L2", "L6", "L7"]
LABELS = {"L2": "L2 (Loop 1)", "L6": "L6 (Loop 2)", "L7": "L7 (Loop 3)"}


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_entity = {row["entity"]: row for row in rows}
    return [by_entity[name] for name in ORDER]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(__file__).resolve().parent / "pcsk9_1459D05_loop_interface_summary.tsv",
    )
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    rows = read_summary(args.summary)
    labels = [LABELS[row["entity"]] for row in rows]
    delta_asa = [float(row["deltaASA_A2_rounded_1dp"]) for row in rows]
    atom_pairs = [int(row["heavy_atom_atom_pair_count_4p5A"]) for row in rows]

    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    bars = ax.barh(labels, delta_asa, color=["#fde68a", "#d1d5db", "#fb923c"])
    ax.invert_yaxis()
    ax.set_xlabel("ΔASA (Å²; Shrake–Rupley, 960 points, 1.4 Å probe)")
    ax.set_title("PCSK9/4OV6 Adnectin interface burial (4.5 Å atom-pair contacts)")
    for bar, value, contacts in zip(bars, delta_asa, atom_pairs):
        ax.text(
            bar.get_width() + 8,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.1f} Å²; {contacts} atom-pair contacts",
            va="center",
            fontsize=8,
        )
    ax.set_xlim(0, max(delta_asa) * 1.45)
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)
    fig.tight_layout()

    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / "Figure_S10D_PCSK9_4OV6_deltaASA.png"
    svg = args.outdir / "Figure_S10D_PCSK9_4OV6_deltaASA.svg"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print(png)
    print(svg)


if __name__ == "__main__":
    main()
