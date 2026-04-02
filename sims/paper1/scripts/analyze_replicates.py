#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "manifest.json"
ANALYSIS = ROOT / "analysis"


def load_xvg(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped[0] in "#@":
            continue
        rows.append([float(field) for field in stripped.split()])
    return rows


def maybe_extract(gmx: str, edr: Path, out_path: Path, labels):
    if out_path.exists():
        return
    proc = subprocess.run(
        [gmx, "energy", "-f", str(edr), "-o", str(out_path)],
        input="\n".join(labels) + "\n",
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"failed to extract {out_path.name}")


def mean(values):
    return statistics.fmean(values) if values else float("nan")


def stddev(values):
    return statistics.pstdev(values) if len(values) > 1 else 0.0


def cycle_jump(delta_p, nstlist: int):
    if nstlist <= 1 or len(delta_p) < nstlist:
        return float("nan")
    jumps = []
    for start in range(0, len(delta_p) - nstlist + 1, nstlist):
        block = delta_p[start:start + nstlist]
        jumps.append(block[-1] - block[0])
    return mean(jumps)


def summarize_condition(rows):
    summary = {
        "runs": len(rows),
        "delta_p_mean_bar_mean": mean([row["delta_p_mean_bar"] for row in rows]),
        "delta_p_mean_bar_std": stddev([row["delta_p_mean_bar"] for row in rows]),
        "delta_p_abs_mean_bar_mean": mean([row["delta_p_abs_mean_bar"] for row in rows]),
    }
    if any(row["ensemble"] == "npt_aniso" for row in rows):
        axes = {"x": 0, "y": 0, "z": 0}
        for row in rows:
            axis = row.get("dominant_axis")
            if axis in axes:
                axes[axis] += 1
        summary["dominant_axis_counts"] = axes
    if any(row["ensemble"] == "npt_semi" for row in rows):
        flags = [row["elongated_z"] for row in rows if row["elongated_z"] in (0, 1)]
        summary["elongated_z_fraction"] = mean(flags)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--gmx", default="gmx")
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    ANALYSIS.mkdir(parents=True, exist_ok=True)

    rows = []
    for condition, records in manifest["conditions"].items():
        for record in records:
            run_dir = ROOT / record["run_dir"]
            edr = run_dir / "run.edr"
            if not edr.exists():
                continue
            pressure = run_dir / "pressure.xvg"
            maybe_extract(args.gmx, edr, pressure, ["Pres-XX", "Pres-YY", "Pres-ZZ"])
            p_rows = load_xvg(pressure)
            delta = [row[3] - 0.5 * (row[1] + row[2]) for row in p_rows]
            out = {
                **record,
                "delta_p_mean_bar": mean(delta),
                "delta_p_abs_mean_bar": mean([abs(value) for value in delta]),
                "delta_p_cycle_jump_bar": cycle_jump(delta, record["nstlist"]),
                "dominant_axis": "",
                "elongated_z": "",
                "final_shape_metric": float("nan"),
            }
            if record["ensemble"] != "nvt":
                box = run_dir / "box.xvg"
                maybe_extract(args.gmx, edr, box, ["Box-X", "Box-Y", "Box-Z"])
                b_rows = load_xvg(box)
                x0, y0, z0 = b_rows[0][1], b_rows[0][2], b_rows[0][3]
                xf, yf, zf = b_rows[-1][1], b_rows[-1][2], b_rows[-1][3]
                if record["ensemble"] == "npt_aniso":
                    dominant = max({"x": xf / x0, "y": yf / y0, "z": zf / z0}, key=lambda key: {"x": xf / x0, "y": yf / y0, "z": zf / z0}[key])
                    out["dominant_axis"] = dominant
                    out["final_shape_metric"] = max(xf / x0, yf / y0, zf / z0)
                else:
                    out["elongated_z"] = int((zf / z0) > ((xf / x0 + yf / y0) / 2.0))
                    out["final_shape_metric"] = (zf / z0) / ((xf / x0 + yf / y0) / 2.0)
            rows.append(out)

    fieldnames = [
        "condition",
        "replicate",
        "ensemble",
        "seed",
        "nsteps",
        "nstlist",
        "scale",
        "rotation",
        "delta_p_mean_bar",
        "delta_p_abs_mean_bar",
        "delta_p_cycle_jump_bar",
        "dominant_axis",
        "elongated_z",
        "final_shape_metric",
        "run_dir",
    ]
    with (ANALYSIS / "per_run.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    by_condition = {}
    for condition in manifest["conditions"]:
        cond_rows = [row for row in rows if row["condition"] == condition]
        by_condition[condition] = summarize_condition(cond_rows)
    (ANALYSIS / "summary.json").write_text(json.dumps(by_condition, indent=2) + "\n")


if __name__ == "__main__":
    main()
