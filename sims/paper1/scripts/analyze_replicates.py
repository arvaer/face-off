#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
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
    clean = [v for v in values if v == v and v not in (float("inf"), float("-inf"))]
    return statistics.pstdev(clean) if len(clean) > 1 else 0.0


def parse_runtime_settings(run_dir: Path, requested_nstlist: int):
    out = {
        "effective_nstlist": requested_nstlist,
        "effective_rlist": float("nan"),
        "dual_pair_list": 0,
        "failure_reason": "",
        "status": "missing",
    }
    log_path = run_dir / "run.log"
    if not log_path.exists():
        return out

    text = log_path.read_text(errors="replace")
    out["status"] = "started"
    changed = re.search(
        r"Changing nstlist from \d+ to (\d+), rlist from [0-9.]+ to ([0-9.]+)",
        text,
    )
    if changed:
        out["effective_nstlist"] = int(changed.group(1))
        out["effective_rlist"] = float(changed.group(2))

    outer = re.search(r"outer list:\s+updated every\s+(\d+)\s+steps,.*rlist\s+([0-9.]+)", text)
    if outer:
        out["effective_nstlist"] = int(outer.group(1))
        out["effective_rlist"] = float(outer.group(2))
    if "inner list:" in text:
        out["dual_pair_list"] = 1

    if "Fatal error:" in text:
        out["status"] = "failed"
        after = text.split("Fatal error:", 1)[1].strip().splitlines()
        if after:
            out["failure_reason"] = after[0].strip()
    elif "Writing final coordinates." in text or "Finished mdrun" in text:
        out["status"] = "complete"
    return out


def cycle_end_minus_start(values, nstlist: int):
    if nstlist <= 1 or len(values) < nstlist:
        return float("nan")
    diffs = []
    for start in range(0, len(values) - nstlist + 1, nstlist):
        block = values[start:start + nstlist]
        diffs.append(block[-1] - block[0])
    return mean(diffs)


def cycle_position_means(values, nstlist: int):
    if nstlist <= 0 or len(values) < nstlist:
        return []
    blocks = []
    for start in range(0, len(values) - nstlist + 1, nstlist):
        blocks.append(values[start:start + nstlist])
    if not blocks:
        return []
    return [mean([block[idx] for block in blocks]) for idx in range(nstlist)]


def summarize_condition(rows):
    complete_rows = [row for row in rows if row["status"] == "complete"]
    summary = {
        "runs": len(complete_rows),
        "total_records": len(rows),
        "failed_runs": sum(row["status"] == "failed" for row in rows),
        "missing_runs": sum(row["status"] == "missing" for row in rows),
        "started_runs": sum(row["status"] == "started" for row in rows),
        "delta_p_mean_bar_mean": mean([row["delta_p_mean_bar"] for row in complete_rows]),
        "delta_p_mean_bar_std": stddev([row["delta_p_mean_bar"] for row in complete_rows]),
        "delta_p_abs_mean_bar_mean": mean([row["delta_p_abs_mean_bar"] for row in complete_rows]),
    }
    if any(row["ensemble"] == "nvt" for row in complete_rows):
        summary["p_parallel_cycle_end_minus_start_bar_mean"] = mean(
            [row["p_parallel_cycle_end_minus_start_bar"] for row in complete_rows]
        )
        summary["p_perp_cycle_end_minus_start_bar_mean"] = mean(
            [row["p_perp_cycle_end_minus_start_bar"] for row in complete_rows]
        )
        summary["delta_p_cycle_end_minus_start_bar_mean"] = mean(
            [row["delta_p_cycle_end_minus_start_bar"] for row in complete_rows]
        )
    if any(row["ensemble"] == "npt_aniso" for row in complete_rows):
        axes = {"x": 0, "y": 0, "z": 0}
        for row in complete_rows:
            axis = row.get("dominant_axis")
            if axis in axes:
                axes[axis] += 1
        summary["dominant_axis_counts"] = axes
    if any(row["ensemble"] == "npt_semi" for row in complete_rows):
        flags = [row["elongated_z"] for row in complete_rows if row["elongated_z"] in (0, 1)]
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
            out = {
                **record,
                "status": "missing",
                "effective_nstlist": record["nstlist"],
                "effective_rlist": float("nan"),
                "dual_pair_list": 0,
                "failure_reason": "",
                "delta_p_mean_bar": float("nan"),
                "delta_p_abs_mean_bar": float("nan"),
                "p_parallel_cycle_end_minus_start_bar": float("nan"),
                "p_perp_cycle_end_minus_start_bar": float("nan"),
                "delta_p_cycle_end_minus_start_bar": float("nan"),
                "p_parallel_cycle_start_bar": float("nan"),
                "p_parallel_cycle_end_bar": float("nan"),
                "p_perp_cycle_start_bar": float("nan"),
                "p_perp_cycle_end_bar": float("nan"),
                "delta_p_cycle_start_bar": float("nan"),
                "delta_p_cycle_end_bar": float("nan"),
                "dominant_axis": "",
                "elongated_z": "",
                "final_shape_metric": float("nan"),
            }
            runtime = parse_runtime_settings(run_dir, record["nstlist"])
            out["status"] = runtime["status"]
            out["effective_nstlist"] = runtime["effective_nstlist"]
            out["effective_rlist"] = runtime["effective_rlist"]
            out["dual_pair_list"] = runtime["dual_pair_list"]
            out["failure_reason"] = runtime["failure_reason"]
            if edr.exists():
                try:
                    pressure = run_dir / "pressure.xvg"
                    maybe_extract(args.gmx, edr, pressure, ["Pres-XX", "Pres-YY", "Pres-ZZ"])
                    p_rows = load_xvg(pressure)
                    p_parallel = [0.5 * (row[1] + row[2]) for row in p_rows]
                    p_perp = [row[3] for row in p_rows]
                    delta = [perp - parallel for parallel, perp in zip(p_parallel, p_perp)]
                    cycle_nstlist = runtime["effective_nstlist"]
                    p_parallel_profile = cycle_position_means(p_parallel, cycle_nstlist)
                    p_perp_profile = cycle_position_means(p_perp, cycle_nstlist)
                    delta_profile = cycle_position_means(delta, cycle_nstlist)
                    out["status"] = "complete"
                    out["delta_p_mean_bar"] = mean(delta)
                    out["delta_p_abs_mean_bar"] = mean([abs(value) for value in delta])
                    out["p_parallel_cycle_end_minus_start_bar"] = cycle_end_minus_start(p_parallel, cycle_nstlist)
                    out["p_perp_cycle_end_minus_start_bar"] = cycle_end_minus_start(p_perp, cycle_nstlist)
                    out["delta_p_cycle_end_minus_start_bar"] = cycle_end_minus_start(delta, cycle_nstlist)
                    out["p_parallel_cycle_start_bar"] = p_parallel_profile[0] if p_parallel_profile else float("nan")
                    out["p_parallel_cycle_end_bar"] = p_parallel_profile[-1] if p_parallel_profile else float("nan")
                    out["p_perp_cycle_start_bar"] = p_perp_profile[0] if p_perp_profile else float("nan")
                    out["p_perp_cycle_end_bar"] = p_perp_profile[-1] if p_perp_profile else float("nan")
                    out["delta_p_cycle_start_bar"] = delta_profile[0] if delta_profile else float("nan")
                    out["delta_p_cycle_end_bar"] = delta_profile[-1] if delta_profile else float("nan")
                    if record["ensemble"] != "nvt":
                        box = run_dir / "box.xvg"
                        maybe_extract(args.gmx, edr, box, ["Box-X", "Box-Y", "Box-Z"])
                        b_rows = load_xvg(box)
                        x0, y0, z0 = b_rows[0][1], b_rows[0][2], b_rows[0][3]
                        xf, yf, zf = b_rows[-1][1], b_rows[-1][2], b_rows[-1][3]
                        if record["ensemble"] == "npt_aniso":
                            dominant = max(
                                {"x": xf / x0, "y": yf / y0, "z": zf / z0},
                                key=lambda key: {"x": xf / x0, "y": yf / y0, "z": zf / z0}[key],
                            )
                            out["dominant_axis"] = dominant
                            out["final_shape_metric"] = max(xf / x0, yf / y0, zf / z0)
                        else:
                            out["elongated_z"] = int((zf / z0) > ((xf / x0 + yf / y0) / 2.0))
                            out["final_shape_metric"] = (zf / z0) / ((xf / x0 + yf / y0) / 2.0)
                except Exception as exc:
                    out["status"] = "failed"
                    out["failure_reason"] = out["failure_reason"] or str(exc)
            rows.append(out)

    fieldnames = [
        "condition",
        "replicate",
        "ensemble",
        "status",
        "failure_reason",
        "seed",
        "nsteps",
        "nstlist",
        "effective_nstlist",
        "effective_rlist",
        "dual_pair_list",
        "scale",
        "rotation",
        "delta_p_mean_bar",
        "delta_p_abs_mean_bar",
        "p_parallel_cycle_end_minus_start_bar",
        "p_perp_cycle_end_minus_start_bar",
        "delta_p_cycle_end_minus_start_bar",
        "p_parallel_cycle_start_bar",
        "p_parallel_cycle_end_bar",
        "p_perp_cycle_start_bar",
        "p_perp_cycle_end_bar",
        "delta_p_cycle_start_bar",
        "delta_p_cycle_end_bar",
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
