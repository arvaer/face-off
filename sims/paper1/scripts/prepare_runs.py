#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = ROOT / "templates"
RUNS = ROOT / "runs"

CONDITIONS = {
    "artifact_nvt": {
        "template": "artifact_nvt.mdp.in",
        "steps_arg": "nvt_steps",
        "ensemble": "nvt",
        "nstlist": 20,
    },
    "safe_nvt": {
        "template": "safe_nvt.mdp.in",
        "steps_arg": "nvt_steps",
        "ensemble": "nvt",
        "nstlist": 20,
    },
    "paper_nvt_rl128_nst25": {
        "template": "paper_nvt_rl128_nst25.mdp.in",
        "steps_arg": "nvt_steps",
        "ensemble": "nvt",
        "nstlist": 25,
    },
    "paper_nvt_rl1422_nst20": {
        "template": "safe_nvt.mdp.in",
        "steps_arg": "nvt_steps",
        "ensemble": "nvt",
        "nstlist": 20,
    },
    "artifact_aniso_pr": {
        "template": "artifact_aniso_pr.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_aniso",
        "nstlist": 20,
    },
    "safe_aniso_pr": {
        "template": "safe_aniso_pr.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_aniso",
        "nstlist": 1,
    },
    "paper_aniso_pr_rl128_nst25": {
        "template": "paper_aniso_pr_rl128_nst25.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_aniso",
        "nstlist": 25,
    },
    "paper_aniso_pr_rl19_nst20": {
        "template": "paper_aniso_pr_rl19_nst20.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_aniso",
        "nstlist": 20,
    },
    "paper_aniso_pr_rl19_nst1": {
        "template": "safe_aniso_pr.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_aniso",
        "nstlist": 1,
    },
    "artifact_semi_pr": {
        "template": "artifact_semi_pr.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_semi",
        "nstlist": 20,
    },
    "safe_semi_ber": {
        "template": "safe_semi_ber.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_semi",
        "nstlist": 1,
    },
    "paper_semi_pr_rl19_nst20": {
        "template": "paper_semi_pr_rl19_nst20.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_semi",
        "nstlist": 20,
    },
    "paper_semi_pr_rl19_nst1": {
        "template": "paper_semi_pr_rl19_nst1.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_semi",
        "nstlist": 1,
    },
    "paper_semi_ber_rl19_nst1": {
        "template": "safe_semi_ber.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_semi",
        "nstlist": 1,
    },
    "paper_semi_cres_rl19_nst1": {
        "template": "paper_semi_cres_rl19_nst1.mdp.in",
        "steps_arg": "npt_steps",
        "ensemble": "npt_semi",
        "nstlist": 1,
    },
}


def stable_seed(condition: str, replicate: int) -> int:
    digest = hashlib.sha256(f"{condition}:{replicate}".encode()).hexdigest()
    return 100000 + int(digest[:8], 16) % 900000000


def npt_variant(replicate: int):
    scales = [0.99, 1.00, 1.01]
    rotations = ["none", "x90", "y90", "z90"]
    scale = scales[(replicate - 1) % len(scales)]
    rotation = rotations[((replicate - 1) // len(scales)) % len(rotations)]
    return scale, rotation


def render_template(name: str, *, nsteps: int, seed: int) -> str:
    text = (TEMPLATES / name).read_text()
    return text.replace("__NSTEPS__", str(nsteps)).replace("__GEN_SEED__", str(seed))


def write_run_script(run_dir: Path):
    script = """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

GMX="${GMX:-gmx}"
ROOT="$(cd ../../.. && pwd)"
OMP_THREADS="${OMP_THREADS:-8}"
export OMP_NUM_THREADS="$OMP_THREADS"

"$GMX" grompp \
  -f run.mdp \
  -c start.gro \
  -p "$ROOT/topology/water_antifreeze.top" \
  -o run.tpr \
  -maxwarn 2

"$GMX" mdrun -deffnm run -ntmpi 1 -ntomp "$OMP_THREADS" -pin on
"""
    path = run_dir / "run.sh"
    path.write_text(script)
    path.chmod(0o755)


def build_start(run_dir: Path, scale: float, rotation: str):
    cmd = [
        "python3",
        str(ROOT / "scripts" / "build_antifreeze_coordinates.py"),
        "--input",
        str(ROOT / "../martini/sim/em.gro"),
        "--output",
        str(run_dir / "start.gro"),
        "--scale",
        f"{scale:.2f}",
        "--rotation",
        rotation,
    ]
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=24)
    parser.add_argument("--nvt-steps", type=int, default=500000)
    parser.add_argument("--npt-steps", type=int, default=2500000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    RUNS.mkdir(parents=True, exist_ok=True)
    manifest = {"replicates": args.replicates, "conditions": {}}

    for condition, meta in CONDITIONS.items():
        cond_dir = RUNS / condition
        if args.overwrite and cond_dir.exists():
            shutil.rmtree(cond_dir)
        cond_dir.mkdir(parents=True, exist_ok=True)
        manifest["conditions"][condition] = []
        for rep in range(1, args.replicates + 1):
            run_dir = cond_dir / f"replicate-{rep:04d}"
            if run_dir.exists() and not args.overwrite:
                if not (run_dir / "metadata.json").exists():
                    raise SystemExit(f"{run_dir} exists without metadata; use --overwrite")
            run_dir.mkdir(parents=True, exist_ok=True)
            scale, rotation = (1.0, "none")
            if meta["ensemble"] != "nvt":
                scale, rotation = npt_variant(rep)
            seed = stable_seed(condition, rep)
            nsteps = getattr(args, meta["steps_arg"])
            (run_dir / "run.mdp").write_text(
                render_template(meta["template"], nsteps=nsteps, seed=seed)
            )
            build_start(run_dir, scale, rotation)
            write_run_script(run_dir)
            record = {
                "condition": condition,
                "replicate": rep,
                "seed": seed,
                "nsteps": nsteps,
                "ensemble": meta["ensemble"],
                "nstlist": meta["nstlist"],
                "scale": scale,
                "rotation": rotation,
                "run_dir": str(run_dir.relative_to(ROOT)),
            }
            (run_dir / "metadata.json").write_text(json.dumps(record, indent=2) + "\n")
            manifest["conditions"][condition].append(record)

    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
