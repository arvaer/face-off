# Paper 1: Martini Water Neighbor-List Artifact

Reproduce the paper's claim that inadequate neighbor-list settings in Martini water with antifreeze create a statistically detectable anisotropy in the apparent pressure tensor, and condition-dependent box-shape drift under anisotropic or semi-isotropic pressure coupling, relative to safe-list controls.

This setup is ensemble-first. It is built to compare distributions across replicate trajectories, not to hunt for a single dramatic run.

## What This Uses

- Base coordinates: `../martini/sim/em.gro`
- Base Martini 2.2 force field: `../martini/topology/martini_v2.2.itp`
- Standard water bead definition: `../martini/topology/w.itp`

The paper's small solvent system uses a `6 x 6 x 6 nm` Martini water box with a `9:1` ratio of water to antifreeze beads (`W:WF`). This directory adds the missing `WF` bead definition and the minimal nonbonded extension needed to run that mixture locally.

## Layout

- `topology/`
- `templates/`
- `scripts/build_antifreeze_coordinates.py`
- `scripts/prepare_runs.py`
- `scripts/analyze_replicates.py`
- `run_all.sh`

## Conditions

- `artifact_nvt`: default-like problematic cutoff handling for pressure-tensor bias measurement.
- `safe_nvt`: widened safe list with dual list disabled.
- `artifact_aniso_pr`: anisotropic PR barostat with problematic list handling.
- `safe_aniso_pr`: anisotropic PR barostat with `nstlist = 1` and wide cutoff.
- `artifact_semi_pr`: semi-isotropic PR barostat with problematic list handling.
- `safe_semi_ber`: semi-isotropic Berendsen control with `nstlist = 1` and wide cutoff.

`NPT` preparations also cycle through:

- scale factors: `0.99`, `1.00`, `1.01`
- orientations: `none`, `x90`, `y90`, `z90`

That mirrors the paper's attempt to avoid interpreting one initial condition as an effect.

## Recommended Workflow

Prepare a 24-replicate ensemble:

```bash
python3 scripts/prepare_runs.py --replicates 24
```

Run everything:

```bash
bash run_all.sh
```

Run only the anisotropic PR comparison:

```bash
bash run_all.sh artifact_aniso_pr safe_aniso_pr
```

Analyze completed runs:

```bash
python3 scripts/analyze_replicates.py
```

Outputs land in:

- `runs/<condition>/replicate-XXXX/`
- `analysis/per_run.csv`
- `analysis/summary.json`

## Notes

- `prepare_runs.py` writes per-run `start.gro` files with the `9:1` `W:WF` mixture.
- `analyze_replicates.py` extracts `Pres-XX`, `Pres-YY`, `Pres-ZZ`, `Box-X`, `Box-Y`, and `Box-Z` from `run.edr` when needed.
- The primary pressure observable is

```text
Delta P = Pzz - (Pxx + Pyy) / 2
```

- The primary NPT observable is condition-dependent final box-shape drift, not a screenshot.
