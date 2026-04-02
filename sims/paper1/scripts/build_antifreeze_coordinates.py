#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path

TOTAL_BEADS = 1530
WF_COUNT = 153


def parse_gro(path: Path):
    lines = path.read_text().splitlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    atom_lines = lines[2:2 + natoms]
    box_line = lines[2 + natoms]
    atoms = []
    for line in atom_lines:
        resid = int(line[:5])
        resname = line[5:10]
        atomname = line[10:15]
        atomnr = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        tail = line[44:]
        atoms.append(
            {
                "resid": resid,
                "resname": resname,
                "atomname": atomname,
                "atomnr": atomnr,
                "xyz": [x, y, z],
                "tail": tail,
            }
        )
    box = [float(field) for field in box_line.split()]
    return title, atoms, box


def format_atom(atom):
    x, y, z = atom["xyz"]
    return (
        f'{atom["resid"]:5d}{atom["resname"]:<5}{atom["atomname"]:>5}'
        f'{atom["atomnr"]:5d}{x:8.3f}{y:8.3f}{z:8.3f}{atom["tail"]}'
    )


def wf_indices(count: int) -> set[int]:
    step = TOTAL_BEADS / count
    return {min(TOTAL_BEADS - 1, int(round(i * step))) for i in range(count)}


def rotate(vec, op: str):
    x, y, z = vec
    if op == "none":
        return [x, y, z]
    if op == "x90":
        return [x, -z, y]
    if op == "y90":
        return [z, y, -x]
    if op == "z90":
        return [-y, x, z]
    raise ValueError(f"unknown rotation {op}")


def transform_atoms(atoms, box, scale: float, rotation: str):
    center = [dim / 2.0 for dim in box[:3]]
    new_box = [dim * scale for dim in box[:3]]
    waters = []
    antifreeze = []
    picked = wf_indices(WF_COUNT)
    for index, atom in enumerate(atoms):
        shifted = [coord - ctr for coord, ctr in zip(atom["xyz"], center)]
        rotated = rotate(shifted, rotation)
        scaled = [coord * scale for coord in rotated]
        wrapped = []
        for coord, dim in zip(scaled, new_box):
            shifted_coord = coord + dim / 2.0
            wrapped.append(shifted_coord % dim)
        updated = dict(atom)
        updated["xyz"] = wrapped
        if index in picked:
            updated["resname"] = "WF"
            updated["atomname"] = "WF"
            antifreeze.append(updated)
        else:
            updated["resname"] = "W"
            updated["atomname"] = "W"
            waters.append(updated)
    out = waters + antifreeze
    for atomnr, atom in enumerate(out, start=1):
        atom["resid"] = atomnr
        atom["atomnr"] = atomnr
    return out, new_box


def write_gro(path: Path, title: str, atoms, box):
    lines = [title, f"{len(atoms):5d}"]
    lines.extend(format_atom(atom) for atom in atoms)
    lines.append(" ".join(f"{value:10.5f}" for value in box))
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("../martini/sim/em.gro"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scale", type=float, default=1.0)
    parser.add_argument(
        "--rotation",
        choices=["none", "x90", "y90", "z90"],
        default="none",
    )
    args = parser.parse_args()

    title, atoms, box = parse_gro(args.input)
    if len(atoms) != TOTAL_BEADS:
        raise SystemExit(f"expected {TOTAL_BEADS} beads, found {len(atoms)}")
    transformed, new_box = transform_atoms(atoms, box, args.scale, args.rotation)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_gro(
        args.output,
        f"{title} | W:WF = 9:1 | scale={args.scale:.2f} | rot={args.rotation}",
        transformed,
        new_box,
    )


if __name__ == "__main__":
    main()
