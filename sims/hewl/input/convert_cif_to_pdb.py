#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


ATOM_FIELDS = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert a simple PDBx/mmCIF atom table to PDB.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--chain", default="A")
    return parser.parse_args()


def tokenize(line: str) -> list[str]:
    return line.split()


def format_atom_name(atom_name: str, element: str) -> str:
    atom_name = atom_name.strip()
    if len(atom_name) == 4:
        return atom_name
    if len(element.strip()) == 1:
        return f" {atom_name:<3}"
    return f"{atom_name:<4}"


def main() -> None:
    args = parse_args()
    lines = args.input.read_text().splitlines()
    in_atom_loop = False
    field_names: list[str] = []
    records: list[str] = []

    for line in lines:
        if line.startswith("loop_"):
            in_atom_loop = False
            field_names = []
            continue
        if line.startswith("_atom_site."):
            in_atom_loop = True
            field_names.append(line.split(".", 1)[1].strip())
            continue
        if in_atom_loop and field_names == ATOM_FIELDS:
            if not line or line.startswith("#"):
                break
            values = tokenize(line)
            if len(values) != len(ATOM_FIELDS):
                continue
            row = dict(zip(ATOM_FIELDS, values))
            if row["group_PDB"] not in {"ATOM", "HETATM"}:
                continue
            if row["auth_asym_id"] != args.chain:
                continue
            alt_id = row["label_alt_id"]
            if alt_id not in {".", "?", "A"}:
                continue
            atom_serial = int(row["id"])
            atom_name = format_atom_name(row["auth_atom_id"], row["type_symbol"])
            residue_name = row["auth_comp_id"][:3]
            chain_id = row["auth_asym_id"][:1]
            residue_seq = int(row["auth_seq_id"])
            insertion_code = " " if row["pdbx_PDB_ins_code"] in {"?", "."} else row["pdbx_PDB_ins_code"][:1]
            x = float(row["Cartn_x"])
            y = float(row["Cartn_y"])
            z = float(row["Cartn_z"])
            occupancy = 1.0 if row["occupancy"] in {"?", "."} else float(row["occupancy"])
            bfactor = 0.0 if row["B_iso_or_equiv"] in {"?", "."} else float(row["B_iso_or_equiv"])
            element = row["type_symbol"].strip().upper()[:2]
            record = (
                f"{row['group_PDB']:<6}{atom_serial:>5} {atom_name}{' ':1}"
                f"{residue_name:>3} {chain_id:1}{residue_seq:>4}{insertion_code:1}   "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}          "
                f"{element:>2}"
            )
            records.append(record)

    if not records:
        raise SystemExit("No ATOM/HETATM records extracted from mmCIF")

    output_lines = [f"REMARK Converted from {args.input.name}", *records, "TER", "END"]
    args.output.write_text("\n".join(output_lines) + "\n")
    print(f"Wrote {len(records)} atoms to {args.output}")


if __name__ == "__main__":
    main()
