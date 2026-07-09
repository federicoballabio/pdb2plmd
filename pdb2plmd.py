#!/usr/bin/env python3
"""
pdb2plmd.py

Convert a general PDB extracted from a simulation into a PLUMED SAXS.cpp ONEBEAD-compatible template PDB,
while preserving the atom order of the selected simulation atoms.

Key design choice:
  The selected atom order is NEVER sorted or rebuilt. Output atoms appear in the same order as the
  ATOM/HETATM records selected from the input PDB. This is essential because PLUMED maps the ATOMS=
  list to TEMPLATE atoms by order.

Required options:
  -i / --input    input PDB
  -o / --output   output PDB

Optional:
  -a / --atoms    atom/order range to keep, e.g. "1-1062" or "1-100,150,200-250". Default: all.
                  The range refers to the 1-based ATOM/HETATM order in the input PDB, not necessarily
                  the PDB serial number.
  -g / --log      write a verbose log. If used without a filename, writes <output>.log

What it does:
  * preserves selected atom order;
  * renumbers output atom serials sequentially from 1;
  * assigns chain IDs if missing, preferentially from CHARMM/CHARMM-GUI segid such as RNAA/PROA;
  * inserts TER records between inferred chains/breaks;
  * renumbers residues sequentially within each output chain to avoid residue-number gaps;
  * converts common nucleic-acid residue names to SAXS.cpp/AMBER-style ONEBEAD names:
      ADE/RA -> A, CYT/RC -> C, GUA/RG -> G, URA/RU -> U;
      DADE/DA -> DA, DCYT/DC -> DC, DGUA/DG -> DG, THY/DT -> DT;
  * adds terminal RNA/DNA suffixes 5/3/T when inferable from atom content:
      C5/U5/A5/G5 for 5'-OH residues; C3/U3/A3/G3 for 3'-OH residues;
      CT/UT/AT/GT for 5'-phosphate with terminal OP3/O3P/HOP3/HP.
  * normalizes common nucleic-acid atom-name variants: * -> ', O1P/O2P/O3P kept as accepted by SAXS.cpp.

Caveats:
  The script cannot know the actual GROMACS atom order unless the input PDB was extracted from the same
  simulation/index group in that order. Use this on a PDB produced from the same TPR/index selection used
  in PLUMED ATOMS=.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass, replace
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

RNA_MAP = {
    "A": "A", "ADE": "A", "RA": "A", "RAD": "A", "R_A": "A",
    "C": "C", "CYT": "C", "RC": "C", "RCY": "C", "R_C": "C",
    "G": "G", "GUA": "G", "RG": "G", "RGU": "G", "R_G": "G",
    "U": "U", "URA": "U", "RU": "U", "URI": "U", "R_U": "U",
}
DNA_MAP = {
    "DA": "DA", "ADE_D": "DA", "DADE": "DA", "DAD": "DA",
    "DC": "DC", "CYT_D": "DC", "DCYT": "DC", "DCD": "DC",
    "DG": "DG", "GUA_D": "DG", "DGUA": "DG", "DGD": "DG",
    "DT": "DT", "THY": "DT", "DTH": "DT", "DTHY": "DT",
}
PROTEIN_NAMES = {
    "ALA","ARG","ASN","ASP","CYS","CYX","GLN","GLU","GLY","HIS","HID","HIE","HIP",
    "HSD","HSE","HSP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}
# Residues that usually should be passed through as small molecules/ions if selected.
PASS_THROUGH_RESNAMES = {
    "NA","K","CL","CLA","SOD","POT","MG","CA","ZN","ZN2","MN","FE","CU","CO","CD","NI",
    "FAD","FMN","HEM","HOH","WAT","TIP3","SOL"
}

PENTOSE_ATOMS = {
    "O5'","C5'","O4'","C4'","O3'","C3'","O2'","C2'","C1'",
    "H5'","H5''","H4'","H3'","H2'","H2''","H2'2","H1'",
    "HO5'","HO3'","HO2'","H5'1","H5'2","HO'2","H2'1","H5T","H3T",
}
BASE_ATOMS = {
    "N1","N2","N3","N4","N6","N7","N9","C2","C4","C5","C6","C7","C8",
    "O2","O4","O6","H1","H2","H3","H5","H6","H8","H21","H22","H41","H42",
    "H61","H62","H71","H72","H73",
}
PHOSPHATE_ATOMS = {"P","OP1","OP2","OP3","O1P","O2P","O3P","HP","HOP3"}
KNOWN_ONEBEAD_NUC_ATOMS = PENTOSE_ATOMS | BASE_ATOMS | PHOSPHATE_ATOMS

@dataclass
class AtomRecord:
    record: str
    input_atom_index: int       # 1-based order among ATOM/HETATM records
    input_serial: int
    atom_name: str
    altloc: str
    resname_orig: str
    chain_orig: str
    resseq_orig: int
    icode_orig: str
    x: float
    y: float
    z: float
    occ: float
    bfac: float
    segid: str
    element: str
    charge: str
    ter_before: bool = False

    # Filled during conversion
    chain_out: str = ""
    resseq_out: int = 0
    resname_out: str = ""
    atom_name_out: str = ""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare a PDB extracted from a simulation for PLUMED SAXS.cpp ONEBEAD, preserving atom order."
    )
    p.add_argument("-i", "--input", required=True, help="Input PDB extracted from the simulation/TPR selection.")
    p.add_argument("-o", "--output", required=True, help="Output SAXS.cpp-compatible template PDB.")
    p.add_argument("-a", "--atoms", default="all",
                   help="1-based ATOM/HETATM order range to keep, e.g. '1-1062' or '1-100,150,200-250'. Default: all.")
    p.add_argument("-g", "--log", nargs="?", const="__AUTO__", default=None,
                   help="Write verbose log. Optional filename; if omitted, writes <output>.log.")
    return p.parse_args()


def parse_range(expr: str, n_atoms: int) -> List[int]:
    expr = (expr or "all").strip().lower()
    if expr in {"all", "*"}:
        return list(range(1, n_atoms + 1))
    selected: List[int] = []
    seen = set()
    for part in expr.split(','):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise SystemExit(f"Invalid atom range component: {part!r}")
        a = int(m.group(1))
        b = int(m.group(2)) if m.group(2) else a
        if a < 1 or b < 1 or b < a:
            raise SystemExit(f"Invalid atom range component: {part!r}")
        if b > n_atoms:
            raise SystemExit(f"Atom range {part!r} exceeds number of ATOM/HETATM records ({n_atoms}).")
        for idx in range(a, b + 1):
            if idx not in seen:
                selected.append(idx)
                seen.add(idx)
    return selected


def safe_int(s: str, default: int = 0) -> int:
    try:
        return int(s.strip())
    except Exception:
        return default


def safe_float(s: str, default: float = 0.0) -> float:
    try:
        return float(s.strip())
    except Exception:
        return default


def infer_element(atom_name: str, element_field: str = "") -> str:
    e = element_field.strip().upper()
    if e:
        return e[:2]
    name = atom_name.strip().replace("'", "").replace("*", "")
    if not name:
        return ""
    # For names such as 1H5, H5', C1', OP1, CL, NA.
    if name[0].isdigit() and len(name) > 1:
        name = name[1:]
    # Two-letter ions/elements when atom name is exactly a two-letter element.
    up = name.upper()
    if up.startswith("CL"):
        return "CL"
    if up.startswith("NA"):
        return "NA"
    if up.startswith("MG"):
        return "MG"
    if up.startswith("ZN"):
        return "ZN"
    if up.startswith("CA") and len(up) <= 2:
        return "CA"
    return up[0]


def parse_pdb(path: str) -> List[AtomRecord]:
    atoms: List[AtomRecord] = []
    ter_pending = False
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "TER":
                ter_pending = True
                continue
            if rec not in {"ATOM", "HETATM"}:
                continue
            padded = line.rstrip('\n')
            padded = padded + " " * max(0, 80 - len(padded))
            atom_name = padded[12:16].strip()
            element = padded[76:78].strip()
            atom = AtomRecord(
                record=rec,
                input_atom_index=len(atoms) + 1,
                input_serial=safe_int(padded[6:11], len(atoms) + 1),
                atom_name=atom_name,
                altloc=padded[16:17].strip(),
                resname_orig=padded[17:20].strip(),
                chain_orig=padded[21:22].strip(),
                resseq_orig=safe_int(padded[22:26], 0),
                icode_orig=padded[26:27].strip(),
                x=safe_float(padded[30:38]),
                y=safe_float(padded[38:46]),
                z=safe_float(padded[46:54]),
                occ=safe_float(padded[54:60], 1.0),
                bfac=safe_float(padded[60:66], 0.0),
                segid=padded[72:76].strip(),
                element=infer_element(atom_name, element),
                charge=padded[78:80].strip(),
                ter_before=ter_pending,
            )
            atoms.append(atom)
            ter_pending = False
    return atoms


def normalize_atom_name(name: str) -> str:
    # SAXS.cpp recognizes apostrophe names. Convert old PDB '*' notation to apostrophe.
    name = name.strip().replace('*', "'")
    # Common aliases from some builders; keep O1P/O2P too because SAXS.cpp accepts them.
    aliases = {
        "O1P": "O1P", "O2P": "O2P", "O3P": "O3P",
        "OP1": "OP1", "OP2": "OP2", "OP3": "OP3",
        "H5'1": "H5'1", "H5'2": "H5'2", "H2'1": "H2'1", "H2'2": "H2'2",
    }
    return aliases.get(name, name)


def base_resname(resname: str) -> str:
    r = resname.strip().upper()
    # Remove common CHARMM terminal prefixes/suffixes without losing explicit PLUMED suffixes.
    r = r.replace("5'", "").replace("3'", "")
    if r in RNA_MAP:
        return RNA_MAP[r]
    if r in DNA_MAP:
        return DNA_MAP[r]
    return r


def is_nucleic_base_name(r: str) -> bool:
    b = base_resname(r)
    return b in {"A","C","G","U","DA","DC","DG","DT"}


def is_rna_base(b: str) -> bool:
    return b in {"A", "C", "G", "U"}


def is_dna_base(b: str) -> bool:
    return b in {"DA", "DC", "DG", "DT"}


def derive_chain_from_segid(segid: str) -> str:
    s = segid.strip()
    if not s:
        return ""
    # CHARMM-GUI frequently uses PROA/RNAA/PROB/RNAB. Last alphanumeric char is usually the chain.
    for ch in reversed(s):
        if ch.isalnum():
            return ch
    return ""


def next_chain_id(used: set) -> str:
    for c in CHAIN_IDS:
        if c not in used:
            used.add(c)
            return c
    raise SystemExit("Too many inferred chains for single-character PDB chain IDs (>62).")


def residue_key(a: AtomRecord) -> Tuple[str, str, int, str, str]:
    return (a.chain_orig, a.segid, a.resseq_orig, a.icode_orig, a.resname_orig)


def group_residues(atoms: Sequence[AtomRecord]) -> List[List[AtomRecord]]:
    residues: List[List[AtomRecord]] = []
    current: List[AtomRecord] = []
    last_key = None
    for a in atoms:
        k = residue_key(a)
        if current and k != last_key:
            residues.append(current)
            current = []
        current.append(a)
        last_key = k
    if current:
        residues.append(current)
    return residues


def residue_atom_names(res_atoms: Sequence[AtomRecord]) -> set:
    return {normalize_atom_name(a.atom_name) for a in res_atoms}


def terminal_suffix(base: str, atom_names: set, is_first_in_chain: bool, is_last_in_chain: bool) -> str:
    if not (is_rna_base(base) or is_dna_base(base)):
        return ""
    # 5'-phosphate terminal with an extra terminal phosphate oxygen/hydrogen.
    if is_first_in_chain and ("OP3" in atom_names or "O3P" in atom_names or "HOP3" in atom_names or "HP" in atom_names):
        return "T"
    # 5'-OH terminal: no phosphate in the residue and terminal sugar H/O marker.
    if is_first_in_chain and not (atom_names & PHOSPHATE_ATOMS):
        if {"H5T", "HO5'"} & atom_names or "O5'" in atom_names:
            return "5"
    # 3'-OH terminal marker.
    if is_last_in_chain and ({"H3T", "HO3'"} & atom_names):
        return "3"
    return ""


def split_into_chains(residues: Sequence[Sequence[AtomRecord]]) -> List[List[List[AtomRecord]]]:
    chains: List[List[List[AtomRecord]]] = []
    current: List[List[AtomRecord]] = []
    prev_first: Optional[AtomRecord] = None
    for res in residues:
        first = res[0]
        new_chain = False
        if not current:
            new_chain = True
        else:
            assert prev_first is not None
            # Explicit TER before this residue.
            if first.ter_before:
                new_chain = True
            # Explicit original chain or segid change.
            elif first.chain_orig and prev_first.chain_orig and first.chain_orig != prev_first.chain_orig:
                new_chain = True
            elif first.segid and prev_first.segid and first.segid != prev_first.segid:
                new_chain = True
            # Missing/ambiguous chain: residue numbering reset or gap. This avoids residue-range gaps in PLUMED.
            elif first.resseq_orig <= prev_first.resseq_orig:
                new_chain = True
            elif first.resseq_orig > prev_first.resseq_orig + 1:
                new_chain = True
        if new_chain:
            if current:
                chains.append(current)
            current = [list(res)]
        else:
            current.append(list(res))
        prev_first = first
    if current:
        chains.append(current)
    return chains


def convert_atoms(atoms: List[AtomRecord]) -> Tuple[List[AtomRecord], List[str]]:
    log: List[str] = []
    residues = group_residues(atoms)
    chains = split_into_chains(residues)
    used_chains = set()
    converted: List[AtomRecord] = []
    log.append(f"Input selected atoms: {len(atoms)}")
    log.append(f"Input selected residues: {len(residues)}")
    log.append(f"Inferred output chains: {len(chains)}")

    for chain_idx, chain_residues in enumerate(chains, start=1):
        first_atom = chain_residues[0][0]
        chain_id = first_atom.chain_orig or derive_chain_from_segid(first_atom.segid)
        if not chain_id or chain_id in used_chains:
            chain_id = next_chain_id(used_chains)
        else:
            used_chains.add(chain_id)
        log.append(f"Chain {chain_idx}: output chain ID {chain_id!r}; residues={len(chain_residues)}; atoms={sum(len(r) for r in chain_residues)}; source first res={first_atom.resname_orig}{first_atom.resseq_orig} segid={first_atom.segid!r}")

        for r_idx, res_atoms in enumerate(chain_residues, start=1):
            orig = res_atoms[0]
            base = base_resname(orig.resname_orig)
            names = residue_atom_names(res_atoms)
            suffix = terminal_suffix(base, names, r_idx == 1, r_idx == len(chain_residues))
            resname_out = base + suffix if (is_rna_base(base) or is_dna_base(base)) else base
            if resname_out == orig.resname_orig.strip():
                log_res_conversion = "kept"
            else:
                log_res_conversion = f"{orig.resname_orig.strip()} -> {resname_out}"
            if is_nucleic_base_name(orig.resname_orig):
                unknown = sorted(n for n in names if n not in KNOWN_ONEBEAD_NUC_ATOMS)
                if unknown:
                    log.append(f"WARNING residue {orig.resname_orig}{orig.resseq_orig} chain {chain_id}: atom names not recognized by SAXS.cpp ONEBEAD nucleic-acid mapping: {','.join(unknown)}")
            elif base not in PROTEIN_NAMES and base not in PASS_THROUGH_RESNAMES:
                log.append(f"WARNING residue {orig.resname_orig}{orig.resseq_orig} chain {chain_id}: residue name after conversion is {base!r}; check SAXS.cpp support if selected.")
            if suffix:
                log.append(f"Terminal inference: chain {chain_id} residue input {orig.resname_orig}{orig.resseq_orig} -> {resname_out} ({'first' if r_idx == 1 else 'last'} residue).")
            elif log_res_conversion != "kept":
                log.append(f"Residue conversion: chain {chain_id} residue input {orig.resname_orig}{orig.resseq_orig} -> {resname_out}.")

            for a in res_atoms:
                converted.append(replace(
                    a,
                    chain_out=chain_id,
                    resseq_out=r_idx,
                    resname_out=resname_out,
                    atom_name_out=normalize_atom_name(a.atom_name),
                ))
    return converted, log


def format_atom_name(name: str, element: str) -> str:
    n = name.strip()
    if len(n) >= 4:
        return n[:4]
    # PDB convention: one-letter elements are right-justified in atom-name field.
    if len(element.strip()) <= 1 and n and not n[0].isdigit():
        return f"{n:>4}"
    return f"{n:<4}"


def format_pdb_atom(a: AtomRecord, out_serial: int) -> str:
    rec = a.record if a.record in {"ATOM", "HETATM"} else "ATOM"
    atom_field = format_atom_name(a.atom_name_out, a.element)
    res_field = f"{a.resname_out:>3}"[-3:]
    chain = (a.chain_out or "A")[:1]
    element = (a.element or infer_element(a.atom_name_out)).upper()[:2]
    # Keep segid in columns 73-76 when available; otherwise use chain-friendly placeholder.
    segid = (a.segid or ("RNA" + chain if is_nucleic_base_name(a.resname_out) else chain))[:4]
    return (f"{rec:<6}{out_serial:5d} {atom_field}{a.altloc[:1]:1s}{res_field:>3s} {chain:1s}"
            f"{a.resseq_out:4d}{a.icode_orig[:1]:1s}   "
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occ:6.2f}{a.bfac:6.2f}"
            f"      {segid:<4s}{element:>2s}{a.charge:>2s}")


def format_ter(serial: int, last_atom: AtomRecord) -> str:
    return f"TER   {serial:5d}      {last_atom.resname_out:>3s} {last_atom.chain_out[:1]:1s}{last_atom.resseq_out:4d}"


def write_pdb(atoms: List[AtomRecord], path: str) -> None:
    with open(path, 'w', encoding='utf-8') as out:
        out.write("REMARK Prepared by pdb2plmd.py for PLUMED SAXS.cpp ONEBEAD\n")
        out.write("REMARK Atom order preserved from selected ATOM/HETATM input order\n")
        prev_atom = None
        for i, a in enumerate(atoms, start=1):
            if prev_atom is not None and a.chain_out != prev_atom.chain_out:
                out.write(format_ter(i, prev_atom) + "\n")
            out.write(format_pdb_atom(a, i) + "\n")
            prev_atom = a
        if atoms:
            out.write(format_ter(len(atoms) + 1, atoms[-1]) + "\n")
        out.write("END\n")


def write_log(log_path: Optional[str], args: argparse.Namespace, all_atoms: List[AtomRecord], selected_indices: List[int], converted: List[AtomRecord], log_lines: List[str]) -> None:
    if not log_path:
        return
    with open(log_path, 'w', encoding='utf-8') as log:
        log.write("pdb2plmd.py verbose log\n")
        log.write("================================\n")
        log.write(f"Input:  {args.input}\n")
        log.write(f"Output: {args.output}\n")
        log.write(f"Atom range expression: {args.atoms}\n")
        log.write(f"Total ATOM/HETATM records in input: {len(all_atoms)}\n")
        log.write(f"Selected atoms: {len(selected_indices)}\n")
        if selected_indices:
            log.write(f"Selected input atom-order range: {selected_indices[0]}..{selected_indices[-1]}\n")
        log.write("\n")
        log.write("Notes:\n")
        log.write("- The -a range refers to input ATOM/HETATM order, not PDB serial.\n")
        log.write("- Output atom order is identical to the selected input atom order.\n")
        log.write("- Output atom serials are renumbered sequentially; this does not change PLUMED atom order.\n")
        log.write("- Residues are renumbered sequentially within each output chain to avoid PLUMED residue-range gaps.\n")
        log.write("\n")
        for line in log_lines:
            log.write(line + "\n")
        log.write("\nFirst 20 atom mapping rows:\n")
        log.write("out_serial\tinput_atom_order\tinput_serial\tchain\tresid\tresname\tatom\telement\n")
        for i, a in enumerate(converted[:20], start=1):
            log.write(f"{i}\t{a.input_atom_index}\t{a.input_serial}\t{a.chain_out}\t{a.resseq_out}\t{a.resname_out}\t{a.atom_name_out}\t{a.element}\n")
        if len(converted) > 20:
            log.write(f"... {len(converted)-20} more atoms not shown ...\n")


def main() -> None:
    args = parse_args()
    if args.log == "__AUTO__":
        args.log = args.output + ".log"

    all_atoms = parse_pdb(args.input)
    if not all_atoms:
        raise SystemExit(f"No ATOM/HETATM records found in {args.input}")
    selected_indices = parse_range(args.atoms, len(all_atoms))
    index_set = set(selected_indices)
    selected_atoms = [a for a in all_atoms if a.input_atom_index in index_set]

    # Preserve original atom order, regardless of range expression order.
    converted, log_lines = convert_atoms(selected_atoms)
    write_pdb(converted, args.output)
    write_log(args.log, args, all_atoms, selected_indices, converted, log_lines)

    print(f"Wrote {args.output} with {len(converted)} atoms.")
    if args.log:
        print(f"Wrote verbose log: {args.log}")
    print("Atom order was preserved from the selected input ATOM/HETATM order.")


if __name__ == "__main__":
    main()
