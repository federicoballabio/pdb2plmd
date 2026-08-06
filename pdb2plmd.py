#!/usr/bin/env python3
"""Prepare an ordered PDB template for PLUMED SAXS.cpp ONEBEAD."""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter
from dataclasses import dataclass, replace
from typing import List, Optional, Sequence, Tuple

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
ION_MAP = {
    "NA": "NA", "SOD": "NA", "NA+": "NA",
    "K": "K", "POT": "K", "K+": "K",
    "CL": "CL", "CLA": "CL", "CL-": "CL",
    "CA": "CAL", "CAL": "CAL", "CA2+": "CAL",
    "MG": "MG", "MG2+": "MG",
    "ZN": "ZN", "ZN2": "ZN", "ZN2+": "ZN",
    "FE2": "FE2", "FE3": "FE3",
    "MN": "MN", "MN2+": "MN",
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

# Glycan names are normalized to the residue names used by SAXS.cpp ONEBEAD.

_GLC_FAMILY = {"GLC", "BGC", "AGLC", "BGLC"}
_GAL_FAMILY = {"GAL", "GLA", "AGAL", "BGAL"}
_MAN_A      = {"MAN", "AMAN"}
_MAN_B      = {"BMA", "BMAN"}
_FUC_FAMILY = {"FUC", "FUL", "FCA", "FCB", "AFUC", "BFUC"}
_NEU_FAMILY = {"SIA", "SLB", "ANE5AC", "BNE5AC",
               "ANE5", "BNE5"}   # CHARMM names as truncated by the PDB writer
_NAG_NAMES  = {"NAG", "NDG", "AGLCNA", "BGLCNA"}
_NGA_NAMES  = {"NGA", "A2G", "AGALNA", "BGALNA"}

# recognised monosaccharides for which no ONEBEAD parameters exist
GLYCAN_UNSUPPORTED = {
    "RAM", "RM4", "XXR", "XYL", "XYS", "XYP", "LXZ", "GCU", "BDP", "GCV", "IDS", "IDR", "SGN", "SUS",
    "UAP", "NGC", "NGE", "RIB", "ARA", "ARB", "AHR", "GLP", "PA1", "GCS",
}

# Residues renamed by a builder to mark a glycosylation site. Only a hydrogen is
# lost to the glycosidic bond, so the bead is the parent amino acid.
GLYCOSYLATION_SITE_MAP = {
    "NLN": "ASN",
    "OLS": "SER",
    "OLT": "THR",
}
UNSUPPORTED_SITE_NAMES = {"OLP"}

# GLYCAM/CHARMM -> PDB chemical component atom names.
GLYCAN_ATOM_MAP_HEX = {
    "N": "N2", "C": "C7", "O": "O7", "CT": "C8",
    "HN": "HN2", "HT1": "H81", "HT2": "H82", "HT3": "H83",
}
GLYCAN_ATOM_MAP_SIA = {
    "N": "N5", "C": "C10", "O": "O10", "CT": "C11",
    "O11": "O1A", "O12": "O1B",
    "HN": "HN5", "HT1": "H111", "HT2": "H112", "HT3": "H113",
}


def canonical_glycan(resname, atom_names):
    """Canonical PDB-CCD code for a monosaccharide, or None."""
    r = resname.strip().upper()
    has_n = any(a[:1] == "N" for a in atom_names)
    if r in _NAG_NAMES:
        return "NAG"
    if r in _NGA_NAMES:
        return "NGA"
    if r in _GLC_FAMILY:
        return "NAG" if has_n else "GLC"
    if r in _GAL_FAMILY:
        return "NGA" if has_n else "GAL"
    if r in _MAN_A:
        return "MAN"
    if r in _MAN_B:
        return "BMA"
    if r in _FUC_FAMILY:
        return "FUC"
    if r in _NEU_FAMILY:
        return "SIA"
    return None


def normalize_glycan_atom_name(name, resname_out):
    n = name.strip()
    if resname_out == "SIA":
        return GLYCAN_ATOM_MAP_SIA.get(n, n)
    if resname_out in {"NAG", "NGA"}:
        return GLYCAN_ATOM_MAP_HEX.get(n, n)
    return n


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
    line_number: int = 0
    ter_before: bool = False

    # Filled during conversion
    chain_out: str = ""
    resseq_out: int = 0
    resname_out: str = ""
    atom_name_out: str = ""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare an ordered PDB template for PLUMED SAXS.cpp ONEBEAD."
    )
    p.add_argument("-i", "--input", required=True, help="Input PDB extracted from the simulation/TPR selection.")
    p.add_argument("-o", "--output", required=True, help="Output SAXS.cpp-compatible template PDB.")
    p.add_argument("-a", "--atoms", default="all",
                   help="1-based ATOM/HETATM order range to keep, e.g. '1-1062' or '1-100,150,200-250'. Default: all.")
    p.add_argument("--model", type=int, default=1,
                   help="MODEL serial to keep from a multi-model file. Default: 1.")
    p.add_argument("--charmm", dest="charmm", action="store_true", default=None,
                   help="Force CHARMM/CHARMM-GUI handling: read the four-character residue "
                        "name from columns 18-21 and take the chain from the segid in columns "
                        "73-76. Auto-detected when column 21 is non-blank and the chain column "
                        "is empty.")
    p.add_argument("--no-charmm", dest="charmm", action="store_false",
                   help="Disable CHARMM handling even if auto-detected.")
    p.add_argument("--split-on-gaps", action="store_true",
                   help="Start a new chain at residue-number resets or gaps.")
    p.add_argument("--drop-solvent", action="store_true",
                   help="Remove water and common crystallisation additives "
                        "(HOH, EDO, GOL, SO4, CIT, EPE, TLA, MES, TRS, PEG, ACT, DMS).")
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


def required_int(s: str, field: str, line_number: int) -> int:
    try:
        return int(s.strip())
    except ValueError as exc:
        raise SystemExit(
            f"Invalid {field} at PDB line {line_number}: {s!r}"
        ) from exc


def required_float(s: str, field: str, line_number: int) -> float:
    try:
        return float(s.strip())
    except ValueError as exc:
        raise SystemExit(
            f"Invalid {field} at PDB line {line_number}: {s!r}"
        ) from exc


def optional_float(s: str, default: float, field: str, line_number: int) -> float:
    if not s.strip():
        return default
    return required_float(s, field, line_number)


def infer_element(
    atom_name: str,
    element_field: str = "",
    resname: str = "",
) -> str:
    e = element_field.strip().upper()
    if e:
        return e[:2]
    name = atom_name.strip().replace("'", "").replace("*", "")
    if not name:
        return ""
    # For names such as 1H5, H5', C1', OP1, CL, NA.
    if name[0].isdigit() and len(name) > 1:
        name = name[1:]
    if resname.strip().upper() in ION_MAP:
        ion = ION_MAP[resname.strip().upper()]
        return "CA" if ion == "CAL" else ion[:2]
    # Common two-letter elements outside protein atom naming.
    up = name.upper()
    if up.startswith("CL"):
        return "CL"
    if up.startswith("NA"):
        return "NA"
    if up.startswith("MG"):
        return "MG"
    if up.startswith("ZN"):
        return "ZN"
    return up[0]


SOLVENT_AND_ADDITIVES = {
    "HOH", "WAT", "TIP3", "SOL", "DOD",
    "EDO", "GOL", "PEG", "PG4", "SO4", "PO4", "CIT", "EPE", "TLA", "MES", "TRS",
    "ACT", "DMS", "IMD", "FMT", "NO3", "MPD", "BME",
}


def detect_charmm(path: str) -> bool:
    """True when the file carries a four-character residue name and no chain id."""
    wide = False
    chain_used = False
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            if len(line) > 20 and line[20:21].strip():
                wide = True
            if len(line) > 21 and line[21:22].strip():
                chain_used = True
    return wide and not chain_used


def parse_pdb(
    path: str,
    model: int = 1,
    drop_solvent: bool = False,
    charmm: bool = False,
) -> Tuple[List[AtomRecord], List[str]]:
    """Read one model and select one coordinate for each alternate location."""
    provisional: List[AtomRecord] = []
    notes: List[str] = []
    ter_pending = False
    current_model: Optional[int] = None
    model_records_seen = False
    target_model_seen = False
    model_sequence = 0
    n_model_skipped = 0
    n_solvent_skipped = 0

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line_number, line in enumerate(fh, start=1):
            rec = line[:6].strip()
            if rec == "MODEL":
                model_records_seen = True
                model_sequence += 1
                fields = line[6:].split()
                current_model = int(fields[0]) if fields and fields[0].isdigit() else model_sequence
                if current_model == model:
                    target_model_seen = True
                ter_pending = False
                continue
            if rec == "ENDMDL":
                current_model = None
                ter_pending = False
                continue
            if rec == "TER":
                if not model_records_seen or current_model == model:
                    ter_pending = True
                continue
            if rec not in {"ATOM", "HETATM"}:
                continue
            if model_records_seen and current_model != model:
                n_model_skipped += 1
                continue

            padded = line.rstrip("\n").ljust(80)
            resn = (padded[17:21] if charmm else padded[17:20]).strip().upper()
            if drop_solvent and resn in SOLVENT_AND_ADDITIVES:
                n_solvent_skipped += 1
                continue

            atom_name = padded[12:16].strip()
            atom = AtomRecord(
                record=rec,
                input_atom_index=0,
                input_serial=required_int(padded[6:11], "atom serial", line_number),
                atom_name=atom_name,
                altloc=padded[16:17].strip(),
                resname_orig=resn,
                chain_orig=padded[21:22].strip(),
                resseq_orig=required_int(padded[22:26], "residue number", line_number),
                icode_orig=padded[26:27].strip(),
                x=required_float(padded[30:38], "x coordinate", line_number),
                y=required_float(padded[38:46], "y coordinate", line_number),
                z=required_float(padded[46:54], "z coordinate", line_number),
                occ=optional_float(padded[54:60], 1.0, "occupancy", line_number),
                bfac=optional_float(padded[60:66], 0.0, "B factor", line_number),
                segid=padded[72:76].strip(),
                element=infer_element(atom_name, padded[76:78], resn),
                charge=padded[78:80].strip(),
                line_number=line_number,
                ter_before=ter_pending,
            )
            provisional.append(atom)
            ter_pending = False

    if model_records_seen and not target_model_seen:
        raise SystemExit(f"MODEL {model} was not found in {path}")
    if n_model_skipped:
        notes.append(
            f"MODEL filter: kept MODEL {model}; dropped {n_model_skipped} atom records."
        )
    if n_solvent_skipped:
        notes.append(f"Solvent filter: dropped {n_solvent_skipped} atom records.")

    groups = {}
    for order, atom in enumerate(provisional):
        key = (
            atom.chain_orig,
            atom.segid,
            atom.resseq_orig,
            atom.icode_orig,
            atom.resname_orig,
            atom.atom_name,
        )
        groups.setdefault(key, []).append((order, atom))

    chosen_orders = set()
    replacements = {}
    for entries in groups.values():
        anchor_order = entries[0][0]
        blank = [entry for entry in entries if not entry[1].altloc]
        candidates = blank or entries
        chosen = min(
            candidates,
            key=lambda entry: (
                -entry[1].occ,
                0 if entry[1].altloc == "A" else 1,
                entry[0],
            ),
        )
        chosen_orders.add(anchor_order)
        replacements[anchor_order] = replace(
            chosen[1],
            altloc="",
            ter_before=any(atom.ter_before for _, atom in entries),
        )

    atoms: List[AtomRecord] = []
    for order, atom in enumerate(provisional):
        if order not in chosen_orders:
            continue
        atoms.append(
            replace(replacements[order], input_atom_index=len(atoms) + 1)
        )

    n_altloc_skipped = len(provisional) - len(atoms)
    if n_altloc_skipped:
        notes.append(
            f"altLoc filter: selected one coordinate per atom; dropped {n_altloc_skipped} records."
        )
    if notes:
        notes.append("The -a range is applied after MODEL, altLoc and solvent filtering.")
    return atoms, notes

def normalize_atom_name(name: str) -> str:
    # SAXS.cpp recognizes apostrophe names.
    name = name.strip().replace('*', "'")
    # Common aliases from some builders; keep O1P/O2P too because SAXS.cpp accepts them.
    aliases = {
        "O1P": "O1P", "O2P": "O2P", "O3P": "O3P",
        "OP1": "OP1", "OP2": "OP2", "OP3": "OP3",
        "H5'1": "H5'1", "H5'2": "H5'2", "H2'1": "H2'1", "H2'2": "H2'2",
    }
    return aliases.get(name, name)


def is_rna_resname_out(resname_out: str) -> bool:
    r = resname_out.strip().upper()
    return r in {"A", "C", "G", "U", "A3", "C3", "G3", "U3", "A5", "C5", "G5", "U5", "AT", "CT", "GT", "UT"}


def normalize_atom_name_for_residue(name: str, resname_out: str) -> str:
    """Normalize residue-dependent atom names."""
    n = normalize_atom_name(name)
    if is_rna_resname_out(resname_out):
        if n == "H2''":
            return "H2'"
        if n == "H2'":
            return "HO2'"
    return n


def base_resname(resname: str) -> str:
    r = resname.strip().upper()
    r = r.replace("5'", "").replace("3'", "")
    if re.fullmatch(r"[ACGU][35T]", r):
        return r[0]
    if re.fullmatch(r"D[ACGT][35T]", r):
        return r[:2]
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


def existing_terminal_suffix(resname: str) -> str:
    r = resname.strip().upper()
    if re.fullmatch(r"[ACGU][35T]", r):
        return r[-1]
    if re.fullmatch(r"D[ACGT][35T]", r):
        return r[-1]
    return ""


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


def split_into_chains(
    residues: Sequence[Sequence[AtomRecord]],
    split_on_gaps: bool = False,
) -> List[List[List[AtomRecord]]]:
    chains: List[List[List[AtomRecord]]] = []
    current: List[List[AtomRecord]] = []
    prev_first: Optional[AtomRecord] = None

    for residue in residues:
        first = residue[0]
        source_id = (first.chain_orig, first.segid)
        previous_id = (
            (prev_first.chain_orig, prev_first.segid)
            if prev_first is not None else source_id
        )
        new_chain = not current
        if current and prev_first is not None:
            if first.ter_before or source_id != previous_id:
                new_chain = True
            elif split_on_gaps and first.resseq_orig != prev_first.resseq_orig + 1:
                new_chain = True

        if new_chain:
            if current:
                chains.append(current)
            current = [list(residue)]
        else:
            current.append(list(residue))
        prev_first = first

    if current:
        chains.append(current)
    return chains


def convert_atoms(
    atoms: List[AtomRecord],
    split_on_gaps: bool = False,
) -> Tuple[List[AtomRecord], List[str]]:
    if not atoms:
        raise SystemExit("The atom selection is empty.")

    log: List[str] = []
    errors: List[str] = []
    residues = group_residues(atoms)
    chains = split_into_chains(residues, split_on_gaps=split_on_gaps)
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

        atom_count = sum(len(residue) for residue in chain_residues)
        log.append(
            f"Chain {chain_idx}: ID {chain_id}; residues={len(chain_residues)}; "
            f"atoms={atom_count}."
        )
        if len(chain_residues) > 9999:
            errors.append(f"Chain {chain_id} has more than 9999 residues.")

        for residue_index, res_atoms in enumerate(chain_residues, start=1):
            orig = res_atoms[0]
            base = base_resname(orig.resname_orig)
            names = residue_atom_names(res_atoms)
            is_first = residue_index == 1
            is_last = residue_index == len(chain_residues)

            if base in UNSUPPORTED_SITE_NAMES:
                errors.append(
                    f"Unsupported glycosylation-site residue {base} at "
                    f"chain {chain_id}, input residue {orig.resseq_orig}."
                )
                continue

            if base in GLYCOSYLATION_SITE_MAP:
                parent = GLYCOSYLATION_SITE_MAP[base]
                log.append(
                    f"Glycosylation site: {base}{orig.resseq_orig} -> {parent}."
                )
                base = parent

            glycan = canonical_glycan(base, names)
            ion = ION_MAP.get(base)
            suffix = existing_terminal_suffix(orig.resname_orig) or terminal_suffix(
                base, names, is_first, is_last
            )

            if glycan is not None:
                resname_out = glycan
                if glycan != base:
                    log.append(
                        f"Glycan: {orig.resname_orig}{orig.resseq_orig} -> {glycan}."
                    )
            elif base in GLYCAN_UNSUPPORTED:
                errors.append(
                    f"Monosaccharide {base} at chain {chain_id}, input residue "
                    f"{orig.resseq_orig} has no ONEBEAD parameters."
                )
                continue
            elif ion is not None:
                if len(res_atoms) != 1:
                    errors.append(
                        f"Ion residue {base}{orig.resseq_orig} contains "
                        f"{len(res_atoms)} atoms; ONEBEAD expects one atom."
                    )
                    continue
                resname_out = ion
            elif is_rna_base(base) or is_dna_base(base):
                resname_out = base + suffix
                unknown = sorted(name for name in names if name not in KNOWN_ONEBEAD_NUC_ATOMS)
                if unknown:
                    errors.append(
                        f"Unknown nucleic-acid atom names in {orig.resname_orig}"
                        f"{orig.resseq_orig}: {','.join(unknown)}."
                    )
                    continue
                if suffix:
                    log.append(
                        f"Terminal residue: {orig.resname_orig}{orig.resseq_orig} "
                        f"-> {resname_out}."
                    )
            elif base in PROTEIN_NAMES:
                resname_out = base
            else:
                errors.append(
                    f"Unsupported residue {orig.resname_orig} at chain {chain_id}, "
                    f"input residue {orig.resseq_orig}."
                )
                continue

            output_names = []
            for atom in res_atoms:
                if glycan is not None:
                    atom_out = normalize_glycan_atom_name(atom.atom_name, resname_out)
                else:
                    atom_out = normalize_atom_name_for_residue(atom.atom_name, resname_out)
                output_names.append(atom_out)

            duplicates = sorted(
                name for name, count in Counter(output_names).items() if count > 1
            )
            if duplicates:
                errors.append(
                    f"Duplicate output atom names in {resname_out}{residue_index} "
                    f"chain {chain_id}: {','.join(duplicates)}."
                )
                continue

            for atom, atom_out in zip(res_atoms, output_names):
                if len(atom_out) > 4:
                    errors.append(
                        f"Atom name {atom_out!r} exceeds four PDB columns at input "
                        f"line {atom.line_number}."
                    )
                    continue
                if atom_out != normalize_atom_name(atom.atom_name):
                    log.append(
                        f"Atom name: {normalize_atom_name(atom.atom_name)} -> {atom_out} "
                        f"in {resname_out}{residue_index}."
                    )
                converted.append(
                    replace(
                        atom,
                        chain_out=chain_id,
                        resseq_out=residue_index,
                        resname_out=resname_out,
                        atom_name_out=atom_out,
                    )
                )

    if errors:
        details = "\n".join(f"- {error}" for error in errors)
        raise SystemExit(f"Conversion stopped:\n{details}")
    if len(converted) > 99999:
        raise SystemExit("The output contains more than 99999 atoms and cannot fit PDB columns.")
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
    element = (
        a.element or infer_element(a.atom_name_out, resname=a.resname_out)
    ).upper()[:2]
    # Keep segid in columns 73-76 when available; otherwise use chain-friendly placeholder.
    segid = (a.segid or ("RNA" + chain if is_nucleic_base_name(a.resname_out) else chain))[:4]
    line = (
        f"{rec:<6}{out_serial:5d} {atom_field}{a.altloc[:1]:1s}"
        f"{res_field:>3s} {chain:1s}{a.resseq_out:4d}{a.icode_orig[:1]:1s}   "
        f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occ:6.2f}{a.bfac:6.2f}"
        f"      {segid:<4s}{element:>2s}{a.charge:>2s}"
    )
    if len(line) != 80:
        raise SystemExit(
            f"PDB field overflow for output atom {out_serial} "
            f"({a.resname_out}{a.resseq_out}/{a.atom_name_out})."
        )
    return line


def format_ter(serial: int, last_atom: AtomRecord) -> str:
    return f"TER   {serial:5d}      {last_atom.resname_out:>3s} {last_atom.chain_out[:1]:1s}{last_atom.resseq_out:4d}"


def write_pdb(atoms: List[AtomRecord], path: str) -> None:
    if not atoms:
        raise SystemExit("Cannot write an empty PDB template.")
    for atom in atoms:
        if len(atom.resname_out) > 3:
            raise SystemExit(f"Residue name {atom.resname_out!r} exceeds three PDB columns.")
        if len(atom.charge) > 2:
            raise SystemExit(f"Charge field {atom.charge!r} exceeds two PDB columns.")
    with open(path, 'w', encoding='utf-8') as out:
        out.write("REMARK Prepared by pdb2plmd\n")
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
        log.write("pdb2plmd verbose log\n")
        log.write("=======================\n")
        log.write(f"Input:  {args.input}\n")
        log.write(f"Output: {args.output}\n")
        log.write(f"Atom range expression: {args.atoms}\n")
        log.write(f"Total ATOM/HETATM records in input: {len(all_atoms)}\n")
        log.write(f"Selected atoms: {len(selected_indices)}\n")
        if selected_indices:
            log.write(f"Selected input atom-order range: {selected_indices[0]}..{selected_indices[-1]}\n")
        log.write("\n")
        log.write("Notes:\n")
        log.write("- The -a range is applied after MODEL, altLoc and solvent filtering.\n")
        log.write("- It refers to ATOM/HETATM order, not PDB serial.\n")
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

    use_charmm = detect_charmm(args.input) if args.charmm is None else args.charmm
    if use_charmm:
        print("CHARMM input detected: reading the four-character residue name from columns 18-21.")
    all_atoms, parse_notes = parse_pdb(
        args.input,
        model=args.model,
        charmm=use_charmm,
        drop_solvent=args.drop_solvent,
    )
    for n in parse_notes:
        print(n)
    if not all_atoms:
        raise SystemExit(f"No ATOM/HETATM records found in {args.input}")
    selected_indices = parse_range(args.atoms, len(all_atoms))
    index_set = set(selected_indices)
    selected_atoms = [a for a in all_atoms if a.input_atom_index in index_set]
    if not selected_atoms:
        raise SystemExit("The atom selection is empty.")

    converted, log_lines = convert_atoms(
        selected_atoms,
        split_on_gaps=args.split_on_gaps,
    )
    log_lines = parse_notes + log_lines
    if any(
        not math.isfinite(value)
        for atom in converted
        for value in (atom.x, atom.y, atom.z, atom.occ, atom.bfac)
    ):
        raise SystemExit("The selected atoms contain non-finite numeric values.")

    n_h = sum(1 for a in converted if a.element == "H")
    if converted and n_h < 0.2 * len(converted):
        msg = (
            f"WARNING only {n_h}/{len(converted)} atoms are hydrogens "
            f"({100.0 * n_h / len(converted):.1f}%). The ONEBEAD template is "
            f"expected to be all-atom."
        )
        print(msg)
        log_lines.append(msg)
    write_pdb(converted, args.output)
    write_log(args.log, args, all_atoms, selected_indices, converted, log_lines)

    print(f"Wrote {args.output} with {len(converted)} atoms.")
    if args.log:
        print(f"Wrote verbose log: {args.log}")
    print("Atom order was preserved from the selected input ATOM/HETATM order.")


if __name__ == "__main__":
    main()
