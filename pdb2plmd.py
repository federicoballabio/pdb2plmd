#!/usr/bin/env python3
"""Prepare an ordered PDB template for PLUMED SAXS.cpp ONEBEAD."""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from collections import Counter
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Sequence, Tuple

VERSION = "v26080804"

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
ION_ELEMENT = {
    "NA": "NA", "K": "K", "CL": "CL", "CAL": "CA", "MG": "MG",
    "ZN": "ZN", "FE2": "FE", "FE3": "FE", "MN": "MN",
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

# GLYCAM terminal caps are separate moieties, not monosaccharide beads.
TERMINAL_CAP_NAMES = {"ROH", "OME", "TBT"}

MAX_REPORTED_ISSUE_GROUPS = 20

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


@dataclass(frozen=True)
class CompatibilityIssue:
    key: str
    message: str
    hint: str = ""


def stop_for_compatibility_issues(issues: Sequence[CompatibilityIssue]) -> None:
    """Stop once, grouping repeated incompatibilities into a concise report."""
    if not issues:
        return

    grouped = {}
    for issue in issues:
        if issue.key not in grouped:
            grouped[issue.key] = [issue, 0]
        grouped[issue.key][1] += 1

    lines = [
        f"ERROR: found {len(issues)} incompatible residue/moiety issue(s) "
        f"in {len(grouped)} group(s)."
    ]
    shown = list(grouped.values())[:MAX_REPORTED_ISSUE_GROUPS]
    for issue, count in shown:
        label = issue.key.split(":", 1)[0]
        repeat = f" ({count} occurrences; first shown)" if count > 1 else ""
        lines.append(f"- [{label}] {issue.message}{repeat}")

    if len(grouped) > len(shown):
        lines.append(
            f"- ... and {len(grouped) - len(shown)} more issue group(s)."
        )

    hints = []
    for issue in issues:
        if issue.hint and issue.hint not in hints:
            hints.append(issue.hint)
    if hints:
        lines.append("")
        lines.extend(f"HINT: {hint}" for hint in hints)

    lines.append("Conversion stopped before writing the output PDB.")
    raise SystemExit("\n".join(lines))


def partial_selection_issues(
    all_atoms: Sequence[AtomRecord], selected_indices: Sequence[int]
) -> List[str]:
    """Return fatal diagnostics when -a keeps only part of an input residue."""
    selected = set(selected_indices)
    issues: List[str] = []
    for residue in group_residues(all_atoms):
        kept = sum(atom.input_atom_index in selected for atom in residue)
        if kept in {0, len(residue)}:
            continue
        first = residue[0]
        chain = first.chain_orig or derive_chain_from_segid(first.segid) or "?"
        omitted = len(residue) - kept
        issues.append(
            f"{first.resname_orig} at chain {chain}, input residue "
            f"{first.resseq_orig}: selected {kept}/{len(residue)} atoms "
            f"and left {omitted} atom(s) outside the requested range."
        )
    if len(issues) > MAX_REPORTED_ISSUE_GROUPS:
        omitted = len(issues) - MAX_REPORTED_ISSUE_GROUPS
        issues = issues[:MAX_REPORTED_ISSUE_GROUPS]
        issues.append(f"... and {omitted} more incomplete selected residues.")
    return issues


def selection_boundary_notes(
    all_atoms: Sequence[AtomRecord], selected_indices: Sequence[int]
) -> Tuple[List[str], bool]:
    """Describe the selected boundary and flag PDB-serial/atom-order mismatches."""
    if not selected_indices:
        return [], False

    def describe(idx: int) -> str:
        atom = all_atoms[idx - 1]
        chain = atom.chain_orig or derive_chain_from_segid(atom.segid) or "?"
        return (
            f"atom-order {idx} -> PDB serial {atom.input_serial}: "
            f"{atom.atom_name} {atom.resname_orig} chain {chain} "
            f"residue {atom.resseq_orig}"
        )

    first_idx = selected_indices[0]
    last_idx = selected_indices[-1]
    notes = [
        f"Selection start: {describe(first_idx)}.",
        f"Selection end: {describe(last_idx)}.",
    ]
    mismatch = (
        all_atoms[first_idx - 1].input_serial != first_idx
        or all_atoms[last_idx - 1].input_serial != last_idx
    )
    if mismatch:
        notes.append(
            "Selection numbering note: -a/--atoms uses 1-based ATOM/HETATM "
            "record order after MODEL/altLoc/solvent filtering, not the PDB "
            "serial field. TER and other non-coordinate records are not counted."
        )
    return notes, mismatch


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare an ordered PDB template for PLUMED SAXS.cpp ONEBEAD."
    )
    p.add_argument("--version", action="version", version=f"pdb2plmd {VERSION}")
    p.add_argument("-i", "--input", required=True, help="Input PDB extracted from the simulation/TPR selection.")
    p.add_argument("-o", "--output", required=True, help="Output SAXS.cpp-compatible template PDB.")
    selection = p.add_mutually_exclusive_group()
    selection.add_argument("-a", "--atoms", default=None,
                           help="1-based ATOM/HETATM record-order range to keep, e.g. "
                                "'1-1062' or '1-100,150,200-250'. TER records are not "
                                "counted. Default: all when neither -a nor -s is given.")
    selection.add_argument("-s", "--serials", default=None,
                           help="PDB atom-serial range to keep, e.g. '1-1069'. Selects "
                                "ATOM/HETATM records whose PDB serial field falls in the "
                                "requested range; TER serials are naturally skipped. Use -a "
                                "for files with repeated/wrapped atom serials.")
    p.add_argument("--model", type=int, default=1,
                   help="MODEL serial to keep from a multi-model file. Default: 1.")
    # Legacy expert overrides retained for backward compatibility. Normal users do not
    # need to specify an input convention: PDB layout is inspected automatically.
    p.add_argument("--charmm", dest="charmm", action="store_true", default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("--no-charmm", dest="charmm", action="store_false",
                   help=argparse.SUPPRESS)
    p.add_argument("--split-on-gaps", action="store_true",
                   help="Start a new chain at residue-number resets or gaps.")
    p.add_argument("-altloc", "--altloc", default="auto", metavar="LABEL",
                   help="Alternate-location handling. Default: auto selects one coherent "
                        "residue-level conformer by completeness then occupancy. Give a "
                        "single label such as A or B to force that label for every residue "
                        "with alternate coordinates.")
    p.add_argument("--drop-solvent", action="store_true",
                   help="Remove water and common crystallisation additives "
                        "(including HOH, EDO, EOH, GOL, SO4, CIT, EPE, TLA, MES, TRS, PEG, ACT, DMS).")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose mode. On success write <output-stem>.log in addition to the PDB. "
                        "On failure a log is written in both standard and verbose modes.")
    # Backward-compatible legacy switch. It is intentionally hidden from --help so
    # the public interface remains standard mode versus -v/--verbose.
    p.add_argument("-g", "--log", nargs="?", const="__AUTO__", default=None,
                   help=argparse.SUPPRESS)
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


def parse_serial_selection(expr: str, all_atoms: Sequence[AtomRecord]) -> List[int]:
    """Return input atom-order indices selected by the PDB serial field.

    Missing integers inside a requested serial interval are allowed because TER and other
    non-coordinate PDB records can consume serial numbers. Repeated coordinate-record serials
    are rejected because a serial-based selection would otherwise be ambiguous.
    """
    expression = (expr or "").strip().lower()
    if expression in {"", "all", "*"}:
        return [atom.input_atom_index for atom in all_atoms]

    intervals: List[Tuple[int, int]] = []
    for part in expression.split(','):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise SystemExit(f"Invalid PDB serial range component: {part!r}")
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        if start < 1 or end < 1 or end < start:
            raise SystemExit(f"Invalid PDB serial range component: {part!r}")
        intervals.append((start, end))
    if not intervals:
        return []

    def requested(serial: int) -> bool:
        return any(start <= serial <= end for start, end in intervals)

    selected_atoms = [atom for atom in all_atoms if requested(atom.input_serial)]
    if not selected_atoms:
        raise SystemExit(
            f"PDB serial selection {expr!r} did not match any ATOM/HETATM records."
        )

    by_serial: Dict[int, List[int]] = {}
    for atom in selected_atoms:
        by_serial.setdefault(atom.input_serial, []).append(atom.input_atom_index)
    duplicated = {serial: idxs for serial, idxs in by_serial.items() if len(idxs) > 1}
    if duplicated:
        preview = sorted(duplicated.items())[:5]
        detail = "; ".join(
            f"serial {serial} occurs at atom-order {','.join(map(str, idxs))}"
            for serial, idxs in preview
        )
        if len(duplicated) > len(preview):
            detail += f"; ... and {len(duplicated) - len(preview)} more repeated serial(s)"
        raise SystemExit(
            "PDB serial selection is ambiguous because the selected ATOM/HETATM records "
            f"contain repeated atom serials ({detail}). Use -a/--atoms record-order "
            "selection for files with repeated or wrapped PDB serials."
        )

    return [atom.input_atom_index for atom in selected_atoms]


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
    "ACT", "DMS", "IMD", "FMT", "NO3", "MPD", "BME", "EOH",
}


@dataclass(frozen=True)
class PdbFormatPreflight:
    sampled_atoms: int
    chain_records: int
    segid_records: int
    four_char_residue_records: int
    style: str


def preflight_pdb_format(path: str, sample_limit: int = 5000) -> PdbFormatPreflight:
    """Inspect enough ATOM/HETATM records to describe the PDB layout.

    Parsing itself remains record-aware, so this classification is informative rather
    than a brittle global switch. Limiting the scan keeps the preflight cheap for very
    large solvated systems.
    """
    sampled = 0
    chain_records = 0
    segid_records = 0
    four_char_records = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            padded = line.rstrip("\n").ljust(80)
            sampled += 1
            chain = padded[21:22].strip()
            segid = padded[72:76].strip()
            if chain:
                chain_records += 1
            if segid:
                segid_records += 1
            # In the CHARMM four-character extension, column 21 belongs to the
            # residue name while the standard PDB chain column (22) is blank.
            if padded[20:21].strip() and not chain:
                four_char_records += 1
            if sampled >= sample_limit:
                break

    if sampled == 0:
        style = "no ATOM/HETATM records in preflight sample"
    elif four_char_records:
        style = "CHARMM/extended PDB style (four-character residue names observed)"
    elif segid_records >= max(1, int(0.8 * sampled)) and chain_records == 0:
        style = "CHARMM/CHARMM-GUI style (SEGID populated; PDB chain column blank)"
    elif chain_records:
        style = "standard PDB chain-column style"
    elif segid_records:
        style = "SEGID-based PDB style"
    else:
        style = "chainless standard PDB style"

    return PdbFormatPreflight(
        sampled_atoms=sampled,
        chain_records=chain_records,
        segid_records=segid_records,
        four_char_residue_records=four_char_records,
        style=style,
    )


def _read_resname(padded: str, charmm_override: Optional[bool]) -> str:
    """Read a 3- or 4-character residue name without requiring a user mode flag."""
    if charmm_override is True:
        field = padded[17:21]
    elif charmm_override is False:
        field = padded[17:20]
    elif padded[20:21].strip() and not padded[21:22].strip():
        field = padded[17:21]
    else:
        field = padded[17:20]
    return field.strip().upper()


def normalize_altloc_request(value: str) -> str:
    """Return 'auto' or one explicit single-character altLoc label."""
    value = (value or "auto").strip()
    if value.lower() == "auto":
        return "auto"
    if len(value) != 1 or value.isspace():
        raise SystemExit(
            "--altloc expects 'auto' or one nonblank PDB altLoc character, "
            "for example --altloc A or --altloc B."
        )
    return value.upper() if value.isalpha() else value


def _choose_duplicate_record(entries: Sequence[Tuple[int, AtomRecord]]) -> Tuple[int, AtomRecord]:
    """Choose deterministically among duplicate records for one atom/altLoc label."""
    return min(entries, key=lambda entry: (-entry[1].occ, entry[0]))


def _select_altloc_for_residue(
    entries: Sequence[Tuple[int, AtomRecord]],
    requested: str,
) -> Tuple[List[Tuple[int, AtomRecord]], Optional[str], List[str]]:
    """Select a coherent residue-level alternate conformation.

    Blank-altLoc records are treated as common atoms. For a nonblank conformer,
    records from other nonblank labels are never mixed into the selected residue.
    """
    by_name = {}
    labels = set()
    for order, atom in entries:
        by_name.setdefault(atom.atom_name, []).append((order, atom))
        if atom.altloc:
            labels.add(atom.altloc)

    if not labels:
        chosen = []
        for atom_entries in by_name.values():
            blank = [entry for entry in atom_entries if not entry[1].altloc]
            chosen.append(_choose_duplicate_record(blank or atom_entries))
        return sorted(chosen), None, []

    labels_sorted = sorted(labels)
    if requested != "auto":
        if requested not in labels:
            first = entries[0][1]
            available = ",".join(labels_sorted)
            raise SystemExit(
                f"Requested altLoc {requested} is unavailable for "
                f"{first.resname_orig} at chain {first.chain_orig or '-'}, input "
                f"residue {first.resseq_orig}. Available alternate locations: {available}."
            )
        selected_label = requested
    else:
        scored = []
        for label in labels_sorted:
            selected_for_label = []
            for atom_entries in by_name.values():
                blank = [entry for entry in atom_entries if not entry[1].altloc]
                labelled = [entry for entry in atom_entries if entry[1].altloc == label]
                candidates = blank or labelled
                if candidates:
                    selected_for_label.append(_choose_duplicate_record(candidates))
            completeness = len(selected_for_label)
            occupancy_sum = sum(atom.occ for _, atom in selected_for_label)
            scored.append((
                -completeness,
                -occupancy_sum,
                0 if label == "A" else 1,
                label,
            ))
        selected_label = min(scored)[3]

    chosen = []
    for atom_entries in by_name.values():
        blank = [entry for entry in atom_entries if not entry[1].altloc]
        labelled = [entry for entry in atom_entries if entry[1].altloc == selected_label]
        candidates = blank or labelled
        if candidates:
            chosen.append(_choose_duplicate_record(candidates))

    return sorted(chosen), selected_label, labels_sorted


def parse_pdb(
    path: str,
    model: int = 1,
    drop_solvent: bool = False,
    charmm: Optional[bool] = None,
    altloc: str = "auto",
) -> Tuple[List[AtomRecord], List[str]]:
    """Read one model and resolve alternate locations at residue level."""
    provisional: List[AtomRecord] = []
    notes: List[str] = []
    ter_pending = False
    current_model: Optional[int] = None
    model_records_seen = False
    target_model_seen = False
    model_sequence = 0
    n_model_skipped = 0
    n_solvent_skipped = 0
    altloc_request = normalize_altloc_request(altloc)

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
            resn = _read_resname(padded, charmm)
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

    residue_groups = {}
    for order, atom in enumerate(provisional):
        key = (
            atom.chain_orig,
            atom.segid,
            atom.resseq_orig,
            atom.icode_orig,
            atom.resname_orig,
        )
        residue_groups.setdefault(key, []).append((order, atom))

    selected_entries: List[Tuple[int, AtomRecord]] = []
    altloc_selections = []
    for entries in residue_groups.values():
        chosen, label, available = _select_altloc_for_residue(entries, altloc_request)
        ter_before = any(atom.ter_before for _, atom in entries)
        if chosen:
            first_order = min(order for order, _ in chosen)
            for order, atom in chosen:
                selected_entries.append((
                    order,
                    replace(
                        atom,
                        altloc="",
                        ter_before=(ter_before and order == first_order),
                    ),
                ))
        if label is not None:
            first = entries[0][1]
            altloc_selections.append(
                f"{first.resname_orig} {first.chain_orig or '-'}{first.resseq_orig}"
                f"{first.icode_orig or ''}: {label} from {','.join(available)}"
            )

    atoms: List[AtomRecord] = []
    for _, atom in sorted(selected_entries, key=lambda entry: entry[0]):
        atoms.append(replace(atom, input_atom_index=len(atoms) + 1))

    n_altloc_skipped = len(provisional) - len(atoms)
    if altloc_selections:
        mode = "automatic" if altloc_request == "auto" else f"forced {altloc_request}"
        notes.append(
            f"altLoc filter: residue-level {mode} selection; "
            f"resolved {len(altloc_selections)} residue(s) and dropped "
            f"{n_altloc_skipped} alternate record(s)."
        )
        preview = altloc_selections[:10]
        notes.extend(f"altLoc selection: {item}." for item in preview)
        if len(altloc_selections) > len(preview):
            notes.append(
                f"altLoc selection: ... and {len(altloc_selections) - len(preview)} more residue(s)."
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


def normalize_atom_names_for_residue(
    res_atoms: Sequence[AtomRecord],
    resname_out: str,
) -> List[str]:
    """Normalize atom names using residue-level context where required."""
    normalized = [normalize_atom_name(atom.atom_name) for atom in res_atoms]
    if not is_rna_resname_out(resname_out):
        return normalized

    source_names = set(normalized)
    legacy_rna_h2_pair = (
        "H2''" in source_names
        and "H2'" in source_names
        and "HO2'" not in source_names
    )
    if not legacy_rna_h2_pair:
        return normalized

    return [
        "H2'" if name == "H2''" else "HO2'" if name == "H2'" else name
        for name in normalized
    ]


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
    errors: List[CompatibilityIssue] = []
    residues = group_residues(atoms)
    chains = split_into_chains(residues, split_on_gaps=split_on_gaps)
    used_chains = set()
    converted: List[AtomRecord] = []
    log.append(f"Input selected atoms: {len(atoms)}")
    log.append(f"Input selected residues: {len(residues)}")
    log.append(f"Inferred output chains: {len(chains)}")
    if len(chains) > len(CHAIN_IDS):
        cause = " with --split-on-gaps" if split_on_gaps else ""
        raise SystemExit(
            f"Inferred {len(chains)} output chains{cause}, but the PDB output format "
            f"supports at most {len(CHAIN_IDS)} single-character chain IDs. "
            "Disable --split-on-gaps if those breaks are not intended, or split the "
            "system into separate templates."
        )

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
            errors.append(
                CompatibilityIssue(
                    "pdb-limit:residues",
                    f"Chain {chain_id} has more than 9999 residues and cannot fit "
                    "the PDB residue-number field.",
                )
            )

        for residue_index, res_atoms in enumerate(chain_residues, start=1):
            orig = res_atoms[0]
            base = base_resname(orig.resname_orig)
            names = residue_atom_names(res_atoms)
            is_first = residue_index == 1
            is_last = residue_index == len(chain_residues)

            if base in UNSUPPORTED_SITE_NAMES:
                errors.append(
                    CompatibilityIssue(
                        f"unsupported-site:{base}",
                        f"Glycosylation-site residue {base} at chain {chain_id}, "
                        f"input residue {orig.resseq_orig} is not compatible with "
                        "ONEBEAD. Hydroxyproline contains OD1, which has no PRO "
                        "LCPO mapping.",
                        "Do not approximate OLP as PRO; a dedicated validated "
                        "ONEBEAD/LCPO mapping is required.",
                    )
                )
                continue

            if base in TERMINAL_CAP_NAMES:
                errors.append(
                    CompatibilityIssue(
                        f"terminal-cap:{base}",
                        f"Terminal capping group {base} at chain {chain_id}, input "
                        f"residue {orig.resseq_orig} is not a monosaccharide and has "
                        "no ONEBEAD bead.",
                        "Remove GLYCAM terminal caps before conversion.",
                    )
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
                unexpected_elements = sorted(
                    {atom.element.upper() for atom in res_atoms}
                    - {"H", "C", "N", "O"}
                )
                if unexpected_elements:
                    errors.append(
                        CompatibilityIssue(
                            f"modified-glycan:{glycan}:{','.join(unexpected_elements)}",
                            f"Supported glycan {orig.resname_orig} at chain {chain_id}, "
                            f"input residue {orig.resseq_orig} contains unexpected "
                            f"element(s) {','.join(unexpected_elements)}. This indicates "
                            "a substituted moiety not covered by the implemented "
                            "monosaccharide bead.",
                            "Remove the unsupported substituent or implement and validate "
                            "a dedicated ONEBEAD parameterization.",
                        )
                    )
                    continue
                resname_out = glycan
                if glycan != base:
                    log.append(
                        f"Glycan: {orig.resname_orig}{orig.resseq_orig} -> {glycan}."
                    )
            elif base in GLYCAN_UNSUPPORTED:
                errors.append(
                    CompatibilityIssue(
                        f"unsupported-glycan:{base}",
                        f"Monosaccharide {base} at chain {chain_id}, input residue "
                        f"{orig.resseq_orig} has no ONEBEAD parameters.",
                        "No chemically similar monosaccharide is substituted "
                        "automatically.",
                    )
                )
                continue
            elif ion is not None:
                if len(res_atoms) != 1:
                    errors.append(
                        CompatibilityIssue(
                            f"ion-atom-count:{base}",
                            f"Ion residue {base} at chain {chain_id}, input residue "
                            f"{orig.resseq_orig} contains {len(res_atoms)} atoms; "
                            "ONEBEAD expects exactly one atom.",
                        )
                    )
                    continue
                resname_out = ion
                if ion != base:
                    log.append(
                        f"Ion: {orig.resname_orig}{orig.resseq_orig} -> {ion}."
                    )
            elif is_rna_base(base) or is_dna_base(base):
                resname_out = base + suffix
                unknown = sorted(
                    name for name in names if name not in KNOWN_ONEBEAD_NUC_ATOMS
                )
                if unknown:
                    errors.append(
                        CompatibilityIssue(
                            f"nucleic-atom-names:{orig.resname_orig}:{','.join(unknown)}",
                            f"Nucleic-acid residue {orig.resname_orig} at chain "
                            f"{chain_id}, input residue {orig.resseq_orig} contains "
                            f"unsupported atom name(s): {','.join(unknown)}.",
                            "Use AMBER OL3/OL15-compatible nucleic-acid atom names.",
                        )
                    )
                    continue
                if suffix:
                    log.append(
                        f"Terminal residue: {orig.resname_orig}{orig.resseq_orig} "
                        f"-> {resname_out}."
                    )
            elif base in PROTEIN_NAMES:
                resname_out = base
            elif base in SOLVENT_AND_ADDITIVES:
                errors.append(
                    CompatibilityIssue(
                        f"unsupported-solvent:{base}",
                        f"Solvent or crystallisation additive {base} at chain "
                        f"{chain_id}, input residue {orig.resseq_orig} has no "
                        "ONEBEAD mapping.",
                        "Use --drop-solvent to remove water and common "
                        "crystallisation additives.",
                    )
                )
                continue
            else:
                errors.append(
                    CompatibilityIssue(
                        f"unsupported-residue:{base}",
                        f"Residue or moiety {orig.resname_orig} at chain {chain_id}, "
                        f"input residue {orig.resseq_orig} has no implemented "
                        "ONEBEAD mapping.",
                        "Remove it from the PLUMED ATOMS/TEMPLATE selection or add "
                        "a validated mapping to both pdb2plmd and SAXS.cpp.",
                    )
                )
                continue

            if glycan is not None:
                output_names = [
                    normalize_glycan_atom_name(atom.atom_name, resname_out)
                    for atom in res_atoms
                ]
            else:
                output_names = normalize_atom_names_for_residue(res_atoms, resname_out)

            duplicates = sorted(
                name for name, count in Counter(output_names).items() if count > 1
            )
            if duplicates:
                errors.append(
                    CompatibilityIssue(
                        f"duplicate-atom-names:{resname_out}:{','.join(duplicates)}",
                        f"Duplicate output atom names in {resname_out}{residue_index} "
                        f"chain {chain_id}: {','.join(duplicates)}.",
                        "Check the input naming and residue-specific atom-name "
                        "conversion before using the template.",
                    )
                )
                continue

            for atom, atom_out in zip(res_atoms, output_names):
                if len(atom_out) > 4:
                    errors.append(
                        CompatibilityIssue(
                            f"atom-name-too-long:{atom_out}",
                            f"Atom name {atom_out!r} exceeds four PDB columns at "
                            f"input line {atom.line_number}.",
                        )
                    )
                    continue
                if atom_out != normalize_atom_name(atom.atom_name):
                    log.append(
                        f"Atom name: {normalize_atom_name(atom.atom_name)} -> {atom_out} "
                        f"in {resname_out}{residue_index}."
                    )
                element_out = atom.element
                if ion is not None:
                    expected_element = ION_ELEMENT[resname_out]
                    if element_out.upper() != expected_element:
                        log.append(
                            f"Ion element: {element_out or '?'} -> {expected_element} "
                            f"in {resname_out}{residue_index}."
                        )
                    element_out = expected_element
                converted.append(
                    replace(
                        atom,
                        chain_out=chain_id,
                        resseq_out=residue_index,
                        resname_out=resname_out,
                        atom_name_out=atom_out,
                        element=element_out,
                    )
                )

    stop_for_compatibility_issues(errors)
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
        out.write(f"REMARK pdb2plmd version {VERSION}\n")
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


def default_log_path(output_path: str) -> str:
    """Return the automatic log name used by both verbose and error runs."""
    p = Path(output_path)
    if p.suffix:
        return str(p.with_suffix(".log"))
    return str(Path(str(p) + ".log"))


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def atomic_write_pdb(atoms: List[AtomRecord], path: str) -> None:
    """Write a complete PDB to a temporary sibling and replace only on success."""
    dest = Path(path)
    tmp = dest.with_name(f".{dest.name}.pdb2plmd.{os.getpid()}.tmp")
    try:
        write_pdb(atoms, str(tmp))
        os.replace(tmp, dest)
    finally:
        try:
            if tmp.exists():
                tmp.unlink()
        except OSError:
            pass


def write_success_log(
    log_path: str,
    args: argparse.Namespace,
    all_atoms: List[AtomRecord],
    selected_indices: List[int],
    converted: List[AtomRecord],
    log_lines: List[str],
) -> None:
    with open(log_path, 'w', encoding='utf-8') as log:
        log.write(f"pdb2plmd {VERSION} verbose log\n")
        log.write("=" * 32 + "\n")
        log.write(f"Timestamp UTC: {utc_timestamp()}\n")
        log.write("Status: OK\n")
        log.write(f"Python: {sys.executable}\n")
        log.write(f"Command: {shlex.join(sys.argv)}\n")
        log.write(f"Input:  {args.input}\n")
        log.write(f"Output: {args.output}\n")
        if args.serials is not None:
            log.write(f"PDB serial range expression: {args.serials}\n")
            log.write("Selection mode: PDB atom serial (-s/--serials)\n")
        else:
            log.write(f"Atom-order range expression: {args.atoms if args.atoms is not None else 'all'}\n")
            log.write("Selection mode: ATOM/HETATM record order (-a/--atoms)\n")
        log.write(f"MODEL: {args.model}\n")
        log.write("Input-format handling: automatic preflight + record-aware parsing\n")
        if args.charmm is not None:
            log.write(f"Legacy format override: {'CHARMM forced' if args.charmm else 'CHARMM disabled'}\n")
        log.write(f"Split on gaps: {args.split_on_gaps}\n")
        log.write(f"AltLoc selection: {args.altloc}\n")
        log.write(f"Drop solvent: {args.drop_solvent}\n")
        log.write(f"Total ATOM/HETATM records after MODEL/altLoc/solvent filtering: {len(all_atoms)}\n")
        log.write(f"Selected atoms: {len(selected_indices)}\n")
        log.write(f"Output atoms: {len(converted)}\n")
        if selected_indices:
            log.write(f"Selected input atom-order range: {selected_indices[0]}..{selected_indices[-1]}\n")
        log.write("\nNotes:\n")
        if args.serials is not None:
            log.write("- The -s range is applied to the PDB atom-serial field after MODEL, altLoc and solvent filtering.\n")
            log.write("- TER and other non-coordinate records are not selected; repeated selected atom serials are rejected as ambiguous.\n")
        else:
            log.write("- The -a range is applied after MODEL, altLoc and solvent filtering.\n")
            log.write("- It refers to ATOM/HETATM order, not PDB serial.\n")
        log.write("- Output atom order is identical to the selected input atom order.\n")
        log.write("- Output atom serials are renumbered sequentially; this does not change PLUMED atom order.\n")
        log.write("- Residues are renumbered sequentially within each output chain to avoid PLUMED residue-range gaps.\n")
        log.write("- Residue/moiety compatibility check: PASS for all implemented converter rules.\n")
        if log_lines:
            log.write("\nConversion details:\n")
            for line in log_lines:
                log.write(line + "\n")


def write_error_log(
    log_path: str,
    args: argparse.Namespace,
    error_text: str,
    context_lines: Sequence[str],
    unexpected_traceback: str = "",
    output_existed_before: bool = False,
) -> str:
    """Write a complete failure report and return the path actually used."""
    requested = Path(log_path)
    candidates = [requested]
    fallback = Path.cwd() / requested.name
    if fallback != requested:
        candidates.append(fallback)

    last_error: Optional[Exception] = None
    for target in candidates:
        try:
            with target.open('w', encoding='utf-8') as log:
                log.write(f"pdb2plmd {VERSION} error log\n")
                log.write("=" * 30 + "\n")
                log.write(f"Timestamp UTC: {utc_timestamp()}\n")
                log.write("Status: ERROR\n")
                log.write(f"Python: {sys.executable}\n")
                log.write(f"Command: {shlex.join(sys.argv)}\n")
                log.write(f"Input:  {args.input}\n")
                log.write(f"Output requested: {args.output}\n")
                log.write(f"Verbose requested: {bool(args.verbose or args.log is not None)}\n")
                log.write(f"Output existed before run: {output_existed_before}\n")
                log.write("Output from this run: NOT WRITTEN\n")
                if output_existed_before:
                    log.write("Existing output file: LEFT UNCHANGED\n")
                if context_lines:
                    log.write("\nContext before failure:\n")
                    for line in context_lines:
                        log.write(line + "\n")
                log.write("\nError/details:\n")
                log.write((error_text or "Unknown conversion error").rstrip() + "\n")
                if unexpected_traceback:
                    log.write("\nUnexpected exception traceback:\n")
                    log.write(unexpected_traceback.rstrip() + "\n")
            return str(target)
        except OSError as exc:
            last_error = exc
    raise OSError(f"could not write error log {requested}: {last_error}")


def run_conversion(args: argparse.Namespace) -> tuple[
    List[AtomRecord], List[int], List[AtomRecord], List[str]
]:
    preflight = preflight_pdb_format(args.input)
    context_lines: List[str] = [
        f"Input format preflight: {preflight.style}.",
        (
            "Preflight sample: "
            f"{preflight.sampled_atoms} atom records; "
            f"chain column populated in {preflight.chain_records}; "
            f"SEGID populated in {preflight.segid_records}; "
            f"four-character residue extension observed in "
            f"{preflight.four_char_residue_records}."
        ),
        "Residue-name parsing: automatic per record (3- or 4-character PDB/CHARMM layout).",
        "Chain handling: standard PDB chain ID when present; otherwise SEGID-derived fallback.",
    ]
    if args.charmm is not None:
        context_lines.append(
            "Legacy format override: "
            + ("CHARMM forced." if args.charmm else "CHARMM handling disabled.")
        )

    all_atoms, parse_notes = parse_pdb(
        args.input,
        model=args.model,
        charmm=args.charmm,
        drop_solvent=args.drop_solvent,
        altloc=args.altloc,
    )
    context_lines.extend(parse_notes)
    if not all_atoms:
        raise SystemExit(f"No ATOM/HETATM records found in {args.input}")

    if args.serials is not None:
        selected_indices = parse_serial_selection(args.serials, all_atoms)
    else:
        atom_expr = args.atoms if args.atoms is not None else "all"
        selected_indices = parse_range(atom_expr, len(all_atoms))
    index_set = set(selected_indices)
    selected_atoms = [a for a in all_atoms if a.input_atom_index in index_set]
    if not selected_atoms:
        raise SystemExit("The atom selection is empty.")

    boundary_notes, serial_mismatch = selection_boundary_notes(
        all_atoms, selected_indices
    )
    context_lines.extend(boundary_notes)
    if (
        args.serials is None
        and serial_mismatch
        and (args.atoms or "all").strip().lower() not in {"all", "*"}
    ):
        print(
            "WARNING: -a/--atoms uses ATOM/HETATM record order, not PDB serials; "
            f"the selected end is {boundary_notes[1].removeprefix('Selection end: ').rstrip('.')}.",
            file=sys.stderr,
        )

    selection_issues = partial_selection_issues(all_atoms, selected_indices)
    if selection_issues:
        lines = [
            "ERROR: atom selection leaves one or more incomplete residues.",
            *(f"- {issue}" for issue in selection_issues),
            "Select complete amino-acid/nucleic-acid residues (and complete "
            "multi-atom moieties) only.",
        ]
        if boundary_notes:
            lines.append("Selection boundary:")
            lines.extend(f"- {note}" for note in boundary_notes)
        raise SystemExit("\n".join(lines))

    try:
        converted, conversion_lines = convert_atoms(
            selected_atoms,
            split_on_gaps=args.split_on_gaps,
        )
    except SystemExit as exc:
        error_text = str(exc) if exc.code not in (None, "") else "Conversion failed."
        if boundary_notes:
            error_text += "\n\nSelection boundary:\n" + "\n".join(
                f"- {note}" for note in boundary_notes
            )
        raise SystemExit(error_text) from None
    context_lines.extend(conversion_lines)

    if any(
        not math.isfinite(value)
        for atom in converted
        for value in (atom.x, atom.y, atom.z, atom.occ, atom.bfac)
    ):
        raise SystemExit("The selected atoms contain non-finite numeric values.")

    n_h = sum(1 for a in converted if a.element == "H")
    if converted and n_h < 0.2 * len(converted):
        msg = (
            f"WARNING: only {n_h}/{len(converted)} atoms are hydrogens "
            f"({100.0 * n_h / len(converted):.1f}%). The ONEBEAD template is "
            f"expected to be all-atom."
        )
        print(msg, file=sys.stderr)
        context_lines.append(msg)

    return all_atoms, selected_indices, converted, context_lines


def main() -> int:
    args = parse_args()
    if args.log is not None:
        # Preserve old -g/--log behavior while keeping it out of the public help.
        args.verbose = True

    auto_log = default_log_path(args.output)
    if args.log in (None, "__AUTO__"):
        log_path = auto_log
    else:
        log_path = args.log

    output_existed_before = Path(args.output).exists()
    context_lines: List[str] = []
    try:
        all_atoms, selected_indices, converted, context_lines = run_conversion(args)
        # In verbose mode prepare the log first. If the final atomic PDB write
        # fails, the exception handler replaces this with an error log, so a
        # failed run never leaves a newly written PDB behind.
        if args.verbose:
            write_success_log(
                log_path, args, all_atoms, selected_indices, converted, context_lines
            )
        atomic_write_pdb(converted, args.output)
        print(
            f"pdb2plmd {VERSION}: OK - wrote {args.output} ({len(converted)} atoms)"
            + (f"; log {log_path}" if args.verbose else "")
        )
        return 0
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
        if code == 0:
            return 0
        error_text = str(exc) if exc.code not in (None, "") else "Conversion failed."
        try:
            actual_log = write_error_log(
                log_path,
                args,
                error_text,
                context_lines,
                output_existed_before=output_existed_before,
            )
            first_error = error_text.splitlines()[0]
            if first_error.upper().startswith("ERROR:"):
                first_error = first_error.split(":", 1)[1].strip()
            print(f"pdb2plmd {VERSION}: ERROR - {first_error}", file=sys.stderr)
            print(f"Details: {actual_log}", file=sys.stderr)
        except OSError as log_exc:
            print(f"pdb2plmd {VERSION}: ERROR - {error_text}", file=sys.stderr)
            print(f"Additionally failed to write the error log: {log_exc}", file=sys.stderr)
        return code or 1
    except Exception as exc:
        tb = traceback.format_exc()
        error_text = f"Unexpected {type(exc).__name__}: {exc}"
        try:
            actual_log = write_error_log(
                log_path,
                args,
                error_text,
                context_lines,
                unexpected_traceback=tb,
                output_existed_before=output_existed_before,
            )
            print(f"pdb2plmd {VERSION}: ERROR - {error_text}", file=sys.stderr)
            print(f"Details: {actual_log}", file=sys.stderr)
        except OSError as log_exc:
            print(f"pdb2plmd {VERSION}: ERROR - {error_text}", file=sys.stderr)
            print(f"Additionally failed to write the error log: {log_exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
