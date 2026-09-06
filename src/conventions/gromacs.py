from __future__ import annotations
from typing import List, Optional, Sequence, Tuple
from ..models import AtomRecord, CompatibilityIssue
from ..pdbio import normalize_atom_names_for_residue
from .amber import resolve_amber_residue

_GMX_CHARMM_RES = {
    "CYS2": "CYX",
    "ASPP": "ASH",
    "GLUP": "GLH",
}
_GMX_CHARMM_UNSUPPORTED = {
    "LSN": "GROMACS/CHARMM LSN is neutral lysine; the current ONEBEAD target has no dedicated neutral-lysine parameter.",
}
_GMX_CHARMM_CAPS = {"ACE", "CT3", "NH2"}




_GMX_OPLS_RES = {
    "CYSH": "CYS",
    "CYS2": "CYX",
    "ASPH": "ASH",
    "GLUH": "GLH",
    "HISD": "HID",
    "HISE": "HIE",
    "HISH": "HIP",
    "LYSH": "LYS",
}
_GMX_OPLS_UNSUPPORTED = {
    "ARGN": "GROMACS/OPLS-AA ARGN is neutral arginine; the current ONEBEAD target has no dedicated neutral-arginine parameter.",
    "LYSN": "GROMACS/OPLS-AA LYSN selects the neutral lysine building block; the current ONEBEAD target has no dedicated neutral-lysine parameter.",
    "AIB": "GROMACS/OPLS-AA AIB is a non-standard amino acid without a validated ONEBEAD parameter.",
    "PGLU": "GROMACS/OPLS-AA PGLU is pyroglutamate and has no validated ONEBEAD/LCPO mapping.",
}
_GMX_OPLS_CAPS = {"ACE", "NAC", "NH2", "NHE"}




_GMX_GROMOS_RES = {
    "CYSH": "CYS",
    "CYS2": "CYX",
    "ASPH": "ASH",
    "GLUH": "GLH",
    "HISA": "HID",
    "HISD": "HID",
    "HISB": "HIE",
    "HISE": "HIE",
    "HISH": "HIP",
    "LYSH": "LYS",
}
_GMX_GROMOS_UNSUPPORTED = {
    "ARGN": "GROMACS/GROMOS ARGN is neutral arginine; the current ONEBEAD target has no dedicated neutral-arginine parameter.",
    "LYSN": "GROMACS/GROMOS LYSN selects the neutral lysine building block; the current ONEBEAD target has no dedicated neutral-lysine parameter.",
    "CYS1": "GROMACS/GROMOS CYS1 is a specialized cysteine building block and is not assigned to the generic CYS/CYX ONEBEAD states.",
    "HIS1": "GROMACS/GROMOS HIS1 is a specialized histidine building block (used by GROMACS special-bond handling) and is not mapped automatically.",
    "HIS2": "GROMACS/GROMOS HIS2 is a specialized histidine building block and is not mapped automatically.",
    "HYP": "GROMACS/GROMOS HYP is hydroxyproline and has no validated PRO ONEBEAD/LCPO mapping.",
}
_GMX_GROMOS_CAPS = {"ACE", "NH2"}


def _issue(prefix: str, r: str, message: str, hint: str = "Implement and validate a dedicated ONEBEAD parameter before using this state."):
    return CompatibilityIssue(f"{prefix}:{r}", message, hint)


def resolve_gromacs_residue(
    raw: str,
    convention: str,
    is_first: bool,
    is_last: bool,
    atom_names: Optional[Sequence[str]] = None,
):
    r = raw.strip().upper()
    notes: List[str] = []
    if convention == "gromacs-amber":
        return resolve_amber_residue(r, is_first, is_last)

    if convention == "gromacs-charmm":
        if r in _GMX_CHARMM_UNSUPPORTED:
            return r, _issue("gromacs-charmm-state", r, _GMX_CHARMM_UNSUPPORTED[r]), notes
        if r in _GMX_CHARMM_CAPS:
            return r, _issue(
                "gromacs-charmm-cap", r,
                f"GROMACS/CHARMM terminal capping group {r} has no ONEBEAD bead.",
                "Remove the cap from the PLUMED template selection or implement a dedicated validated bead.",
            ), notes
        if r in _GMX_CHARMM_RES:
            parent = _GMX_CHARMM_RES[r]
            notes.append(f"GROMACS/CHARMM residue state: {r} -> {parent}.")
            return parent, None, notes
        return r, None, notes

    if convention == "gromacs-oplsaa":
        if r in _GMX_OPLS_UNSUPPORTED:
            return r, _issue("gromacs-oplsaa-state", r, _GMX_OPLS_UNSUPPORTED[r]), notes
        if r in _GMX_OPLS_CAPS:
            return r, _issue(
                "gromacs-oplsaa-cap", r,
                f"GROMACS/OPLS-AA terminal/capping group {r} has no ONEBEAD bead.",
                "Remove the cap from the PLUMED template selection or implement a dedicated validated bead.",
            ), notes
        if r in _GMX_OPLS_RES:
            parent = _GMX_OPLS_RES[r]
            notes.append(f"GROMACS/OPLS-AA residue state: {r} -> {parent}.")
            return parent, None, notes
        return r, None, notes

    if convention == "gromacs-gromos":
        if r in _GMX_GROMOS_UNSUPPORTED:
            return r, _issue("gromacs-gromos-state", r, _GMX_GROMOS_UNSUPPORTED[r]), notes
        if r in _GMX_GROMOS_CAPS:
            return r, _issue(
                "gromacs-gromos-cap", r,
                f"GROMACS/GROMOS terminal/capping group {r} has no ONEBEAD bead.",
                "Remove the cap from the PLUMED template selection or implement a dedicated validated bead.",
            ), notes
        if r in _GMX_GROMOS_RES:
            parent = _GMX_GROMOS_RES[r]
            notes.append(f"GROMACS/GROMOS residue state: {r} -> {parent}.")
            return parent, None, notes
        return r, None, notes

    return r, None, notes


def normalize_gromacs_protein_atom_names(
    res_atoms: Sequence[AtomRecord], resname_out: str, is_first: bool, is_last: bool
) -> List[str]:
    names = normalize_atom_names_for_residue(res_atoms, resname_out)
    out: List[str] = []
    for name in names:
        n = name
        if n == "HN":
            n = "H"
        if is_first:
            n = {"HT1":"H1", "HT2":"H2", "HT3":"H3"}.get(n, n)
        if is_last:
            n = {
                "O1":"O", "OT1":"O", "O'":"O", "OC1":"O",
                "O2":"OXT", "OT2":"OXT", "OT":"OXT", "O''":"OXT", "OC2":"OXT",
            }.get(n, n)
        if resname_out == "ILE":
            n = {
                "CD":"CD1",
                "HD1":"HD11", "HD2":"HD12", "HD3":"HD13",
            }.get(n, n)
        out.append(n)
    return out
