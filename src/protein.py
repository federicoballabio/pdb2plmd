from __future__ import annotations
from typing import Sequence, Tuple, List
from .models import AtomRecord

PROTEIN_REQUIRED_HEAVY = {
    "ALA": {"N","CA","C","O","CB"},
    "ARG": {"N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2"},
    "ASN": {"N","CA","C","O","CB","CG","OD1","ND2"},
    "ASP": {"N","CA","C","O","CB","CG","OD1","OD2"},
    "CYS": {"N","CA","C","O","CB","SG"},
    "GLN": {"N","CA","C","O","CB","CG","CD","OE1","NE2"},
    "GLU": {"N","CA","C","O","CB","CG","CD","OE1","OE2"},
    "GLY": {"N","CA","C","O"},
    "HIS": {"N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"},
    "ILE": {"N","CA","C","O","CB","CG1","CG2","CD1"},
    "LEU": {"N","CA","C","O","CB","CG","CD1","CD2"},
    "LYS": {"N","CA","C","O","CB","CG","CD","CE","NZ"},
    "MET": {"N","CA","C","O","CB","CG","SD","CE"},
    "PHE": {"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ"},
    "PRO": {"N","CA","C","O","CB","CG","CD"},
    "SER": {"N","CA","C","O","CB","OG"},
    "THR": {"N","CA","C","O","CB","OG1","CG2"},
    "TRP": {"N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"},
    "TYR": {"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH"},
    "VAL": {"N","CA","C","O","CB","CG1","CG2"},
}

_HIS_STATES = {"HIS","HID","HIE","HIP","HSD","HSE","HSP"}


def resolve_protein_state(base: str) -> Tuple[str, str, List[str]]:
    b = base.upper()
    warnings: List[str] = []
    if b == "GLH":
        warnings.append("Warning: GLH is approximated as GLU for SAXS ONEBEAD.")
        return "GLU", "GLU", warnings
    if b == "ASH":
        warnings.append("Warning: ASH is approximated as ASP for SAXS ONEBEAD.")
        return "ASP", "ASP", warnings
    if b in _HIS_STATES:
        return b, "HIS", warnings
    if b == "CYX":
        return "CYX", "CYS", warnings
    return b, b, warnings


def _completeness_name(atom: AtomRecord, chemistry: str, input_convention: str = "generic") -> str:
    name = atom.atom_name.strip().replace("*", "'")
    if input_convention == "glycam" and chemistry == "ASN" and name == "Cg":
        return "CG"
    if chemistry == "ILE" and name == "CD":
        return "CD1"
    if name == "OT1":
        return "O"
    if name == "OT2":
        return "OXT"
    if input_convention.startswith("gromacs"):
        if name in {"O1","O'","OC1"}:
            return "O"
        if name in {"O2","OT","O''","OC2"}:
            return "OXT"
    return name


def protein_missing_heavy(res_atoms: Sequence[AtomRecord], chemistry: str, input_convention: str = "generic") -> List[str]:
    required = PROTEIN_REQUIRED_HEAVY.get(chemistry)
    if required is None:
        return []
    present = {
        _completeness_name(a, chemistry, input_convention)
        for a in res_atoms
        if a.element.upper() not in {"H","D","T"}
    }
    return sorted(required - present)
