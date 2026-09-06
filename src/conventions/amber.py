from __future__ import annotations
from typing import List, Optional, Tuple
from ..models import CompatibilityIssue

_STANDARD = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","ILE","LEU",
    "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
}
_STATE = {"HID":"HID","HIE":"HIE","HIP":"HIP","CYX":"CYX"}
_AMBER_TERMINAL_PARENT = {
    **{f"N{x}": x for x in _STANDARD},
    **{f"C{x}": x for x in _STANDARD},
    "NHID":"HID", "NHIE":"HIE", "NHIP":"HIP", "NCYX":"CYX",
    "CHID":"HID", "CHIE":"HIE", "CHIP":"HIP", "CCYX":"CYX",
}
_UNSUPPORTED_STATES = {
    "CYM": "AMBER CYM is a deprotonated cysteine thiolate; the current ONEBEAD target has no dedicated CYM parameter.",
    "LYN": "AMBER LYN is neutral lysine; the current ONEBEAD target has no dedicated neutral-lysine parameter.",
}
_CAPS = {"ACE", "NME", "NHE"}
_MONOMER_NUC = {"AN","CN","GN","UN","DAN","DCN","DGN","DTN","RAN","RCN","RGN","RUN"}


def resolve_amber_residue(raw: str, is_first: bool, is_last: bool) -> Tuple[str, Optional[CompatibilityIssue], List[str]]:
    r = raw.strip().upper()
    notes: List[str] = []
    if r in _UNSUPPORTED_STATES:
        return r, CompatibilityIssue(
            f"amber-state:{r}", _UNSUPPORTED_STATES[r],
            "Implement and validate a dedicated ONEBEAD parameter before using this state.",
        ), notes
    if r in _CAPS:
        return r, CompatibilityIssue(
            f"amber-cap:{r}",
            f"AMBER terminal capping group {r} has no ONEBEAD bead.",
            "Remove the cap from the PLUMED template selection or implement a dedicated validated bead.",
        ), notes
    if r in _MONOMER_NUC:
        return r, CompatibilityIssue(
            f"amber-nucleic-monomer:{r}",
            f"AMBER nucleic-acid monomer building block {r} contains both 5'-OH and 3'-OH terminal chemistry, while the current SAXS ONEBEAD target has no dedicated monomer bead/state.",
            "Do not collapse the monomer to a 5' or 3' terminal bead without a dedicated validated ONEBEAD parameter/state.",
        ), notes
    if r in _AMBER_TERMINAL_PARENT:
        expected_n = r.startswith("N") and not r.startswith("NME")
        expected_c = r.startswith("C")
        if expected_n and not is_first:
            return r, CompatibilityIssue(
                f"amber-terminal-position:{r}",
                f"AMBER N-terminal building block {r} occurs away from the first residue of its inferred chain.",
                "Check chain/TER handling or correct the source residue naming.",
            ), notes
        if expected_c and not is_last:
            return r, CompatibilityIssue(
                f"amber-terminal-position:{r}",
                f"AMBER C-terminal building block {r} occurs away from the last residue of its inferred chain.",
                "Check chain/TER handling or correct the source residue naming.",
            ), notes
        parent = _AMBER_TERMINAL_PARENT[r]
        notes.append(f"AMBER terminal building block: {r} -> {parent}.")
        return parent, None, notes
    return r, None, notes
