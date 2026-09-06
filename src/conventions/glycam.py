from __future__ import annotations
from dataclasses import dataclass
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from ..models import CompatibilityIssue

_DATA = Path(__file__).resolve().parents[1] / "data" / "GLYCAM_CODE_INDEX_C4.tsv"

_LINKAGE_POSITIONS = {
    "0": "terminal/unsubstituted",
    "P": "2,3,4,6", "Q": "3,4,6", "R": "2,4,6", "S": "2,3,6",
    "T": "2,3,4", "U": "4,6", "V": "3,6", "W": "3,4",
    "X": "2,6", "Y": "2,4", "Z": "2,3",
}

@dataclass(frozen=True)
class GlycamResidueState:
    code: str
    ccd_comp_id: str
    anomer: str
    config: str
    sugar: str
    ring: str
    linkage_code: str
    method: str
    onebead_target: str
    identity: str
    status: str

    @property
    def linkage_description(self) -> str:
        if self.linkage_code in _LINKAGE_POSITIONS:
            return _LINKAGE_POSITIONS[self.linkage_code]
        if self.linkage_code.isdigit():
            return self.linkage_code
        return "unexpanded:" + self.linkage_code


def _load_index() -> Dict[str, GlycamResidueState]:
    out: Dict[str, GlycamResidueState] = {}
    with _DATA.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            state = GlycamResidueState(
                code=row["code"], ccd_comp_id=row["ccd_comp_id"], anomer=row["anomer"],
                config=row["config"], sugar=row["sugar"], ring=row["ring"],
                linkage_code=row["linkage"], method=row["method"],
                onebead_target=row["onebead_target"], identity=row["identity"],
                status=row["status"],
            )
            out[state.code] = state
    return out

GLYCAM_CODE_INDEX = _load_index()



_GLYCAM_DEFERRED_ATTACHMENT_STATES = {
    "ZOLS": "GLYCAM ZOLS is a zwitterionic terminal O-linked serine variant; the current ONEBEAD target has no dedicated terminal-state bead.",
    "ZOLT": "GLYCAM ZOLT is a zwitterionic terminal O-linked threonine variant; the current ONEBEAD target has no dedicated terminal-state bead.",
}


def resolve_glycam_residue(raw: str) -> Tuple[str, Optional[CompatibilityIssue], List[str], Optional[GlycamResidueState]]:
    r = raw.strip()
    notes: List[str] = []
    ru = r.upper()
    if ru in _GLYCAM_DEFERRED_ATTACHMENT_STATES:
        return r, CompatibilityIssue(
            f"glycam-attachment-state:{ru}", _GLYCAM_DEFERRED_ATTACHMENT_STATES[ru],
            "Do not flatten this terminal state to OLS/OLT without a dedicated validated ONEBEAD state.",
        ), notes, None


    state = GLYCAM_CODE_INDEX.get(r)
    if state is None:
        return raw, None, notes, None
    if state.status != "SUPPORTED_ONEBEAD" or not state.onebead_target:
        ident = state.identity or f"{state.anomer}-{state.config}-{state.sugar} ({state.ring})"
        return raw, CompatibilityIssue(
            f"glycam-unsupported:{state.code}",
            f"GLYCAM residue {state.code} is a recognized {ident} state (linkage code {state.linkage_code}) but has no validated ONEBEAD glycan parameter.",
            "No chemically similar monosaccharide is substituted automatically.",
        ), notes, state
    notes.append(
        "GLYCAM residue state: "
        f"{state.code} -> {state.onebead_target}; identity={state.identity or state.sugar}; "
        f"anomer={state.anomer}; config={state.config}; ring={state.ring}; "
        f"linkage_code={state.linkage_code}; linkage={state.linkage_description}; "
        f"CCD={state.ccd_comp_id or 'unresolved'}."
    )
    return state.onebead_target, None, notes, state



def normalize_glycam_protein_atom_names(res_atoms, resname_out: str) -> List[str]:
    out: List[str] = []
    for atom in res_atoms:
        n = atom.atom_name.strip().replace("*", "'")
        if resname_out == "ASN" and n == "Cg":
            n = "CG"
        out.append(n)
    return out
