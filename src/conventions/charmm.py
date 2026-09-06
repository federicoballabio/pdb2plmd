from __future__ import annotations
from typing import Sequence, Tuple, List
from ..models import AtomRecord, PdbFormatPreflight

_CHARMM_HIS = {"HSD","HSE","HSP"}
_CHARMM_CARBS = {
    "AGLC","BGLC","AGLCNA","BGLCNA","AGAL","BGAL","AGALNA","BGALNA",
    "AMAN","BMAN","AFUC","BFUC","ANE5","BNE5","ANE5AC","BNE5AC",
}
_CHARMM_TERMINAL_ATOMS = {"HT1","HT2","HT3","OT1","OT2"}


def detect_charmm(atoms: Sequence[AtomRecord], preflight: PdbFormatPreflight) -> Tuple[bool, List[str]]:
    evidence_parts: List[str] = []
    resnames = {a.resname_orig.strip().upper() for a in atoms}
    hist = sorted(resnames & _CHARMM_HIS)
    carbs = sorted(resnames & _CHARMM_CARBS)




    if preflight.four_char_residue_records and carbs:
        evidence_parts.append(
            f"four-character residue names observed in {preflight.four_char_residue_records} sampled record(s)"
        )
    if preflight.sampled_atoms and preflight.chain_records == 0 and preflight.segid_records >= max(1, int(0.8 * preflight.sampled_atoms)):
        evidence_parts.append("SEGID-based coordinate layout with an empty standard chain column")
    atomnames = {a.atom_name.strip().upper() for a in atoms}
    terminal = sorted(atomnames & _CHARMM_TERMINAL_ATOMS)
    if hist:
        evidence_parts.append("CHARMM histidine residue name(s): " + ",".join(hist))
    if carbs:
        evidence_parts.append("CHARMM carbohydrate residue name(s): " + ",".join(carbs))
    if terminal:
        evidence_parts.append("CHARMM terminal atom name(s): " + ",".join(terminal))
    if not evidence_parts:
        return False, []
    return True, ["Auto-detected CHARMM/CHARMM-GUI input: " + "; ".join(evidence_parts) + "."]
