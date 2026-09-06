from __future__ import annotations
from .constants import (
    _GLC_FAMILY, _GAL_FAMILY, _MAN_A, _MAN_B, _FUC_FAMILY, _NEU_FAMILY,
    _NAG_NAMES, _NGA_NAMES, GLYCAN_ATOM_MAP_HEX, GLYCAN_ATOM_MAP_SIA,
    _GLYCAN_REQUIRED_HEAVY, _GLYCAN_ALLOWED_HEAVY,
)






_GLYCAM_NAC_ATOM_MAP = {
    "C2N": "C7", "O2N": "O7", "CME": "C8", "H2N": "HN2",
    "H1M": "H81", "H2M": "H82", "H3M": "H83",
}
_GLYCAM_SIA_ATOM_MAP = {
    "C5N": "C10", "O5N": "O10", "CME": "C11", "H5N": "HN5",
    "H1M": "H111", "H2M": "H112", "H3M": "H113",
}
_GLYCAM_HYDROXY_H = {
    "H1O": "HO1", "H2O": "HO2", "H3O": "HO3", "H4O": "HO4",
    "H6O": "HO6", "H7O": "HO7", "H8O": "HO8", "H9O": "HO9",
}


def has_n_acetyl_signature(atom_names) -> bool:
    names = {str(a).strip().upper() for a in atom_names}
    charmm = {"N", "C", "O", "CT"}.issubset(names)
    ccd = {"N2", "C7", "O7", "C8"}.issubset(names)
    glycam = {"N2", "C2N", "O2N", "CME"}.issubset(names)
    return charmm or ccd or glycam


def canonical_glycan(resname, atom_names):
    r = resname.strip().upper()
    has_n_acetyl = has_n_acetyl_signature(atom_names)
    if r in _NAG_NAMES:
        return "NAG"
    if r in _NGA_NAMES:
        return "NGA"
    if r in _GLC_FAMILY:
        return "NAG" if has_n_acetyl else "GLC"
    if r in _GAL_FAMILY:
        return "NGA" if has_n_acetyl else "GAL"
    if r in _MAN_A:
        return "MAN"
    if r in _MAN_B:
        return "BMA"
    if r in _FUC_FAMILY:
        return "FUC"
    if r in _NEU_FAMILY:
        return "SIA"
    return None


def normalize_glycan_atom_name(name, resname_out, input_convention="generic"):
    n = name.strip()
    if input_convention == "glycam":
        n = _GLYCAM_HYDROXY_H.get(n, n)
        if resname_out in {"NAG", "NGA"}:
            n = _GLYCAM_NAC_ATOM_MAP.get(n, n)
        elif resname_out == "SIA":
            n = _GLYCAM_SIA_ATOM_MAP.get(n, n)

    if resname_out == "SIA":
        return GLYCAN_ATOM_MAP_SIA.get(n, n)
    if resname_out in {"NAG", "NGA"}:
        return GLYCAN_ATOM_MAP_HEX.get(n, n)
    return n


def validate_glycan_heavy_atoms(res_atoms, glycan, input_convention="generic"):
    normalized = [normalize_glycan_atom_name(a.atom_name, glycan, input_convention) for a in res_atoms]
    heavy = {
        name for atom, name in zip(res_atoms, normalized)
        if atom.element.upper() not in {"H", "D", "T"}
    }
    required = _GLYCAN_REQUIRED_HEAVY[glycan]
    allowed = _GLYCAN_ALLOWED_HEAVY[glycan]
    missing = sorted(required - heavy)
    extra = sorted(heavy - allowed)
    fatal = []
    readiness = []
    if missing:
        readiness.append((
            "incomplete-glycan",
            "missing required heavy atom(s): " + ",".join(missing),
        ))
    if extra:
        fatal.append((
            "modified-glycan-heavy-atoms",
            "unexpected heavy atom name(s): " + ",".join(extra),
        ))
    notes = []
    if glycan != "SIA" and "O1" in heavy:
        notes.append(
            "Free anomeric O1 is present; accepted as a free/reducing-end state, "
            "but this residue state is not the covalently linked donor state used "
            "for the glycan parameter training at a glycosidic junction."
        )
    if glycan == "SIA" and "O2" in heavy:
        notes.append("Free C2 O2 is present; accepted as an optional free-state atom for SIA.")
    return fatal, readiness, notes
