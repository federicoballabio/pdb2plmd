from __future__ import annotations
from typing import List, Sequence, Tuple, Optional
import re
from .models import AtomRecord, CompatibilityIssue
from .constants import PENTOSE_ATOMS, BASE_ATOMS, PHOSPHATE_ATOMS, KNOWN_ONEBEAD_NUC_ATOMS
from .pdbio import base_resname, normalize_atom_names_for_residue

_CHARMM_BASES = {"ADE":"A", "CYT":"C", "GUA":"G", "THY":"T", "URA":"U"}
_AMBER_OLD_BASES = {"ADE":"A", "CYT":"C", "GUA":"G", "THY":"T", "URA":"U"}
_AMBER_LEGACY_ATOM_MAP = {"O1'":"O4'", "OA":"OP1", "OB":"OP2"}
_THYMINE_MAP = {"C5M":"C7", "H51":"H71", "H52":"H72", "H53":"H73"}

_BASE_REQUIRED_HEAVY = {
    "A": {"N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"},
    "C": {"N1","C2","O2","N3","C4","N4","C5","C6"},
    "G": {"N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"},
    "U": {"N1","C2","O2","N3","C4","O4","C5","C6"},
    "DA": {"N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"},
    "DC": {"N1","C2","O2","N3","C4","N4","C5","C6"},
    "DG": {"N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"},
    "DT": {"N1","C2","O2","N3","C4","O4","C5","C6","C7"},
}
_RNA_SUGAR_REQUIRED = {"O5'","C5'","O4'","C4'","O3'","C3'","O2'","C2'","C1'"}
_DNA_SUGAR_REQUIRED = {"O5'","C5'","O4'","C4'","O3'","C3'","C2'","C1'"}


def resolve_nucleic_base(resname_orig: str, atom_names, convention: str) -> Tuple[Optional[str], Optional[CompatibilityIssue], List[str]]:
    raw = resname_orig.strip().upper()
    names = {n.strip().replace("*", "'") for n in atom_names}
    notes: List[str] = []

    if convention in {"amber", "gromacs-amber"}:



        m = re.fullmatch(r"R([ACGU])(?:[35])?", raw)
        if m:
            base = m.group(1)
            notes.append(f"Legacy AMBER RNA residue: {raw} -> {base}.")
            return base, None, notes




        if raw in _AMBER_OLD_BASES:
            has_o2 = "O2'" in names
            if raw in {"ADE","CYT","GUA"}:
                base = {"ADE":"A","CYT":"C","GUA":"G"}[raw]
                resolved = base if has_o2 else "D" + base
                notes.append(f"Legacy AMBER residue: {raw} -> {resolved} from O2' chemistry.")
                return resolved, None, notes
            if raw == "THY":
                if has_o2:
                    return None, CompatibilityIssue(
                        "unsupported-amber-ribothymidine",
                        "AMBER THY with O2' present is a ribothymidine-like state not supported by the current ONEBEAD target.",
                    ), notes
                notes.append("Legacy AMBER residue: THY -> DT.")
                return "DT", None, notes
            if raw == "URA":
                if not has_o2:
                    return None, CompatibilityIssue(
                        "unsupported-amber-deoxyuridine",
                        "AMBER URA without O2' is a deoxyuridine-like state not supported by the current ONEBEAD target.",
                    ), notes
                notes.append("Legacy AMBER residue: URA -> U.")
                return "U", None, notes

    if convention in {"charmm", "gromacs-charmm"} and raw in _CHARMM_BASES:
        has_o2 = "O2'" in names
        if raw in {"ADE","CYT","GUA"}:
            base = {"ADE":"A","CYT":"C","GUA":"G"}[raw]
            return (base if has_o2 else "D" + base), None, notes
        if raw == "THY":
            if has_o2:
                return None, CompatibilityIssue(
                    "unsupported-charmm-ribothymidine",
                    "CHARMM THY with O2' present is a ribothymidine-like state not supported by the current ONEBEAD target.",
                ), notes
            return "DT", None, notes
        if raw == "URA":
            if not has_o2:
                return None, CompatibilityIssue(
                    "unsupported-charmm-deoxyuridine",
                    "CHARMM URA without O2' is a deoxyuridine-like state not supported by the current ONEBEAD target.",
                ), notes
            return "U", None, notes
    generic = base_resname(raw)
    if generic in {"A","C","G","U","DA","DC","DG","DT"}:
        return generic, None, notes
    return None, None, notes


def source_terminal_suffix(resname_orig: str, convention: str) -> str:
    raw = resname_orig.strip().upper()
    if convention in {"amber", "gromacs-amber"} and re.fullmatch(r"R[ACGU][35]", raw):
        return raw[-1]
    return ""


def normalize_nucleic_atom_names(res_atoms: Sequence[AtomRecord], resname_out: str, convention: str) -> List[str]:
    names = normalize_atom_names_for_residue(res_atoms, resname_out)
    if convention in {"amber", "gromacs-amber"}:



        names = [_AMBER_LEGACY_ATOM_MAP.get(n, n) for n in names]
    if convention in {"charmm", "amber", "gromacs-charmm", "gromacs-amber", "gromacs-gromos"} and resname_out.startswith("DT"):
        names = [_THYMINE_MAP.get(n, n) for n in names]
    return names


def validate_nucleic(res_atoms: Sequence[AtomRecord], base: str, suffix: str, output_names: Sequence[str]):
    fatal: List[CompatibilityIssue] = []
    readiness: List[str] = []
    unknown = sorted({n for n in output_names if n not in KNOWN_ONEBEAD_NUC_ATOMS})
    if unknown:
        first = res_atoms[0]
        fatal.append(CompatibilityIssue(
            f"nucleic-atom-names:{first.resname_orig}:{','.join(unknown)}",
            f"Nucleic-acid residue {first.resname_orig}, input residue {first.resseq_orig} contains unsupported atom name(s): {','.join(unknown)}.",
            "Use AMBER OL3/OL15-compatible target atom names after source-convention conversion.",
        ))
        return fatal, readiness

    sugar = {n for n in output_names if n in PENTOSE_ATOMS}
    base_atoms = {n for n in output_names if n in BASE_ATOMS}
    phosphate = {n for n in output_names if n in PHOSPHATE_ATOMS}
    if not sugar:
        fatal.append(CompatibilityIssue("nucleic-empty-sugar", "Nucleotide has no atoms assignable to the ONEBEAD sugar bead."))
    if not base_atoms:
        fatal.append(CompatibilityIssue("nucleic-empty-base", "Nucleotide has no atoms assignable to the ONEBEAD base bead."))
    if suffix != "5" and not phosphate:
        fatal.append(CompatibilityIssue("nucleic-empty-phosphate", "Nucleotide has no atoms assignable to the ONEBEAD phosphate bead; only an identified 5'-OH terminal state may omit it."))
    if fatal:
        return fatal, readiness

    heavy = {
        name for atom, name in zip(res_atoms, output_names)
        if atom.element.upper() not in {"H","D","T"}
    }
    sugar_required = _DNA_SUGAR_REQUIRED if base.startswith("D") else _RNA_SUGAR_REQUIRED
    missing_sugar = sorted(sugar_required - heavy)
    missing_base = sorted(_BASE_REQUIRED_HEAVY[base] - heavy)
    missing_phos = []
    if suffix != "5":
        if "P" not in heavy:
            missing_phos.append("P")
        if not ({"OP1","O1P"} & heavy):
            missing_phos.append("OP1/O1P")
        if not ({"OP2","O2P"} & heavy):
            missing_phos.append("OP2/O2P")
    if missing_sugar or missing_base or missing_phos:
        parts = []
        if missing_sugar: parts.append("sugar=" + ",".join(missing_sugar))
        if missing_base: parts.append("base=" + ",".join(missing_base))
        if missing_phos: parts.append("phosphate=" + ",".join(missing_phos))
        readiness.append("missing required nucleotide heavy atom(s): " + "; ".join(parts))
    return fatal, readiness
