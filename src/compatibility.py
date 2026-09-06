from __future__ import annotations
from typing import Sequence, Tuple, List, Optional
from .models import AtomRecord, PdbFormatPreflight, CompatibilityReport
from .conventions import detect_charmm
from .conventions.glycam import GLYCAM_CODE_INDEX

_AMBER_STANDARD = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","ILE","LEU",
    "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
}
_AMBER_TERMINAL = {
    *{f"N{x}" for x in _AMBER_STANDARD},
    *{f"C{x}" for x in _AMBER_STANDARD},
    "NHID","NHIE","NHIP","NCYX","CHID","CHIE","CHIP","CCYX",
    "CYM","LYN",
}
_AMBER_NUCLEIC = {
    "RA","RC","RG","RU","RA5","RC5","RG5","RU5","RA3","RC3","RG3","RU3",
    "AN","CN","GN","UN","DAN","DCN","DGN","DTN","RAN","RCN","RGN","RUN",
}
_GMX_CHARMM_STRONG = {"ASPP","GLUP","LSN","CT3"}
_GMX_OPLS_STRONG = {"AIB","PGLU","NAC"}
_GMX_OPLS_GROMOS_SHARED = {"CYSH","ASPH","GLUH","HISD","HISE","HISH","LYSH"}
_EXPLICIT_ALIPHATIC_H = {
    "HA","HA1","HA2","HA3","HB1","HB2","HB3","HG1","HG2","HG3",
    "HD1","HD2","HD3","HE1","HE2","HE3",
}
_GLYCAM_STRONG_GROUPS = (
    {"C2N","O2N","CME"},
    {"C5N","O5N","CME"},
)
_GLYCAM_HYDROXY_H = {f"H{i}O" for i in range(1, 10)}


def _residue_groups(atoms: Sequence[AtomRecord]):
    groups = {}
    for atom in atoms:
        key = (atom.chain_orig, atom.segid, atom.resseq_orig, atom.icode_orig, atom.resname_orig)
        groups.setdefault(key, []).append(atom)
    return list(groups.values())


def _gromacs_atom_signatures(atoms: Sequence[AtomRecord]) -> List[str]:
    hits: List[str] = []
    for residue in _residue_groups(atoms):
        raw = residue[0].resname_orig.strip().upper()
        names = {a.atom_name.strip().upper() for a in residue}
        if {"OC1","OC2"}.issubset(names):
            hits.append(f"{raw}:OC1/OC2")
        if raw == "ILE" and "CD" in names:
            hd = sorted({"HD1","HD2","HD3"} & names)
            hits.append("ILE:CD" + ("/" + "-".join(hd) if hd else ""))
    return sorted(set(hits))


def _opls_all_atom_signatures(atoms: Sequence[AtomRecord]) -> List[str]:
    hits: List[str] = []
    for residue in _residue_groups(atoms):
        raw = residue[0].resname_orig.strip().upper()
        if raw not in _GMX_OPLS_GROMOS_SHARED:
            continue
        names = {a.atom_name.strip().upper() for a in residue}
        explicit = sorted(names & _EXPLICIT_ALIPHATIC_H)
        if len(explicit) >= 2:
            hits.append(f"{raw}:{','.join(explicit[:4])}")
    return hits


def _glycam_signatures(atoms: Sequence[AtomRecord]) -> List[str]:
    code_groups = []
    for residue in _residue_groups(atoms):
        code = residue[0].resname_orig.strip()
        if code not in GLYCAM_CODE_INDEX:
            continue
        names = {a.atom_name.strip().upper() for a in residue}
        lowercase = any(c.islower() for c in code)
        strong_atoms = any(group.issubset(names) for group in _GLYCAM_STRONG_GROUPS)
        hydroxyl_h = len(names & _GLYCAM_HYDROXY_H) >= 2
        code_groups.append((code, lowercase, strong_atoms or hydroxyl_h))
    strong = sorted({code for code, lowercase, atom_signature in code_groups if lowercase or atom_signature})
    if strong:
        return strong
    distinct = sorted({code for code, _, _ in code_groups})
    if len(distinct) >= 2 and len(code_groups) >= 2:
        return distinct
    return []


def _amber_signatures(atoms: Sequence[AtomRecord]) -> List[str]:
    resnames = {a.resname_orig.strip().upper() for a in atoms}
    hits = sorted((resnames & _AMBER_TERMINAL) | (resnames & _AMBER_NUCLEIC))
    return hits


def _specific_auto_candidates(atoms: Sequence[AtomRecord]) -> List[Tuple[str, str]]:
    resnames = {a.resname_orig.strip().upper() for a in atoms}
    candidates: List[Tuple[str, str]] = []
    glycam = _glycam_signatures(atoms)
    if glycam:
        candidates.append(("glycam", "multiple/strong GLYCAM residue and atom signature(s): " + ",".join(glycam[:8])))
    gmx_charmm = sorted(resnames & _GMX_CHARMM_STRONG)
    if gmx_charmm:
        candidates.append(("gromacs-charmm", "GROMACS/CHARMM residue signature(s): " + ",".join(gmx_charmm)))
    opls_unique = sorted(resnames & _GMX_OPLS_STRONG)
    opls_all_atom = _opls_all_atom_signatures(atoms)
    if opls_unique:
        candidates.append(("gromacs-oplsaa", "GROMACS/OPLS-AA residue signature(s): " + ",".join(opls_unique)))
    elif opls_all_atom:
        candidates.append(("gromacs-oplsaa", "OPLS/GROMOS-shared residue naming with explicit aliphatic hydrogens: " + "; ".join(opls_all_atom[:4])))
    amber = _amber_signatures(atoms)
    if amber:
        gmx_atoms = _gromacs_atom_signatures(atoms)
        if gmx_atoms:
            candidates.append(("gromacs-amber", "AMBER-family residue signature(s) plus GROMACS atom naming: " + ",".join(amber[:8]) + "; " + "; ".join(gmx_atoms[:4])))
        else:
            candidates.append(("amber", "AMBER-family residue signature(s): " + ",".join(amber[:8])))
    return candidates


def resolve_input_convention(
    requested: str,
    atoms: Sequence[AtomRecord],
    preflight: PdbFormatPreflight,
    legacy_charmm: Optional[bool] = None,
) -> Tuple[str, List[str]]:
    if legacy_charmm is not None:
        resolved = "charmm" if legacy_charmm else "generic"
        return resolved, [f"Legacy format override resolved input convention as {resolved}."]
    explicit = {
        "charmm": "CHARMM/CHARMM-GUI",
        "amber": "AMBER/LEaP",
        "gromacs-amber": "GROMACS using an AMBER-family force field",
        "gromacs-charmm": "GROMACS using a CHARMM-family force field",
        "gromacs-oplsaa": "GROMACS OPLS-AA/L all-atom nomenclature",
        "gromacs-gromos": "a recovered GROMACS GROMOS-family nomenclature",
        "glycam": "GLYCAM 06j-1/AMBER carbohydrate nomenclature",
        "generic": "generic",
    }
    if requested in explicit:
        return requested, [f"Input convention explicitly set to {explicit[requested]}."]
    candidates = _specific_auto_candidates(atoms)
    names = sorted({name for name, _ in candidates})
    if len(names) > 1:
        details = "; ".join(f"{name}: {evidence}" for name, evidence in candidates)
        raise SystemExit(
            "Auto-detection found conflicting source-nomenclature signatures. "
            f"Candidates: {', '.join(names)}. Evidence: {details}. "
            "Set --input-convention explicitly."
        )
    if len(names) == 1:
        name = names[0]
        evidence = next(item for candidate, item in candidates if candidate == name)
        return name, [f"Auto-detected {name} input from {evidence}."]
    resnames = {a.resname_orig.strip().upper() for a in atoms}
    shared = sorted(resnames & (_GMX_OPLS_GROMOS_SHARED | {"CYS2"}))
    gmx_atoms = _gromacs_atom_signatures(atoms)
    if shared or len(gmx_atoms) >= 2:
        evidence = []
        if shared:
            evidence.append("shared GROMACS force-field residue state(s): " + ",".join(shared))
        if len(gmx_atoms) >= 2:
            evidence.append("multiple GROMACS-style atom signatures: " + "; ".join(gmx_atoms))
        raise SystemExit(
            "Auto-detection recognized GROMACS-style nomenclature but could not determine the force-field family unambiguously. "
            + "Evidence: " + "; ".join(evidence) + ". "
            "Set --input-convention to gromacs-amber, gromacs-charmm, gromacs-oplsaa or gromacs-gromos."
        )
    detected, evidence = detect_charmm(atoms, preflight)
    if detected:
        return "charmm", evidence
    return "generic", ["No unambiguous force-field-specific signature detected; using generic input convention."]


def build_report(
    requested: str,
    resolved: str,
    evidence: List[str],
    converted: Sequence[AtomRecord],
    readiness_issues: List[str],
    warnings: List[str],
    has_nucleic: bool,
) -> CompatibilityReport:
    report = CompatibilityReport(
        requested_convention=requested,
        resolved_convention=resolved,
        convention_evidence=list(evidence),
        readiness_issues=list(readiness_issues),
        warnings=list(warnings),
        has_nucleic=has_nucleic,
    )
    n_h = sum(1 for a in converted if a.element.upper() == "H")
    if converted and n_h == 0:
        msg = "Warning: Structure contains no explicit hydrogens."
        if msg not in report.warnings:
            report.warnings.append(msg)
        if has_nucleic:
            nuc = "Warning: Hydrogen-free nucleic-acid ONEBEAD has not been quantitatively validated."
            if nuc not in report.warnings:
                report.warnings.append(nuc)
    return report
