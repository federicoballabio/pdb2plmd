from __future__ import annotations
from collections import Counter
from dataclasses import replace
from typing import List, Tuple
from .models import AtomRecord, CompatibilityIssue
from .constants import (
    CHAIN_IDS, PROTEIN_NAMES, ION_MAP, ION_ELEMENT, UNSUPPORTED_ION_NAMES,
    GLYCAN_UNSUPPORTED, GLYCOSYLATION_SITE_MAP, UNSUPPORTED_SITE_NAMES,
    TERMINAL_CAP_NAMES, SOLVENT_AND_ADDITIVES,
)
from .glycan import canonical_glycan, normalize_glycan_atom_name, validate_glycan_heavy_atoms
from .protein import resolve_protein_state, protein_missing_heavy
from .nucleic import resolve_nucleic_base, source_terminal_suffix, normalize_nucleic_atom_names, validate_nucleic
from .pdbio import (
    group_residues, split_into_chains, derive_chain_from_segid, next_chain_id,
    base_resname, residue_atom_names, existing_terminal_suffix, terminal_suffix,
    normalize_atom_names_for_residue, normalize_atom_name,
)
from .selection import stop_for_compatibility_issues
from .conventions import (
    resolve_amber_residue, resolve_gromacs_residue, normalize_gromacs_protein_atom_names,
    resolve_glycam_residue, normalize_glycam_protein_atom_names,
)


def convert_atoms(
    atoms: List[AtomRecord],
    split_on_gaps: bool = False,
    input_convention: str = "generic",
) -> Tuple[List[AtomRecord], List[str], List[str], List[str], bool]:
    if not atoms:
        raise SystemExit("The atom selection is empty.")

    log: List[str] = []
    errors: List[CompatibilityIssue] = []
    readiness_issues: List[str] = []
    warnings: List[str] = []
    has_nucleic = False
    residues = group_residues(atoms)
    chains = split_into_chains(residues, split_on_gaps=split_on_gaps)
    used_chains = set()
    converted: List[AtomRecord] = []
    log.append(f"Input selected atoms: {len(atoms)}")
    log.append(f"Input selected residues: {len(residues)}")
    log.append(f"Output chains: {len(chains)}")
    if input_convention == "gromacs-gromos":
        readiness_issues.append(
            "GROMOS-family source uses a united-atom representation with implicit nonpolar hydrogens; the resulting mixed explicit/implicit-hydrogen bead centers have not been quantitatively validated for SAXS ONEBEAD"
        )
        log.append("GROMOS readiness policy: recognized/mappable structures are not ONEBEAD_SAXS_READY until united-atom quantitative validation is completed.")
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
        log.append(f"Chain {chain_idx}: ID {chain_id}; residues={len(chain_residues)}; atoms={atom_count}.")
        if len(chain_residues) > 9999:
            errors.append(CompatibilityIssue(
                "pdb-limit:residues",
                f"Chain {chain_id} has more than 9999 residues and cannot fit the PDB residue-number field.",
            ))

        for residue_index, res_atoms in enumerate(chain_residues, start=1):
            orig = res_atoms[0]
            raw = orig.resname_orig.strip().upper()
            names = residue_atom_names(res_atoms)
            is_first = residue_index == 1
            is_last = residue_index == len(chain_residues)

            source_raw = orig.resname_orig.strip() if input_convention == "glycam" else raw
            convention_notes: List[str] = []
            convention_error = None
            glycam_state = None
            if input_convention == "amber":
                raw, convention_error, convention_notes = resolve_amber_residue(raw, is_first, is_last)
            elif input_convention.startswith("gromacs-"):
                raw, convention_error, convention_notes = resolve_gromacs_residue(raw, input_convention, is_first, is_last, sorted(names))
            elif input_convention == "glycam":
                raw, convention_error, convention_notes, glycam_state = resolve_glycam_residue(orig.resname_orig.strip())
            if convention_error is not None:
                errors.append(convention_error)
                continue
            if raw != source_raw:
                log.append(f"Source residue: {source_raw}{orig.resseq_orig} -> {raw} ({input_convention}).")
            log.extend(convention_notes)
            base = base_resname(raw)

            if raw in UNSUPPORTED_SITE_NAMES or base in UNSUPPORTED_SITE_NAMES:
                name = raw if raw in UNSUPPORTED_SITE_NAMES else base
                errors.append(CompatibilityIssue(
                    f"unsupported-site:{name}",
                    f"Glycosylation-site residue {name} at chain {chain_id}, input residue {orig.resseq_orig} is not compatible with ONEBEAD. Hydroxyproline contains OD1, which has no PRO LCPO mapping.",
                    "Do not approximate OLP as PRO; a dedicated validated ONEBEAD/LCPO mapping is required.",
                ))
                continue

            if raw in TERMINAL_CAP_NAMES or base in TERMINAL_CAP_NAMES:
                name = raw if raw in TERMINAL_CAP_NAMES else base
                errors.append(CompatibilityIssue(
                    f"terminal-cap:{name}",
                    f"Terminal capping group {name} at chain {chain_id}, input residue {orig.resseq_orig} is not a monosaccharide and has no ONEBEAD bead.",
                    "Remove GLYCAM terminal caps before conversion.",
                ))
                continue

            if raw in GLYCOSYLATION_SITE_MAP or base in GLYCOSYLATION_SITE_MAP:
                key = raw if raw in GLYCOSYLATION_SITE_MAP else base
                parent = GLYCOSYLATION_SITE_MAP[key]
                log.append(f"Glycosylation site: {key}{orig.resseq_orig} -> {parent}.")
                base = parent
                raw = parent

            glycan = canonical_glycan(raw, names)
            if glycan is None and raw != base:
                glycan = canonical_glycan(base, names)
            ion = ION_MAP.get(raw) or ION_MAP.get(base)
            nuc_base, nuc_error, nuc_notes = resolve_nucleic_base(orig.resname_orig, names, input_convention)
            if nuc_error is not None:
                errors.append(nuc_error)
                continue
            log.extend(nuc_notes)

            suffix = ""
            resname_out = ""
            output_names: List[str]

            if glycan is not None:
                unexpected_elements = sorted({a.element.upper() for a in res_atoms} - {"H","C","N","O"})
                if unexpected_elements:
                    errors.append(CompatibilityIssue(
                        f"modified-glycan:{glycan}:{','.join(unexpected_elements)}",
                        f"Supported glycan {orig.resname_orig} at chain {chain_id}, input residue {orig.resseq_orig} contains unexpected element(s) {','.join(unexpected_elements)}. This indicates a substituted moiety not covered by the implemented monosaccharide bead.",
                        "Remove the unsupported substituent or implement and validate a dedicated ONEBEAD parameterization.",
                    ))
                    continue
                resname_out = glycan
                if glycan != raw:
                    log.append(f"Glycan: {orig.resname_orig}{orig.resseq_orig} -> {glycan}.")
                fatal_g, ready_g, glycan_notes = validate_glycan_heavy_atoms(res_atoms, resname_out, input_convention)
                for issue_key, issue_text in fatal_g:
                    errors.append(CompatibilityIssue(
                        f"{issue_key}:{resname_out}",
                        f"Supported glycan {orig.resname_orig} at chain {chain_id}, input residue {orig.resseq_orig}: {issue_text}.",
                        "Do not silently repair or approximate substituted glycan chemistry; provide a supported monosaccharide.",
                    ))
                if fatal_g:
                    continue
                for issue_key, issue_text in ready_g:
                    readiness_issues.append(
                        f"{resname_out} chain {chain_id} residue {orig.resseq_orig}: {issue_text}."
                    )
                for note in glycan_notes:
                    log.append(f"Glycan state: {orig.resname_orig}{orig.resseq_orig}: {note}")
                output_names = [normalize_glycan_atom_name(a.atom_name, resname_out, input_convention) for a in res_atoms]

            elif raw in GLYCAN_UNSUPPORTED or base in GLYCAN_UNSUPPORTED:
                name = raw if raw in GLYCAN_UNSUPPORTED else base
                errors.append(CompatibilityIssue(
                    f"unsupported-glycan:{name}",
                    f"Monosaccharide {name} at chain {chain_id}, input residue {orig.resseq_orig} has no ONEBEAD parameters.",
                    "No chemically similar monosaccharide is substituted automatically.",
                ))
                continue

            elif ion is not None:
                if len(res_atoms) != 1:
                    errors.append(CompatibilityIssue(
                        f"ion-atom-count:{raw}",
                        f"Ion residue {raw} at chain {chain_id}, input residue {orig.resseq_orig} contains {len(res_atoms)} atoms; ONEBEAD expects exactly one atom.",
                    ))
                    continue
                resname_out = ion
                if ion != raw:
                    log.append(f"Ion: {orig.resname_orig}{orig.resseq_orig} -> {ion}.")
                output_names = normalize_atom_names_for_residue(res_atoms, resname_out)

            elif raw in UNSUPPORTED_ION_NAMES or base in UNSUPPORTED_ION_NAMES:
                name = raw if raw in UNSUPPORTED_ION_NAMES else base
                errors.append(CompatibilityIssue(
                    f"unsupported-ion:{name}",
                    f"Recognized monoatomic ion {orig.resname_orig} ({UNSUPPORTED_ION_NAMES[name]}) at chain {chain_id}, input residue {orig.resseq_orig} has no parameter in the frozen SAXS ONEBEAD ion set.",
                    "Do not substitute a different ion automatically; remove it from the SAXS selection or add and validate a dedicated SAXS ONEBEAD ion parameter.",
                ))
                continue

            elif nuc_base is not None:
                has_nucleic = True
                suffix = (existing_terminal_suffix(orig.resname_orig) or
                          source_terminal_suffix(orig.resname_orig, input_convention) or
                          terminal_suffix(nuc_base, names, is_first, is_last))
                resname_out = nuc_base + suffix
                output_names = normalize_nucleic_atom_names(res_atoms, resname_out, input_convention)
                fatal_n, ready_n = validate_nucleic(res_atoms, nuc_base, suffix, output_names)
                if fatal_n:
                    errors.extend(fatal_n)
                    continue
                for issue in ready_n:
                    readiness_issues.append(
                        f"{resname_out} chain {chain_id} residue {orig.resseq_orig}: {issue}."
                    )
                if suffix:
                    log.append(f"Terminal residue: {orig.resname_orig}{orig.resseq_orig} -> {resname_out}.")

            elif base in PROTEIN_NAMES or raw in PROTEIN_NAMES:
                protein_name = raw if raw in PROTEIN_NAMES else base
                resname_out, chemistry, protein_warnings = resolve_protein_state(protein_name)
                for msg in protein_warnings:
                    if msg not in warnings:
                        warnings.append(msg)
                    log.append(msg)
                if resname_out != protein_name:
                    log.append(f"Protein state: {orig.resname_orig}{orig.resseq_orig} -> {resname_out}.")
                missing = protein_missing_heavy(res_atoms, chemistry, input_convention)
                if missing:
                    readiness_issues.append(
                        f"{orig.resname_orig} chain {chain_id} residue {orig.resseq_orig}: missing required heavy atom(s): {','.join(missing)}."
                    )
                if input_convention.startswith("gromacs-"):
                    output_names = normalize_gromacs_protein_atom_names(res_atoms, resname_out, is_first, is_last)
                elif input_convention == "glycam":
                    output_names = normalize_glycam_protein_atom_names(res_atoms, resname_out)
                else:
                    output_names = normalize_atom_names_for_residue(res_atoms, resname_out)

            elif raw in SOLVENT_AND_ADDITIVES or base in SOLVENT_AND_ADDITIVES:
                name = raw if raw in SOLVENT_AND_ADDITIVES else base
                errors.append(CompatibilityIssue(
                    f"unsupported-solvent:{name}",
                    f"Solvent or crystallisation additive {name} at chain {chain_id}, input residue {orig.resseq_orig} has no ONEBEAD mapping.",
                    "Use --drop-solvent to remove water and common crystallisation additives.",
                ))
                continue

            else:
                errors.append(CompatibilityIssue(
                    f"unsupported-residue:{raw}",
                    f"Residue or moiety {orig.resname_orig} at chain {chain_id}, input residue {orig.resseq_orig} has no implemented ONEBEAD mapping.",
                    "Remove it from the PLUMED ATOMS/TEMPLATE selection or add a validated mapping to both pdb2plmd and SAXS.cpp.",
                ))
                continue

            duplicates = sorted(name for name, count in Counter(output_names).items() if count > 1)
            if duplicates:
                errors.append(CompatibilityIssue(
                    f"duplicate-atom-names:{resname_out}:{','.join(duplicates)}",
                    f"Duplicate output atom names in {resname_out}{residue_index} chain {chain_id}: {','.join(duplicates)}.",
                    "Check the input naming and residue-specific atom-name conversion before using the template.",
                ))
                continue

            for atom, atom_out in zip(res_atoms, output_names):
                if len(atom_out) > 4:
                    errors.append(CompatibilityIssue(
                        f"atom-name-too-long:{atom_out}",
                        f"Atom name {atom_out!r} exceeds four PDB columns at input line {atom.line_number}.",
                    ))
                    continue
                if atom_out != normalize_atom_name(atom.atom_name):
                    log.append(f"Atom name: {normalize_atom_name(atom.atom_name)} -> {atom_out} in {resname_out}{residue_index}.")
                element_out = atom.element
                if ion is not None:
                    expected_element = ION_ELEMENT[resname_out]
                    if element_out.upper() != expected_element:
                        log.append(f"Ion element: {element_out or '?'} -> {expected_element} in {resname_out}{residue_index}.")
                    element_out = expected_element
                converted.append(replace(
                    atom,
                    chain_out=chain_id,
                    resseq_out=residue_index,
                    resname_out=resname_out,
                    atom_name_out=atom_out,
                    element=element_out,
                    source_convention=("glycam" if glycam_state is not None else ""),
                    chemical_identity=(glycam_state.identity if glycam_state is not None else ""),
                    chemical_state=(
                        f"anomer={glycam_state.anomer};config={glycam_state.config};ring={glycam_state.ring};linkage={glycam_state.linkage_code}"
                        if glycam_state is not None else ""
                    ),
                    glycan_anomer=(glycam_state.anomer if glycam_state is not None else ""),
                    glycan_config=(glycam_state.config if glycam_state is not None else ""),
                    glycan_ring=(glycam_state.ring if glycam_state is not None else ""),
                    glycan_linkage=(glycam_state.linkage_code if glycam_state is not None else ""),
                    glycan_ccd=(glycam_state.ccd_comp_id if glycam_state is not None else ""),
                    mapping_provenance=("GLYCAM_CODE_INDEX_C4.tsv" if glycam_state is not None else ""),
                ))

    stop_for_compatibility_issues(errors)
    if len(converted) > 99999:
        raise SystemExit("The output contains more than 99999 atoms and cannot fit PDB columns.")
    return converted, log, readiness_issues, warnings, has_nucleic
