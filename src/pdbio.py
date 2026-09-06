from __future__ import annotations
import re
from dataclasses import replace
from typing import List, Optional, Sequence, Tuple
from .models import AtomRecord, PdbFormatPreflight
from .constants import (
    CHAIN_IDS, RNA_MAP, DNA_MAP, ION_MAP, PHOSPHATE_ATOMS, SOLVENT_AND_ADDITIVES,
)

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

    if name[0].isdigit() and len(name) > 1:
        name = name[1:]
    if resname.strip().upper() in ION_MAP:
        ion = ION_MAP[resname.strip().upper()]
        return "CA" if ion == "CAL" else ion[:2]

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


def preflight_pdb_format(path: str, sample_limit: int = 5000) -> PdbFormatPreflight:
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


def _read_resname(
    padded: str,
    charmm_override: Optional[bool],
    preserve_case: bool = False,
) -> str:
    if charmm_override is True:
        field = padded[17:21]
    elif charmm_override is False:
        field = padded[17:20]
    elif padded[20:21].strip() and not padded[21:22].strip():
        field = padded[17:21]
    else:
        field = padded[17:20]
    value = field.strip()
    return value if preserve_case else value.upper()


def normalize_altloc_request(value: str) -> str:
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
    return min(entries, key=lambda entry: (-entry[1].occ, entry[0]))


def _select_altloc_for_residue(
    entries: Sequence[Tuple[int, AtomRecord]],
    requested: str,
) -> Tuple[List[Tuple[int, AtomRecord]], Optional[str], List[str]]:
    by_name = {}
    labels = set()
    for order, atom in entries:
        by_name.setdefault(atom.atom_name, []).append((order, atom))
        if atom.altloc:
            labels.add(atom.altloc)

    if not labels:



        return sorted(entries), None, []

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
    model: Optional[int] = None,
    drop_solvent: bool = False,
    charmm: Optional[bool] = None,
    altloc: str = "auto",
    preserve_resname_case: bool = False,
) -> Tuple[List[AtomRecord], List[str]]:
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
                if model is None and not target_model_seen:


                    model = current_model
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
            resn = _read_resname(padded, charmm, preserve_case=preserve_resname_case)
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
    if model_records_seen:
        if n_model_skipped:
            notes.append(
                f"MODEL filter: kept MODEL {model}; dropped {n_model_skipped} atom records."
            )
        else:
            notes.append(f"MODEL filter: kept MODEL {model}.")
    if n_solvent_skipped:
        notes.append(f"Solvent filter: dropped {n_solvent_skipped} atom records.")




    residue_groups: List[List[Tuple[int, AtomRecord]]] = []
    current_entries: List[Tuple[int, AtomRecord]] = []
    last_key = None
    for order, atom in enumerate(provisional):
        key = (
            atom.chain_orig,
            atom.segid,
            atom.resseq_orig,
            atom.icode_orig,
            atom.resname_orig,
        )
        if current_entries and (atom.ter_before or key != last_key):
            residue_groups.append(current_entries)
            current_entries = []
        current_entries.append((order, atom))
        last_key = key
    if current_entries:
        residue_groups.append(current_entries)

    selected_entries: List[Tuple[int, AtomRecord]] = []
    altloc_selections = []
    for entries in residue_groups:
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

    name = name.strip().replace('*', "'")

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
        if current and (a.ter_before or k != last_key):
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

    if is_first_in_chain and ("OP3" in atom_names or "O3P" in atom_names or "HOP3" in atom_names or "HP" in atom_names):
        return "T"

    if is_first_in_chain and not (atom_names & PHOSPHATE_ATOMS):
        if {"H5T", "HO5'"} & atom_names or "O5'" in atom_names:
            return "5"

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
