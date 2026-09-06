from __future__ import annotations
import re
from typing import Dict, List, Sequence, Tuple
from .models import AtomRecord, CompatibilityIssue
from .constants import MAX_REPORTED_ISSUE_GROUPS
from .pdbio import group_residues, derive_chain_from_segid

def stop_for_compatibility_issues(issues: Sequence[CompatibilityIssue]) -> None:
    if not issues:
        return

    grouped = {}
    for issue in issues:
        if issue.key not in grouped:
            grouped[issue.key] = [issue, 0]
        grouped[issue.key][1] += 1

    lines = [
        f"ERROR: found {len(issues)} incompatible residue/moiety issue(s) "
        f"in {len(grouped)} group(s)."
    ]
    shown = list(grouped.values())[:MAX_REPORTED_ISSUE_GROUPS]
    for issue, count in shown:
        label = issue.key.split(":", 1)[0]
        repeat = f" ({count} occurrences; first shown)" if count > 1 else ""
        lines.append(f"- [{label}] {issue.message}{repeat}")

    if len(grouped) > len(shown):
        lines.append(
            f"- ... and {len(grouped) - len(shown)} more issue group(s)."
        )

    hints = []
    for issue in issues:
        if issue.hint and issue.hint not in hints:
            hints.append(issue.hint)
    if hints:
        lines.append("")
        lines.extend(f"HINT: {hint}" for hint in hints)

    lines.append("Conversion stopped before writing the output PDB.")
    raise SystemExit("\n".join(lines))


def partial_selection_issues(
    all_atoms: Sequence[AtomRecord], selected_indices: Sequence[int]
) -> List[str]:
    selected = set(selected_indices)
    issues: List[str] = []
    for residue in group_residues(all_atoms):
        kept = sum(atom.input_atom_index in selected for atom in residue)
        if kept in {0, len(residue)}:
            continue
        first = residue[0]
        chain = first.chain_orig or derive_chain_from_segid(first.segid) or "?"
        omitted = len(residue) - kept
        issues.append(
            f"{first.resname_orig} at chain {chain}, input residue "
            f"{first.resseq_orig}: selected {kept}/{len(residue)} atoms "
            f"and left {omitted} atom(s) outside the requested range."
        )
    if len(issues) > MAX_REPORTED_ISSUE_GROUPS:
        omitted = len(issues) - MAX_REPORTED_ISSUE_GROUPS
        issues = issues[:MAX_REPORTED_ISSUE_GROUPS]
        issues.append(f"... and {omitted} more incomplete selected residues.")
    return issues


def selection_boundary_notes(
    all_atoms: Sequence[AtomRecord], selected_indices: Sequence[int]
) -> Tuple[List[str], bool]:
    if not selected_indices:
        return [], False

    def describe(idx: int) -> str:
        atom = all_atoms[idx - 1]
        chain = atom.chain_orig or derive_chain_from_segid(atom.segid) or "?"
        return (
            f"atom-order {idx} -> PDB serial {atom.input_serial}: "
            f"{atom.atom_name} {atom.resname_orig} chain {chain} "
            f"residue {atom.resseq_orig}"
        )

    first_idx = selected_indices[0]
    last_idx = selected_indices[-1]
    notes = [
        f"Selection start: {describe(first_idx)}.",
        f"Selection end: {describe(last_idx)}.",
    ]
    mismatch = (
        all_atoms[first_idx - 1].input_serial != first_idx
        or all_atoms[last_idx - 1].input_serial != last_idx
    )
    if mismatch:
        notes.append(
            "Selection numbering note: -a/--atoms uses 1-based ATOM/HETATM "
            "record order after MODEL/altLoc/solvent filtering, not the PDB "
            "serial field. TER and other non-coordinate records are not counted."
        )
    return notes, mismatch


def parse_range(expr: str, n_atoms: int) -> List[int]:
    expr = (expr or "all").strip().lower()
    if expr in {"all", "*"}:
        return list(range(1, n_atoms + 1))
    selected: List[int] = []
    seen = set()
    for part in expr.split(','):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise SystemExit(f"Invalid atom range component: {part!r}")
        a = int(m.group(1))
        b = int(m.group(2)) if m.group(2) else a
        if a < 1 or b < 1 or b < a:
            raise SystemExit(f"Invalid atom range component: {part!r}")
        if b > n_atoms:
            raise SystemExit(f"Atom range {part!r} exceeds number of ATOM/HETATM records ({n_atoms}).")
        for idx in range(a, b + 1):
            if idx not in seen:
                selected.append(idx)
                seen.add(idx)
    return selected


def parse_serial_selection(expr: str, all_atoms: Sequence[AtomRecord]) -> List[int]:
    expression = (expr or "").strip().lower()
    if expression in {"", "all", "*"}:
        return [atom.input_atom_index for atom in all_atoms]

    intervals: List[Tuple[int, int]] = []
    for part in expression.split(','):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise SystemExit(f"Invalid PDB serial range component: {part!r}")
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        if start < 1 or end < 1 or end < start:
            raise SystemExit(f"Invalid PDB serial range component: {part!r}")
        intervals.append((start, end))
    if not intervals:
        return []

    def requested(serial: int) -> bool:
        return any(start <= serial <= end for start, end in intervals)

    selected_atoms = [atom for atom in all_atoms if requested(atom.input_serial)]
    if not selected_atoms:
        raise SystemExit(
            f"PDB serial selection {expr!r} did not match any ATOM/HETATM records."
        )

    by_serial: Dict[int, List[int]] = {}
    for atom in selected_atoms:
        by_serial.setdefault(atom.input_serial, []).append(atom.input_atom_index)
    duplicated = {serial: idxs for serial, idxs in by_serial.items() if len(idxs) > 1}
    if duplicated:
        preview = sorted(duplicated.items())[:5]
        detail = "; ".join(
            f"serial {serial} occurs at atom-order {','.join(map(str, idxs))}"
            for serial, idxs in preview
        )
        if len(duplicated) > len(preview):
            detail += f"; ... and {len(duplicated) - len(preview)} more repeated serial(s)"
        raise SystemExit(
            "PDB serial selection is ambiguous because the selected ATOM/HETATM records "
            f"contain repeated atom serials ({detail}). Use -a/--atoms record-order "
            "selection for files with repeated or wrapped PDB serials."
        )

    return [atom.input_atom_index for atom in selected_atoms]
