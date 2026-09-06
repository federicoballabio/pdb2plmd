from __future__ import annotations
import argparse
import math
import sys
import traceback
from pathlib import Path
from typing import List
from ._version import VERSION
from .models import AtomRecord, CompatibilityReport
from .pdbio import preflight_pdb_format, parse_pdb
from .selection import parse_serial_selection, parse_range, selection_boundary_notes, partial_selection_issues
from .compatibility import resolve_input_convention, build_report
from .convert import convert_atoms
from .writer import default_log_path, write_success_log, atomic_write_pdb, write_error_log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare and validate an ordered PDB template for PLUMED SAXS ... ONEBEAD.")
    p.add_argument("--version", action="version", version=f"pdb2plmd {VERSION}")
    p.add_argument("-i", "--input", required=True, help="Input PDB extracted from the simulation/TPR selection.")
    p.add_argument("-o", "--output", required=True, help="Output SAXS ONEBEAD-compatible target PDB.")
    selection = p.add_mutually_exclusive_group()
    selection.add_argument("-a", "--atoms", default=None, help="1-based ATOM/HETATM record-order range to keep, e.g. '1-1062' or '1-100,150,200-250'. TER records are not counted. Default: all when neither -a nor -s is given.")
    selection.add_argument("-s", "--serials", default=None, help="PDB atom-serial range to keep, e.g. '1-1069'. Use -a for files with repeated/wrapped atom serials.")
    p.add_argument("--model", type=int, default=None, help="MODEL serial to keep from a multi-model file. Default: first MODEL encountered.")
    p.add_argument("--input-convention", choices=("auto","generic","charmm","amber","gromacs-amber","gromacs-charmm","gromacs-oplsaa","gromacs-gromos","glycam"), default="auto", help="Source nomenclature convention. auto uses conservative signature-based detection and stops on recognized ambiguity. Default: auto.")
    p.add_argument("--charmm", dest="charmm", action="store_true", default=None, help=argparse.SUPPRESS)
    p.add_argument("--no-charmm", dest="charmm", action="store_false", help=argparse.SUPPRESS)
    p.add_argument("--split-on-gaps", action="store_true", help="Start a new chain at residue-number resets or gaps.")
    p.add_argument("-altloc", "--altloc", default="auto", metavar="LABEL", help="Alternate-location handling. Default: residue-level automatic selection; give A/B/etc. to force one label.")
    p.add_argument("--drop-solvent", action="store_true", help="Remove water and common crystallisation additives.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose mode. Print detailed terminal status and write <output-stem>.log on success. On failure a log is written in both standard and verbose modes.")
    p.add_argument("-g", "--log", nargs="?", const="__AUTO__", default=None, help=argparse.SUPPRESS)
    return p.parse_args()


def run_conversion(args: argparse.Namespace) -> tuple[List[AtomRecord], List[int], List[AtomRecord], CompatibilityReport, List[str]]:
    preflight = preflight_pdb_format(args.input)
    parse_charmm = args.charmm
    if parse_charmm is None:
        if args.input_convention == "charmm":
            parse_charmm = True
        elif args.input_convention == "generic":
            parse_charmm = False
    context_lines: List[str] = [
        f"Input format preflight: {preflight.style}.",
        f"Preflight sample: {preflight.sampled_atoms} atom records; chain column populated in {preflight.chain_records}; SEGID populated in {preflight.segid_records}; four-character residue extension observed in {preflight.four_char_residue_records}.",
        "Residue-name parsing: record-aware 3-/4-character PDB/CHARMM handling.",
        "Chain handling: standard PDB chain ID when present; otherwise SEGID-derived fallback.",
    ]
    if args.charmm is not None:
        context_lines.append("Legacy format override: " + ("CHARMM forced." if args.charmm else "CHARMM handling disabled."))

    all_atoms, parse_notes = parse_pdb(
        args.input,
        model=args.model,
        charmm=parse_charmm,
        drop_solvent=args.drop_solvent,
        altloc=args.altloc,
        preserve_resname_case=(args.input_convention in {"glycam", "auto"}),
    )
    context_lines.extend(parse_notes)
    if not all_atoms:
        raise SystemExit(f"No ATOM/HETATM records found in {args.input}")

    resolved_convention, evidence = resolve_input_convention(args.input_convention, all_atoms, preflight, args.charmm)
    context_lines.extend(evidence)
    if resolved_convention == "glycam":
        context_lines.append("GLYCAM residue-name parsing: source case preserved because GLYCAM residue codes are case-sensitive.")

    if args.serials is not None:
        selected_indices = parse_serial_selection(args.serials, all_atoms)
    else:
        selected_indices = parse_range(args.atoms if args.atoms is not None else "all", len(all_atoms))
    selected_indices = sorted(set(selected_indices))
    index_set = set(selected_indices)
    selected_atoms = [a for a in all_atoms if a.input_atom_index in index_set]
    if not selected_atoms:
        raise SystemExit("The atom selection is empty.")

    boundary_notes, serial_mismatch = selection_boundary_notes(all_atoms, selected_indices)
    context_lines.extend(boundary_notes)
    if args.serials is None and serial_mismatch and (args.atoms or "all").strip().lower() not in {"all", "*"}:
        print("WARNING: -a/--atoms uses ATOM/HETATM record order, not PDB serials; " + f"the selected end is {boundary_notes[1].removeprefix('Selection end: ').rstrip('.')}.", file=sys.stderr)

    selection_issues = partial_selection_issues(all_atoms, selected_indices)
    if selection_issues:
        lines = ["ERROR: atom selection leaves one or more incomplete residues.", *(f"- {issue}" for issue in selection_issues), "Select complete amino-acid/nucleic-acid residues (and complete multi-atom moieties) only."]
        if boundary_notes:
            lines.append("Selection boundary:")
            lines.extend(f"- {note}" for note in boundary_notes)
        raise SystemExit("\n".join(lines))

    try:
        converted, conversion_lines, readiness_issues, conversion_warnings, has_nucleic = convert_atoms(
            selected_atoms,
            split_on_gaps=args.split_on_gaps,
            input_convention=resolved_convention,
        )
    except SystemExit as exc:
        error_text = str(exc) if exc.code not in (None, "") else "Conversion failed."
        if boundary_notes:
            error_text += "\n\nSelection boundary:\n" + "\n".join(f"- {note}" for note in boundary_notes)
        raise SystemExit(error_text) from None
    context_lines.extend(conversion_lines)

    if any(not math.isfinite(value) for atom in converted for value in (atom.x, atom.y, atom.z, atom.occ, atom.bfac)):
        raise SystemExit("The selected atoms contain non-finite numeric values.")

    report = build_report(
        args.input_convention,
        resolved_convention,
        evidence,
        converted,
        readiness_issues,
        conversion_warnings,
        has_nucleic,
    )
    for warning in report.warnings:
        print(warning, file=sys.stderr)
    for issue in report.readiness_issues:
        print(f"Warning: {issue}", file=sys.stderr)
    return all_atoms, selected_indices, converted, report, context_lines


def main() -> int:
    args = parse_args()
    if args.log is not None:
        args.verbose = True
    auto_log = default_log_path(args.output)
    log_path = auto_log if args.log in (None, "__AUTO__") else args.log
    output_existed_before = Path(args.output).exists()
    context_lines: List[str] = []
    try:
        all_atoms, selected_indices, converted, report, context_lines = run_conversion(args)
        if args.verbose:
            write_success_log(log_path, args, all_atoms, selected_indices, converted, report, context_lines)
        atomic_write_pdb(converted, args.output)
        if args.verbose:
            print(
                f"pdb2plmd {VERSION}: OK - wrote {args.output} ({len(converted)} atoms); "
                f"ONEBEAD_SAXS_READY={report.onebead_saxs_ready}; log {log_path}"
            )
        else:
            print(f"pdb2plmd {VERSION}: PDB converted.")
        return 0
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
        if code == 0:
            return 0
        error_text = str(exc) if exc.code not in (None, "") else "Conversion failed."
        try:
            actual_log = write_error_log(log_path, args, error_text, context_lines, output_existed_before=output_existed_before)
            if args.verbose:
                first_error = error_text.splitlines()[0]
                if first_error.upper().startswith("ERROR:"):
                    first_error = first_error.split(":", 1)[1].strip()
                print(f"pdb2plmd {VERSION}: ERROR - {first_error}", file=sys.stderr)
                print(f"Details: {actual_log}", file=sys.stderr)
            else:
                print(f"pdb2plmd {VERSION}: PDB conversion failed.", file=sys.stderr)
        except OSError as log_exc:
            print(f"pdb2plmd {VERSION}: PDB conversion failed.", file=sys.stderr)
            print(f"Error log could not be written: {log_exc}", file=sys.stderr)
        return code or 1
    except Exception as exc:
        tb = traceback.format_exc()
        error_text = f"Unexpected {type(exc).__name__}: {exc}"
        try:
            actual_log = write_error_log(log_path, args, error_text, context_lines, unexpected_traceback=tb, output_existed_before=output_existed_before)
            if args.verbose:
                print(f"pdb2plmd {VERSION}: ERROR - {error_text}", file=sys.stderr)
                print(f"Details: {actual_log}", file=sys.stderr)
            else:
                print(f"pdb2plmd {VERSION}: PDB conversion failed.", file=sys.stderr)
        except OSError as log_exc:
            print(f"pdb2plmd {VERSION}: PDB conversion failed.", file=sys.stderr)
            print(f"Error log could not be written: {log_exc}", file=sys.stderr)
        return 1
