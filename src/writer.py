from __future__ import annotations
import argparse
import os
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence
from ._version import VERSION
from .models import AtomRecord, CompatibilityReport
from .pdbio import infer_element, is_nucleic_base_name


def format_atom_name(name: str, element: str) -> str:
    n = name.strip()
    if len(n) >= 4:
        return n[:4]
    if len(element.strip()) <= 1 and n and not n[0].isdigit():
        return f"{n:>4}"
    return f"{n:<4}"


def format_pdb_atom(a: AtomRecord, out_serial: int) -> str:
    rec = a.record if a.record in {"ATOM", "HETATM"} else "ATOM"
    atom_field = format_atom_name(a.atom_name_out, a.element)
    res_field = f"{a.resname_out:>3}"[-3:]
    chain = (a.chain_out or "A")[:1]
    element = (a.element or infer_element(a.atom_name_out, resname=a.resname_out)).upper()[:2]
    segid = (a.segid or ("RNA" + chain if is_nucleic_base_name(a.resname_out) else chain))[:4]
    line = (
        f"{rec:<6}{out_serial:5d} {atom_field}{a.altloc[:1]:1s}"
        f"{res_field:>3s} {chain:1s}{a.resseq_out:4d}{a.icode_orig[:1]:1s}   "
        f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occ:6.2f}{a.bfac:6.2f}"
        f"      {segid:<4s}{element:>2s}{a.charge:>2s}"
    )
    if len(line) != 80:
        raise SystemExit(f"PDB field overflow for output atom {out_serial} ({a.resname_out}{a.resseq_out}/{a.atom_name_out}).")
    return line


def format_ter(serial: int, last_atom: AtomRecord) -> str:
    return f"TER   {serial:5d}      {last_atom.resname_out:>3s} {last_atom.chain_out[:1]:1s}{last_atom.resseq_out:4d}"


def write_pdb(atoms: List[AtomRecord], path: str) -> None:
    if not atoms:
        raise SystemExit("Cannot write an empty PDB template.")
    for atom in atoms:
        if len(atom.resname_out) > 3:
            raise SystemExit(f"Residue name {atom.resname_out!r} exceeds three PDB columns.")
        if len(atom.charge) > 2:
            raise SystemExit(f"Charge field {atom.charge!r} exceeds two PDB columns.")
    with open(path, "w", encoding="utf-8") as out:
        out.write("REMARK Prepared by pdb2plmd\n")
        out.write(f"REMARK pdb2plmd version {VERSION}\n")
        out.write("REMARK Atom order preserved from selected ATOM/HETATM input order\n")
        prev_atom = None
        for i, a in enumerate(atoms, start=1):
            if prev_atom is not None and a.chain_out != prev_atom.chain_out:
                out.write(format_ter(i, prev_atom) + "\n")
            out.write(format_pdb_atom(a, i) + "\n")
            prev_atom = a
        out.write(format_ter(len(atoms) + 1, atoms[-1]) + "\n")
        out.write("END\n")


def default_log_path(output_path: str) -> str:
    p = Path(output_path)
    return str(p.with_suffix(".log")) if p.suffix else str(Path(str(p) + ".log"))


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def atomic_write_pdb(atoms: List[AtomRecord], path: str) -> None:
    dest = Path(path)
    tmp = dest.with_name(f".{dest.name}.pdb2plmd.{os.getpid()}.tmp")
    try:
        write_pdb(atoms, str(tmp))
        os.replace(tmp, dest)
    finally:
        try:
            if tmp.exists():
                tmp.unlink()
        except OSError:
            pass


def _separator() -> str:
    return "=" * 67


def _write_notes(log, args: argparse.Namespace) -> None:
    log.write("Notes:\n")
    if args.serials is not None:
        log.write("- The -s range is applied to the PDB atom-serial field after MODEL, altLoc and solvent filtering.\n")
        log.write("- TER and other non-coordinate records are not selected; repeated selected atom serials are rejected as ambiguous.\n")
    else:
        log.write("- The -a range is applied after MODEL, altLoc and solvent filtering.\n")
        log.write("- It refers to ATOM/HETATM order, not PDB serial.\n")
    log.write("- Output atom order is identical to the selected input atom order.\n")
    log.write("- Output atom serials are renumbered sequentially; this does not change PLUMED atom order.\n")
    log.write("- Residues are renumbered sequentially within each output chain to avoid PLUMED residue-range gaps.\n")


def write_success_log(
    log_path: str,
    args: argparse.Namespace,
    all_atoms: List[AtomRecord],
    selected_indices: List[int],
    converted: List[AtomRecord],
    report: CompatibilityReport,
    log_lines: List[str],
) -> None:
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(f"pdb2plmd {VERSION} log\n")
        log.write(_separator() + "\n")
        log.write(f"Timestamp UTC: {utc_timestamp()}\n")
        log.write(f"Python: {sys.executable}\n")
        log.write(f"Command: {shlex.join(sys.argv)}\n")
        log.write(f"Input:  {args.input}\n")
        log.write(f"Output: {args.output}\n")
        if args.serials is not None:
            log.write(f"PDB serial range expression: {args.serials}\n")
            log.write("Selection mode: PDB atom serial (-s/--serials)\n")
        else:
            log.write(f"Atom-order range expression: {args.atoms if args.atoms is not None else 'all'}\n")
            log.write("Selection mode: ATOM/HETATM record order (-a/--atoms)\n")
        log.write(f"MODEL: {args.model if args.model is not None else 'first encountered'}\n")
        log.write(f"Input convention option: {args.input_convention}\n")
        if args.charmm is not None:
            log.write(f"Legacy format override: {'CHARMM forced' if args.charmm else 'CHARMM disabled'}\n")
        log.write(f"Split on gaps: {args.split_on_gaps}\n")
        log.write(f"AltLoc selection: {args.altloc}\n")
        log.write(f"Drop solvent: {args.drop_solvent}\n")
        log.write(f"Total ATOM/HETATM records after MODEL/altLoc/solvent filtering: {len(all_atoms)}\n")
        log.write(f"Selected atoms: {len(selected_indices)}\n")
        log.write(f"Output atoms: {len(converted)}\n")
        if selected_indices:
            log.write(f"Selected input atom-order range: {selected_indices[0]}..{selected_indices[-1]}\n")
        log.write(_separator() + "\n")
        log.write("Compatibility:\n")
        log.write(f"PDB_PARSEABLE {report.pdb_parseable}\n")
        log.write(f"MOLINFO_COMPATIBLE {report.molinfo_compatible}\n")
        log.write(f"ONEBEAD_MAPPABLE {report.onebead_mappable}\n")
        log.write(f"ONEBEAD_SAXS_READY {report.onebead_saxs_ready}\n")
        if report.readiness_issues:
            log.write("Readiness issues:\n")
            for item in report.readiness_issues:
                log.write(f"- {item}\n")
        if report.warnings:
            log.write("Warnings:\n")
            for item in report.warnings:
                log.write(f"- {item}\n")
        log.write(_separator() + "\n")
        log.write("Conversion details:\n")
        for line in log_lines:
            log.write(line + "\n")
        log.write(_separator() + "\n")
        _write_notes(log, args)


def write_error_log(
    log_path: str,
    args: argparse.Namespace,
    error_text: str,
    context_lines: Sequence[str],
    unexpected_traceback: str = "",
    output_existed_before: bool = False,
) -> str:
    requested = Path(log_path)
    candidates = [requested]
    fallback = Path.cwd() / requested.name
    if fallback != requested:
        candidates.append(fallback)
    last_error: Optional[Exception] = None
    for target in candidates:
        try:
            with target.open("w", encoding="utf-8") as log:
                log.write(f"pdb2plmd {VERSION} log\n")
                log.write(_separator() + "\n")
                log.write(f"Timestamp UTC: {utc_timestamp()}\n")
                log.write(f"Python: {sys.executable}\n")
                log.write(f"Command: {shlex.join(sys.argv)}\n")
                log.write(f"Input:  {args.input}\n")
                log.write(f"Output: {args.output}\n")
                log.write(_separator() + "\n")
                log.write("Conversion failed:\n")
                log.write((error_text or "Unknown conversion error").rstrip() + "\n")
                if output_existed_before:
                    log.write("Existing output file was left unchanged.\n")
                else:
                    log.write("Output PDB was not written.\n")
                if context_lines:
                    log.write(_separator() + "\n")
                    log.write("Context before failure:\n")
                    for line in context_lines:
                        log.write(line + "\n")
                if unexpected_traceback:
                    log.write(_separator() + "\n")
                    log.write("Unexpected exception traceback:\n")
                    log.write(unexpected_traceback.rstrip() + "\n")
            return str(target)
        except OSError as exc:
            last_error = exc
    raise OSError(f"could not write error log {requested}: {last_error}")
