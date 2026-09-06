from __future__ import annotations
from dataclasses import dataclass, field
from typing import List

@dataclass
class AtomRecord:
    record: str
    input_atom_index: int
    input_serial: int
    atom_name: str
    altloc: str
    resname_orig: str
    chain_orig: str
    resseq_orig: int
    icode_orig: str
    x: float
    y: float
    z: float
    occ: float
    bfac: float
    segid: str
    element: str
    charge: str
    line_number: int = 0
    ter_before: bool = False
    chain_out: str = ""
    resseq_out: int = 0
    resname_out: str = ""
    atom_name_out: str = ""



    source_convention: str = ""
    chemical_identity: str = ""
    chemical_state: str = ""
    glycan_anomer: str = ""
    glycan_config: str = ""
    glycan_ring: str = ""
    glycan_linkage: str = ""
    glycan_ccd: str = ""
    mapping_provenance: str = ""

@dataclass(frozen=True)
class CompatibilityIssue:
    key: str
    message: str
    hint: str = ""

@dataclass(frozen=True)
class PdbFormatPreflight:
    sampled_atoms: int
    chain_records: int
    segid_records: int
    four_char_residue_records: int
    style: str

@dataclass
class CompatibilityReport:
    requested_convention: str
    resolved_convention: str
    convention_evidence: List[str] = field(default_factory=list)
    readiness_issues: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    has_nucleic: bool = False

    @property
    def pdb_parseable(self) -> str:
        return "PASS"

    @property
    def molinfo_compatible(self) -> str:
        return "PASS"

    @property
    def onebead_mappable(self) -> str:
        return "PASS"

    @property
    def onebead_saxs_ready(self) -> str:
        if self.readiness_issues:
            return "FAIL"
        if self.warnings:
            return "PASS_WITH_WARNING"
        return "PASS"
