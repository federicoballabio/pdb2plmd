from __future__ import annotations

CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

RNA_MAP = {
    "A": "A", "ADE": "A", "RA": "A", "RAD": "A", "R_A": "A",
    "C": "C", "CYT": "C", "RC": "C", "RCY": "C", "R_C": "C",
    "G": "G", "GUA": "G", "RG": "G", "RGU": "G", "R_G": "G",
    "U": "U", "URA": "U", "RU": "U", "URI": "U", "R_U": "U",
}
DNA_MAP = {
    "DA": "DA", "ADE_D": "DA", "DADE": "DA", "DAD": "DA",
    "DC": "DC", "CYT_D": "DC", "DCYT": "DC", "DCD": "DC",
    "DG": "DG", "GUA_D": "DG", "DGUA": "DG", "DGD": "DG",
    "DT": "DT", "THY": "DT", "DTH": "DT", "DTHY": "DT",
}
PROTEIN_NAMES = {
    "ALA","ARG","ASN","ASP","ASH","CYS","CYX","GLN","GLU","GLH","GLY",
    "HIS","HID","HIE","HIP","HSD","HSE","HSP","ILE","LEU","LYS","MET",
    "PHE","PRO","SER","THR","TRP","TYR","VAL"
}
ION_MAP = {
    "NA": "NA", "SOD": "NA", "NA+": "NA",
    "K": "K", "POT": "K", "K+": "K",
    "CL": "CL", "CLA": "CL", "CL-": "CL",
    "CA": "CAL", "CAL": "CAL", "CA2+": "CAL",
    "MG": "MG", "MG2+": "MG",
    "ZN": "ZN", "ZN2": "ZN", "ZN2+": "ZN",
    "FE2": "FE2", "FE3": "FE3",
    "MN": "MN", "MN2+": "MN",
}
ION_ELEMENT = {
    "NA": "NA", "K": "K", "CL": "CL", "CAL": "CA", "MG": "MG",
    "ZN": "ZN", "FE2": "FE", "FE3": "FE", "MN": "MN",
}

UNSUPPORTED_ION_NAMES = {
    "LI": "lithium", "LIT": "lithium",
    "CS": "cesium", "CES": "cesium",
    "RB": "rubidium", "RUB": "rubidium",
    "BA": "barium", "BAR": "barium",
    "CD": "cadmium", "CD2": "cadmium",
    "CU": "copper", "CU1": "copper(I)", "CU2+": "copper(II)",
}

PENTOSE_ATOMS = {
    "O5'","C5'","O4'","C4'","O3'","C3'","O2'","C2'","C1'",
    "H5'","H5''","H4'","H3'","H2'","H2''","H2'2","H1'",
    "HO5'","HO3'","HO2'","H5'1","H5'2","HO'2","H2'1","H5T","H3T",
}
BASE_ATOMS = {
    "N1","N2","N3","N4","N6","N7","N9","C2","C4","C5","C6","C7","C8",
    "O2","O4","O6","H1","H2","H3","H5","H6","H8","H21","H22","H41","H42",
    "H61","H62","H71","H72","H73",
}
PHOSPHATE_ATOMS = {"P","OP1","OP2","OP3","O1P","O2P","O3P","HP","HOP3"}
KNOWN_ONEBEAD_NUC_ATOMS = PENTOSE_ATOMS | BASE_ATOMS | PHOSPHATE_ATOMS

_GLC_FAMILY = {"GLC", "BGC", "AGLC", "BGLC"}
_GAL_FAMILY = {"GAL", "GLA", "AGAL", "BGAL"}
_MAN_A      = {"MAN", "AMAN"}
_MAN_B      = {"BMA", "BMAN"}
_FUC_FAMILY = {"FUC", "FUL", "FCA", "FCB", "AFUC", "BFUC"}
_NEU_FAMILY = {"SIA", "SLB", "ANE5AC", "BNE5AC", "ANE5", "BNE5"}
_NAG_NAMES  = {"NAG", "NDG", "AGLCNA", "BGLCNA"}
_NGA_NAMES  = {"NGA", "A2G", "AGALNA", "BGALNA"}

GLYCAN_UNSUPPORTED = {
    "RAM", "RM4", "XXR", "XYL", "XYS", "XYP", "LXZ", "GCU", "BDP", "GCV", "IDS", "IDR", "SGN", "SUS",
    "UAP", "NGC", "NGE", "RIB", "ARA", "ARB", "AHR", "GLP", "PA1", "GCS",
}
GLYCOSYLATION_SITE_MAP = {"NLN": "ASN", "OLS": "SER", "OLT": "THR"}
UNSUPPORTED_SITE_NAMES = {"OLP"}
TERMINAL_CAP_NAMES = {"ROH", "OME", "TBT"}
MAX_REPORTED_ISSUE_GROUPS = 20

GLYCAN_ATOM_MAP_HEX = {
    "N": "N2", "C": "C7", "O": "O7", "CT": "C8",
    "HN": "HN2", "HT1": "H81", "HT2": "H82", "HT3": "H83",
}
GLYCAN_ATOM_MAP_SIA = {
    "N": "N5", "C": "C10", "O": "O10", "CT": "C11",
    "O11": "O1A", "O12": "O1B",
    "HN": "HN5", "HT1": "H111", "HT2": "H112", "HT3": "H113",
}

_GLYCAN_REQUIRED_HEAVY = {
    "GLC": {"C1","C2","C3","C4","C5","C6","O2","O3","O4","O5","O6"},
    "GAL": {"C1","C2","C3","C4","C5","C6","O2","O3","O4","O5","O6"},
    "MAN": {"C1","C2","C3","C4","C5","C6","O2","O3","O4","O5","O6"},
    "BMA": {"C1","C2","C3","C4","C5","C6","O2","O3","O4","O5","O6"},
    "FUC": {"C1","C2","C3","C4","C5","C6","O2","O3","O4","O5"},
    "NAG": {"C1","C2","C3","C4","C5","C6","C7","C8","N2","O3","O4","O5","O6","O7"},
    "NGA": {"C1","C2","C3","C4","C5","C6","C7","C8","N2","O3","O4","O5","O6","O7"},
    "SIA": {"C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","N5","O1A","O1B","O4","O6","O7","O8","O9","O10"},
}
_GLYCAN_ALLOWED_HEAVY = {
    bead: required | ({"O1"} if bead != "SIA" else {"O2"})
    for bead, required in _GLYCAN_REQUIRED_HEAVY.items()
}

SOLVENT_AND_ADDITIVES = {
    "HOH", "WAT", "TIP3", "SOL", "DOD",
    "EDO", "GOL", "PEG", "PG4", "SO4", "PO4", "CIT", "EPE", "TLA", "MES", "TRS",
    "ACT", "DMS", "IMD", "FMT", "NO3", "MPD", "BME", "EOH",
}
