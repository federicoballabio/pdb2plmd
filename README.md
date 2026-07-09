# pdb2plmd

Prepare a PDB file for PLUMED-SAS-ONEBEAD calculations.

The main goal is to preserve the atom order from a PDB extracted from a
simulation or from a selected atom group. This is required because PLUMED maps
the atoms in the ATOMS list to the atoms in the TEMPLATE PDB by order.

## Main features

- Preserves the selected ATOM/HETATM order.
- Renumbers output atom serials from 1.
- Assigns chain IDs if missing.
- Can infer chain IDs from CHARMM-GUI segment IDs such as RNAA, RNAB, PROA, PROB.
- Inserts TER records between inferred chains or breaks.
- Renumbers residues sequentially within each output chain.
- Converts common RNA residue names to SAXS.cpp/AMBER-style names.
- Converts common DNA residue names to SAXS.cpp/AMBER-style names.
- Adds RNA/DNA terminal suffixes when detectable.
- Normalizes old nucleic-acid atom names using * to apostrophe notation.
- Optionally writes a verbose conversion log.

## Requirements

Python 3.8 or newer is recommended.

No external Python packages are required.

## Usage

Basic usage:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb
```

Select a range of atoms by input ATOM/HETATM order:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-1062
```

Select multiple ranges:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-100,150,200-250
```

Write a verbose log using the default log name:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-1062 -g
```

Write a verbose log using an explicit log name:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-1062 -g conversion.log
```

## Options

```text
-i, --input   Input PDB file.
-o, --output  Output PDB file for PLUMED SAXS.cpp TEMPLATE.
-a, --atoms   Atom range to keep. Default: all.
-g, --log     Write a verbose log. Optional file name.
```

## Important note about atom order

The `-a` option uses the 1-based order of ATOM/HETATM records in the input PDB.
It does not use the PDB atom serial number.

For example, `-a 1-1062` means: keep the first 1062 ATOM/HETATM records found in
the input file.

This is intentional. For PLUMED SAXS.cpp ONEBEAD, the output TEMPLATE PDB must
have the same atom order as the atoms selected by PLUMED in the ATOMS field.

## Recommended workflow

1. Extract the SAXS atom group from the simulation using the same order used by
   the PLUMED ATOMS selection.
2. Run `pdb2plmd.py` on that extracted PDB.
3. Use the generated PDB as the PLUMED SAXS.cpp TEMPLATE file.
4. Use the same atom count and order in the PLUMED ATOMS field.

## Nucleic-acid naming

The script converts common RNA names:

```text
ADE -> A
CYT -> C
GUA -> G
URA -> U
RA  -> A
RC  -> C
RG  -> G
RU  -> U
```

The script converts common DNA names:

```text
DADE -> DA
DCYT -> DC
DGUA -> DG
THY  -> DT
DA   -> DA
DC   -> DC
DG   -> DG
DT   -> DT
```

Terminal suffixes are added when detectable:

```text
C5, U5, A5, G5   5-prime hydroxyl terminal RNA residue
C3, U3, A3, G3   3-prime hydroxyl terminal RNA residue
CT, UT, AT, GT   5-prime phosphorylated terminal RNA residue
DC5, DG5, DA5, DT5   5-prime hydroxyl terminal DNA residue
DC3, DG3, DA3, DT3   3-prime hydroxyl terminal DNA residue
DCT, DGT, DAT, DTT   5-prime phosphorylated terminal DNA residue
```

## Output

The output PDB contains:

- REMARK lines documenting that the file was prepared by `pdb2plmd.py`.
- ATOM/HETATM records in the same selected order as the input.
- Sequential output atom serials.
- Sequential residue numbers within each output chain.
- TER records between chains and at the end.
- END at the end of the file.

## Validation checks

After conversion, check the atom count:

```bash
grep -E "^(ATOM|HETATM)" template_saxs.pdb | wc -l
```

Check the first and last atoms:

```bash
grep -E "^(ATOM|HETATM)" template_saxs.pdb | head
grep -E "^(ATOM|HETATM)" template_saxs.pdb | tail
```

If a log was generated, inspect warnings:

```bash
grep WARNING template_saxs.pdb.log
```

## Example

```bash
python3 pdb2plmd.py \
  -i template_AA.pdb \
  -o template_AA_saxs.pdb \
  -a 1-1062 \
  -g
```

Use in PLUMED:

```plumed
MOLINFO STRUCTURE=template_AA_saxs.pdb

saxsdata: SAXS ...
  ATOMS=1-1062
  ONEBEAD
  TEMPLATE=template_AA_saxs.pdb
... SAXS
```
