# pdb2plmd

Prepare a PDB file for PLUMED.

The main goal is to preserve the atom order from a PDB extracted from a
simulation or from a selected atom group. PLUMED maps the atoms in `ATOMS` to
the atoms in the `TEMPLATE` PDB by order.

## Main features

- Preserves the selected ATOM/HETATM order.
- Renumbers output atom serials from 1.
- Assigns chain IDs if missing.
- Can infer chain IDs from CHARMM-GUI segment IDs.
- Inserts TER records between chains.
- Renumbers residues sequentially within each output chain.
- Converts supported protein glycosylation-site aliases.
- Converts supported glycan residue and atom names.
- Converts common RNA and DNA residue names.
- Adds RNA/DNA terminal suffixes when detectable.
- Selects one alternate location per atom.
- Rejects unsupported residues and duplicate output atom names.
- Optionally writes a conversion log.

## Requirements

Python 3.8 or newer is recommended.

No external Python packages are required.

## Usage

Basic usage:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb
```

Select atoms by ATOM/HETATM order:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-1062
```

Select multiple ranges:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-100,150,200-250
```

Write a log:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -g
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -g conversion.log
```

## Options

```text
-i, --input          Input PDB file.
-o, --output         Output PDB file.
-a, --atoms          Atom range to keep. Default: all.
--model              MODEL serial to keep. Default: 1.
--charmm             Force CHARMM four-character residue-name parsing.
--no-charmm          Disable automatic CHARMM parsing.
--split-on-gaps      Split chains at residue-number resets or gaps.
--drop-solvent       Remove water and common crystallisation additives.
-g, --log            Write a log. Optional file name.
```

## Atom order

The `-a` option uses the 1-based ATOM/HETATM order after MODEL, alternate
location and optional solvent filtering. It does not use the PDB atom serial
number.

The script does not sort atoms. The output `TEMPLATE` must contain the same
atoms, in the same order, as the PLUMED `ATOMS` selection.

## Recommended workflow

1. Extract the SAXS atom group from the simulation in the order used by PLUMED.
2. Run `pdb2plmd.py` on that PDB.
3. Check the output atom count and conversion log.
4. Use the generated PDB as the SAXS.cpp `TEMPLATE`.

## Glycans

Glycans are kept as one residue per monosaccharide. Common PDB and force-field
names are converted as follows:

```text
GLC, BGC, AGLC, BGLC       -> GLC or NAG when an amide nitrogen is present
GAL, GLA, AGAL, BGAL       -> GAL or NGA when an amide nitrogen is present
MAN, AMAN                  -> MAN
BMA, BMAN                  -> BMA
FUC, FUL, FCA, FCB         -> FUC
AFUC, BFUC                 -> FUC
NAG, NDG, AGLCNA, BGLCNA  -> NAG
NGA, A2G, AGALNA, BGALNA  -> NGA
SIA, SLB, ANE5AC, BNE5AC  -> SIA
```

The CHARMM-truncated names `ANE5` and `BNE5` are also accepted for Neu5Ac.
Force-field atom names in N-acetyl groups and sialic acid are converted to the
corresponding PDB chemical-component names.

Supported glycosylation-site aliases are:

```text
NLN -> ASN
OLS -> SER
OLT -> THR
```

Monosaccharides without ONEBEAD parameters stop the conversion. These include
rhamnose, xylose, uronic acids, sulfated sugars, Neu5Gc, pentoses and free
glucosamine. `OLP` is also rejected because hydroxyproline has no dedicated
ONEBEAD residue mapping.

## Nucleic-acid residue naming

Common RNA names are converted to `A`, `C`, `G` and `U`. Common DNA names are
converted to `DA`, `DC`, `DG` and `DT`.

Terminal suffixes are added when detectable:

```text
A5, C5, G5, U5       5-prime hydroxyl terminal RNA residue
A3, C3, G3, U3       3-prime hydroxyl terminal RNA residue
AT, CT, GT, UT       5-prime phosphorylated terminal RNA residue
DA5, DC5, DG5, DT5   5-prime hydroxyl terminal DNA residue
DA3, DC3, DG3, DT3   3-prime hydroxyl terminal DNA residue
DAT, DCT, DGT, DTT   5-prime phosphorylated terminal DNA residue
```

## CHARMM RNA 2-prime hydrogen naming

For RNA residues, the script converts:

```text
H2'' -> H2'
H2'  -> HO2'
```

This conversion is not applied to DNA residues.

## Output

The output PDB contains:

- `REMARK Prepared by pdb2plmd`;
- ATOM/HETATM records in the selected input order;
- sequential atom and residue numbers;
- TER records between chains and at the end;
- a final END record.

## Validation

Check the atom count:

```bash
grep -E "^(ATOM|HETATM)" template_saxs.pdb | wc -l
```

Check the first and last atoms:

```bash
grep -E "^(ATOM|HETATM)" template_saxs.pdb | head
grep -E "^(ATOM|HETATM)" template_saxs.pdb | tail
```

Run the tests:

```bash
python3 -m unittest discover -s tests -v
```

## Example

```bash
python3 pdb2plmd.py \
  -i template_AA.pdb \
  -o template_AA_saxs.pdb \
  -a 1-1062 \
  -g
```

```plumed
MOLINFO STRUCTURE=template_AA_saxs.pdb

saxsdata: SAXS ...
  ATOMS=1-1062
  ONEBEAD
  TEMPLATE=template_AA_saxs.pdb
... SAXS
```
