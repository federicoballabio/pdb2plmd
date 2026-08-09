# pdb2plmd

Prepare an ordered PDB template for PLUMED calculations.

PLUMED maps the atoms in `ATOMS` to the atoms in `TEMPLATE` by order. The
converter therefore preserves the selected `ATOM`/`HETATM` order while
normalizing residue, atom and chain names required by ONEBEAD.

## Main features

- Preserves the selected `ATOM`/`HETATM` order.
- Renumbers output atom serials from 1.
- Assigns chain IDs if missing and can infer them from CHARMM-GUI segment IDs.
- Detects CHARMM four-character residue names and reads them from columns 18-21.
- Inserts `TER` records and renumbers residues sequentially within each chain.
- Converts supported glycosylation-site, glycan, ion, RNA and DNA names.
- Adds RNA/DNA terminal suffixes when detectable.
- Resolves alternate locations at residue level. Default `--altloc auto` selects one coherent conformer by completeness, then occupancy; `--altloc A/B/...` can force one label globally.
- Reports incompatible residues or moieties with concise, actionable errors.
- Rejects `-a`/`-s` selections that keep only part of a residue or multi-atom moiety, and warns when the template appears not to be all-atom.
- Embeds the converter version in the PDB and exposes it with `--version`.
- Uses a two-mode file interface: standard success writes only the PDB; `-v` success writes PDB + log; any conversion error writes an error log.

## Requirements

Python 3.8 or newer. No external Python packages are required.

## Usage

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb
```

Select atoms by filtered `ATOM`/`HETATM` order:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-1062
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -a 1-100,150,200-250
```

Select by the PDB atom-serial field instead:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -s 1-1069
```

Verbose run:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb -v
```

Check the exact converter version:

```bash
python3 pdb2plmd.py --version
# pdb2plmd v26080804
```

For `template_saxs.pdb`, the automatic log name is `template_saxs.log`.

## Options

```text
-i, --input          Input PDB file.
-o, --output         Output PDB file.
-a, --atoms          Select by 1-based ATOM/HETATM record order. Default: all.
-s, --serials        Select by the PDB atom-serial field. Mutually exclusive with -a.
--model              MODEL serial to keep. Default: 1.
--split-on-gaps      Split chains at residue-number resets or gaps.
-altloc, --altloc    Select alternate locations. Default: auto.
--drop-solvent       Remove water and common crystallisation additives, including EOH ethanol.
-v, --verbose        On success write a verbose log in addition to the PDB.
--version             Print the exact pdb2plmd version and exit.
```

## Atom order and selection

`-a` selects by 1-based `ATOM`/`HETATM` record order after MODEL, alternate-location resolution and optional solvent filtering. TER and other non-coordinate records are not counted.

`-s` selects by the PDB atom-serial field. This is convenient when the desired boundary is known from the PDB itself; TER serials are naturally skipped. If selected coordinate records contain repeated or wrapped PDB serials, `-s` stops with an error and `-a` should be used instead. `-a` and `-s` are mutually exclusive.

The script does not sort atoms. The output `TEMPLATE` keeps the selected input atom order. A selection that cuts through a residue or multi-atom moiety is rejected before conversion.

## Alternate locations

By default, `--altloc auto` resolves alternate coordinates once per residue. Blank-altLoc atoms are retained as common atoms. For each residue with alternate coordinates, the converter chooses one nonblank label using this deterministic order:

1. most complete residue;
2. highest summed occupancy when completeness is equal;
3. label `A` on an exact tie;
4. lexical label order as the final tie-break.

The converter never mixes atoms from different nonblank altLoc labels in one output residue. The output altLoc column is blank because the converted PDB represents one resolved structure.

A label can be forced globally, for example:

```bash
python3 pdb2plmd.py -i input.pdb -o template.pdb --altloc B
```

`-altloc B` is accepted as an alias. If a residue has alternate coordinates but the requested label is unavailable, conversion stops instead of silently falling back to another conformer.

## Compatibility diagnostics

Known incompatible residues and moieties stop the conversion before the output
PDB is written. Fatal checks include:

- unsupported or unknown residue names;
- monosaccharides without ONEBEAD parameters;
- GLYCAM terminal caps;
- `OLP` hydroxyproline;
- supported glycans carrying unexpected elements, indicating an unsupported
  substituted moiety;
- ion residues containing more than one atom;
- unsupported nucleic-acid atom names;
- duplicate atom names after normalization;
- PDB field overflows.

Repeated failures are grouped by issue and residue name. At most 20 issue groups
are printed, followed by relevant hints. An RCSB structure containing many
waters therefore produces one grouped water error rather than hundreds of
nearly identical lines.

Warnings do not stop conversion. They are printed to standard error and copied
to the verbose log when `-v` is used. The converter checks the residue/moiety rules implemented
here; the installed `SAXS.cpp` remains the final authority for exact ONEBEAD and
LCPO compatibility.

## Solvent and crystallisation additives

`--drop-solvent` removes the implemented water/additive residue names before compatibility checks. `EOH` is included as ethanol. `EOH` is not given a ONEBEAD mapping: without `--drop-solvent`, it is reported as an unsupported solvent/additive and conversion stops.

## Glycans

Glycans are kept as one residue per monosaccharide. Common PDB and force-field
names are converted as follows:

```text
GLC, BGC, AGLC, BGLC        -> GLC or NAG when an amide nitrogen is present
GAL, GLA, AGAL, BGAL        -> GAL or NGA when an amide nitrogen is present
MAN, AMAN                   -> MAN
BMA, BMAN                   -> BMA
FUC, FUL, FCA, FCB          -> FUC
AFUC, BFUC                  -> FUC
NAG, NDG, AGLCNA, BGLCNA   -> NAG
NGA, A2G, AGALNA, BGALNA   -> NGA
SIA, SLB, ANE5AC, BNE5AC   -> SIA
```

CHARMM writes six-character names that the PDB format truncates to four, so
`BGLCNA` arrives as `BGLC` and `AGLCNA` as `AGLC`. These are indistinguishable
from glucose by name alone and are resolved by the amide nitrogen. The same
truncation gives `ANE5` and `BNE5` for Neu5Ac; both are accepted.

Force-field atom names in N-acetyl groups and sialic acid are converted to the
corresponding PDB chemical-component names. A free reducing-end `O1` is accepted.

Monosaccharides without ONEBEAD parameters stop the conversion. These include
rhamnose, xylose, uronic acids, sulfated sugars, Neu5Gc, pentoses and free
glucosamine. No chemically similar sugar is substituted automatically.

GLYCAM terminal capping groups (`ROH`, `OME`, `TBT`) are not monosaccharides and
have no bead. They are rejected and must be removed before conversion.

## Junctions

A glycosylation-site residue is renamed to its parent amino acid:

```text
NLN -> ASN
OLS -> SER
OLT -> THR
```

This is correct for the implemented ONEBEAD model. All heavy atoms of `NLN`,
`OLS` and `OLT` are present in `ASN`, `SER` and `THR`, so the LCPO lookup remains
valid after renaming.

The glycan bead needs no extra correction. Its parameters were derived from
glycans simulated while bonded to a `GLY-ASN-GLY` or `ALA-THR-ALA` tripeptide,
which was removed before the form factor was computed. The reducing
monosaccharide therefore carries no `O1` in the parameter set, matching its
state in a glycoprotein.

On the protein side, one hydrogen attached to `ND2` in ASN, `OG` in SER or `OG1`
in THR is absent at the junction.

`OLP` is rejected. Hydroxyproline contains `OD1`, which PRO does not have, so
renaming it would leave an atom without an LCPO entry.

## Ions

Monoatomic ions are converted to the residue names expected by `SAXS.cpp`:

```text
NA, SOD, NA+     -> NA
K, POT, K+       -> K
CL, CLA, CL-     -> CL
CA, CAL, CA2+    -> CAL
MG, MG2+         -> MG
ZN, ZN2, ZN2+    -> ZN
FE2              -> FE2
FE3              -> FE3
MN, MN2+         -> MN
```

Calcium is written as `CAL` rather than `CA` to avoid confusion with the protein
alpha-carbon atom name. The element field is normalized to the corresponding
chemical element. An ion residue must contain exactly one atom.

`SAXS.cpp` excludes these ions from the SASA calculation.

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

For RNA residues, the legacy CHARMM pair `H2''` + `H2'` is converted to the canonical pair `H2'` + `HO2'` only when the residue contains both legacy names and does not already contain `HO2'`. Already-canonical `H2'` + `HO2'` residues are left unchanged. This makes re-conversion idempotent. The rule is not applied to DNA residues.

## Example

```bash
python3 pdb2plmd.py \
  -i template_AA.pdb \
  -o template_AA_saxs.pdb \
  -a 1-1062 \
  -v
```

```plumed
MOLINFO STRUCTURE=template_AA_saxs.pdb

SAXS ...
  LABEL=saxsdata
  ATOMS=1-1062
  ONEBEAD
  TEMPLATE=template_AA_saxs.pdb
... SAXS
```


## Automatic PDB format handling

Input layout is detected automatically. Users do not need to specify CHARMM or CHARMM-GUI mode. The converter performs a lightweight preflight, reads standard three-character and CHARMM-style four-character residue names automatically, and uses the PDB chain ID when present or a SEGID-derived fallback when needed. Detection details are written to the verbose log. Legacy `--charmm` / `--no-charmm` overrides are retained only for backward compatibility and are hidden from normal help.
