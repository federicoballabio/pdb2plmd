# pdb2plmd

Prepare and validate an ordered PDB template for PLUMED `SAXS ONEBEAD` calculations.

## Main features

- Preserves selected `ATOM`/`HETATM` order and coordinates.
- Renumbers output atom serials from 1 and residues sequentially within each chain.
- Preserves PDB chain IDs and can recover CHARMM-GUI chains from SEGID fields.
- Resolves alternate locations at residue level.
- Supports proteins, RNA, DNA, the validated glycan classes and the ion set implemented by the current `SAXS.cpp`.
- Converts supported CHARMM, AMBER, GLYCAM and force-field-qualified GROMACS nomenclatures.
- Separates successful mapping from quantitative ONEBEAD SAXS readiness.
- Accepts heavy-atom-complete structures without explicit hydrogens with a warning.
- Flags missing required heavy atoms as not SAXS-ready and rejects unsupported chemical states and incomplete multi-atom selections.
- Writes concise diagnostics and an error log when conversion fails.

## Requirements

Python 3.8 or newer. No external Python packages are required.

Keep `pdb2plmd.py` and the `src` folder together.

## Usage

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb
```

For a source convention that should not be guessed automatically, specify it explicitly:

```bash
python3 pdb2plmd.py -i input.pdb -o template_saxs.pdb --input-convention amber -v
```

Check the exact converter version with:

```bash
python3 pdb2plmd.py --version
```

Current version:

```text
pdb2plmd v2609061930
```

## Options

```text
-i, --input          Input PDB file.
-o, --output         Output PDB file.
-a, --atoms          Select by 1-based ATOM/HETATM record order.
-s, --serials        Select by the PDB atom-serial field.
--model              MODEL serial to keep. Default: first MODEL.
--input-convention   Source nomenclature convention.
--split-on-gaps      Split chains at residue-number resets or gaps.
-altloc, --altloc    Alternate-location handling. Default: auto.
--drop-solvent       Remove water and common crystallisation additives.
-v, --verbose        Print detailed terminal status and write a log on success.
--version            Print the exact converter version.
```

Accepted source-convention values are:

```text
auto
generic
charmm
amber
gromacs-amber
gromacs-charmm
gromacs-oplsaa
gromacs-gromos
glycam
```

`auto` uses conservative PDB signatures and resolves a source convention only when the evidence is unambiguous. It can recognize CHARMM/CHARMM-GUI layouts, distinctive AMBER residue forms, strong or multiple GLYCAM residue/atom signatures, and selected force-field-qualified GROMACS signatures. A single ambiguous residue-name match is not sufficient. If GROMACS-style naming is detected but the force-field family cannot be distinguished safely, conversion stops and asks for an explicit `--input-convention`. GROMOS united-atom input without a distinctive source signature should also be selected explicitly.

## Atom order and selection

`-a` selects by 1-based `ATOM`/`HETATM` record order after MODEL, alternate-location and optional solvent handling. `TER` and other non-coordinate records are not counted.

`-s` selects by the PDB atom-serial field. Use `-a` for files with repeated or wrapped atom serials. `-a` and `-s` are mutually exclusive.

The converter does not sort atoms. A selection that cuts through a residue or multi-atom moiety is rejected.

## Alternate locations

By default, `--altloc auto` selects one coherent nonblank alternate-location label per residue using completeness first and occupancy second. Blank-altLoc atoms are retained as common atoms. The output altLoc field is blank.

A label can be forced globally with `--altloc A`, `--altloc B` or another label. If the requested label is unavailable for a residue with alternate coordinates, conversion stops.

## Compatibility levels

Verbose output reports:

```text
PDB_PARSEABLE
MOLINFO_COMPATIBLE
ONEBEAD_MAPPABLE
ONEBEAD_SAXS_READY
```

A structure can be mappable but not suitable for quantitative ONEBEAD SAXS. In particular, recognized GROMOS united-atom input is deliberately reported as not SAXS-ready because its implicit-hydrogen representation has not been quantitatively validated for the current ONEBEAD model.

A heavy-atom-complete structure without explicit hydrogens can be reported as `PASS_WITH_WARNING`. Missing required heavy atoms give `ONEBEAD_SAXS_READY=FAIL`.

## Proteins

The 20 standard amino acids are supported. Histidine tautomer/protonation names used by the supported source conventions are retained or translated according to the current ONEBEAD state definitions. `HIP`/`HSP` use the protonated-histidine bead and `CYX` uses the disulfide-cysteine bead.

The accepted SAXS approximations are:

```text
GLH -> GLU
ASH -> ASP
NLN -> ASN
OLS -> SER
OLT -> THR
```

Unsupported protonation or modification states are not replaced by a chemically similar standard residue.

## Nucleic acids

RNA and DNA are converted to the strict residue and atom vocabulary used by the current ONEBEAD implementation. Supported source translations include CHARMM and AMBER/GROMACS-AMBER aliases, terminal residue forms and legacy atom names.

Ambiguous older `ADE`, `CYT` and `GUA` names are resolved as RNA or DNA from the presence of `O2'`. `THY` is accepted only for deoxy thymine and `URA` only for ribose uracil.

AMBER monomer building-block names without a dedicated current ONEBEAD state remain unsupported.

## Glycans

The current ONEBEAD glycan classes are:

```text
FUC
MAN
BMA
GAL
GLC
NAG
NGA
SIA
```

Supported CHARMM and common PDB aliases are converted to these targets. Ambiguous glucose- or galactose-derived names are promoted to `NAG` or `NGA` only when the required N-acetyl chemistry is present.

`glycam` mode uses the bundled GLYCAM code table. GLYCAM residue-name case is preserved because it carries stereochemical information. Recognized sugars outside the eight validated ONEBEAD classes are rejected rather than substituted.

## Ions

The converter supports the ion aliases implemented by the current `SAXS.cpp` ion resolver. A monoatomic ion residue must contain exactly one coordinate atom. Recognized ions without a current ONEBEAD parameter are reported explicitly and are not substituted.

## Solvent

`--drop-solvent` removes the implemented water and common crystallisation-additive residue names before compatibility checks. Without this option, unsupported solvent or additive residues stop conversion.

## Output and logs

Without `-v`, a successful conversion writes only the output PDB and reports:

```text
pdb2plmd v2609061930: PDB converted.
```

With `-v`, the same conversion also writes `<output-stem>.log` and reports the output path, atom count, `ONEBEAD_SAXS_READY` status and log path on screen.

Without `-v`, a failed conversion reports:

```text
pdb2plmd v2609061930: PDB conversion failed.
```

An error log is written on failure in both modes. With `-v`, the terminal also shows the first error and the error-log path. No converted PDB is written on failure. Warnings remain visible on standard error when relevant.

The output PDB records the exact `pdb2plmd` version used to generate it.

## Examples

`examples/run_examples.py` uses 10 real structures from the Protein Data Bank. The exact PDB files are bundled in `examples/fixtures`. Successful reference outputs are stored in `examples/expected`.

The current set is:

```text
1CRN  protein
1UBQ  protein with crystallographic water
1CLL  calmodulin with calcium and ethanol
3Q2W  cadherin with glycans and calcium
1BNA  DNA duplex
6M0J  glycoprotein complex with NAG, zinc and chloride
6VSB  SARS-CoV-2 spike with NAG
1EHZ  tRNA with modified nucleotides
4HHB  hemoglobin with heme
4COF  GABAA receptor with benzamidine
```

Some examples are expected to convert, including structures that remain `ONEBEAD_SAXS_READY=FAIL` because deposited coordinates are incomplete. Other examples are expected to fail conversion because they contain chemistry outside the current ONEBEAD model.

Run all examples from the package directory with:

```bash
python3 examples/run_examples.py
```

Generated PDB and log files are written to `examples/output`.