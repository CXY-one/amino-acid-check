[简体中文](README.md)

# Amino Acid Check

Amino Acid Check helps you inspect side chains before modeling an enzyme. It extracts local PDB structures for His, Asp, Glu, Asn, and Gln with their neighboring residues. Use the offline 3D viewer to compare candidates, inspect the evidence, measure distances, and record your review.

## Download and run

Download the 64-bit Windows ZIP from [GitHub Releases](https://github.com/CXY-one/amino-acid-check/releases/latest). Extract the entire archive on Windows 10 or 11 and run `Amino Acid Check.exe`. Keep `_internal` next to the executable. The ZIP includes the runtime. An existing PyMOL installation can optionally open the current candidate as sticks.

The interface includes Simplified Chinese and English. Go to **Settings → Interface Language**, select English, and click **Apply**. Use **Open Language Folder** to view or add JSON translations.

## Inspect a structure

1. Drop a `.pdb` or `.ent` file into the window, choose an output location, and run the analysis. Settings control the residue types, target chains, and neighborhood distances.
2. Filter the results by site, residue, status, decision, or pKa. Compare `HID`, `HIE`, and `HIP` for His, protonation states for Asp/Glu, and amide orientations for Asn/Gln. Hover over a cell to read its full text.
3. Use the 3D viewer to inspect the local environment. Drag with the left button to rotate, drag with the middle button to pan, and scroll to zoom. Turn on distance measurement and select two atoms.
4. Confirm a candidate to save your decision and advance to the next site. Export an Excel workbook with columns for site, residue, status, and decision, on a separate sheet for each chain. You can reopen a result folder to resume review.

Automatic assignments use local structure rules and are intended to support model preparation. The pKa column contains optional PROPKA estimates for His, Asp, and Glu. A dash indicates a missing assignment or estimate.

## Command line

The Windows ZIP includes `hischeck-cli.exe`. This example screens His in chain A and writes all three candidates and preview images where their geometry permits. Use a new output directory.

```powershell
.\hischeck-cli.exe .\input.pdb --chain A --residue HIS --all-states --render -o .\output
```

Add `--propka` for pKa estimates or repeat `--residue` to include other residue types. `summary.csv` lists the assignments. Each site directory has an `evidence.json` file and gains a `review.json` file when you confirm a selection.

See the [developer guide](docs/DEVELOPMENT.md) for setup, testing, translations, and builds. The original source is under the [MIT License](LICENSE). See the [third-party notices](packaging/windows/THIRD-PARTY.txt) for bundled components.
