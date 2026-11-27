[简体中文](README.md)

# Amino Acid Check

Amino Acid Check helps you inspect side chains before modeling an enzyme. It extracts local PDB structures for His, Asp, Glu, Asn, and Gln with their neighboring residues. Use the offline 3D viewer to compare candidates, inspect the evidence, measure distances, and record your review.

## Download and run

Download the 64-bit Windows ZIP from [GitHub Releases](https://github.com/CXY-one/amino-acid-check/releases/latest). Extract the entire archive on Windows 10 or 11 and run `Amino Acid Check.exe`. Keep `_internal` next to the executable. The ZIP includes the runtime. An existing PyMOL installation can optionally open the current candidate as sticks.

The interface includes Simplified Chinese and English. Go to **Settings → Interface Language**, select English, and click **Apply**. Use **Open Language Folder** to view or add JSON translations.

The app checks GitHub for the latest stable release in the background at startup. You can also click **Check for Updates** in Settings. When an update is available, **Download Update** opens its release page. Download and extract the new ZIP to start using it.

## Inspect a structure

1. Drop a `.pdb` or `.ent` file into the window, choose an output location, and run the analysis. Settings control the residue types, target chains, and neighborhood distances.
2. The progress bar counts completed sites and saved images in separate stages. It also shows when PROPKA is estimating pKa.
3. Filter the results by site, residue, status, decision, or pKa. Compare `HID`, `HIE`, and `HIP` for His, protonation states for Asp/Glu, and amide orientations for Asn/Gln. Hover over a cell to read its full text.
4. Use the 3D viewer to inspect the local environment. Drag with the left button to rotate, drag with the middle button to pan, and scroll to zoom. Turn on distance measurement and select two atoms.
5. Confirm a candidate to save your decision and advance to the next site. Export an Excel workbook with columns for site, residue, status, and decision, on a separate sheet for each chain. You can reopen a result folder to resume review.

The desktop app saves candidate PNGs with the offline 3D renderer by default. It chooses a viewpoint that aims to keep the target and key contacts visible. If WebGL is unavailable, it uses a static preview. Turn off **Export rendered images** in Settings to shorten processing time and continue inspecting structures in 3D. Metals such as Mo, Rh, and Os are preserved in the local PDBs and marked in 3D. For metals without element-specific distance rules, automatic His assignments use explicit `LINK` or `CONECT` records and flag unconfirmed close contacts for review.

Automatic assignments use local structure rules and are intended to support model preparation. The pKa column contains optional PROPKA estimates for His, Asp, and Glu. A dash indicates a missing assignment or estimate.

## Command line

The Windows ZIP includes `hischeck-cli.exe`. This example screens His in chain A and writes all three candidates and preview images where their geometry permits. Use a new output directory.

```powershell
.\hischeck-cli.exe .\input.pdb --chain A --residue HIS --all-states --render -o .\output
```

Add `--render` creates static PNGs for headless use. `--progress-json` emits JSON progress lines prefixed with `hischeck-progress `. `--propka` for pKa estimates or repeat `--residue` to include other residue types. `summary.csv` lists the assignments. Each site directory has an `evidence.json` file and gains a `review.json` file when you confirm a selection.

See the [developer guide](docs/DEVELOPMENT.md) for setup, testing, translations, and builds. The original source is under the [MIT License](LICENSE). See the [third-party notices](packaging/windows/THIRD-PARTY.txt) for bundled components.
