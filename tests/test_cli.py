import csv
import io
import json
import math
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from PIL import Image, ImageColor
import gemmi

from hischeck.__main__ import main
from hischeck.analysis import H_NAMES
from hischeck.compare import manual_labels
from hischeck.preview import CARBOXYLATE
from hischeck.review import cancel_review, confirm_review, ensure_variant, export_workbook, read_reviews, review_exports
from hischeck.structure import RING, Structure
from hischeck.view_bonds import display_bonds


def atom_line(serial, name, residue, chain, number, xyz, element, *, het=False, alt=" ", occ=1.0, insertion=" "):
    x, y, z = xyz
    field = name.ljust(4) if len(element) == 2 or len(name) == 4 else " " + name.ljust(3)
    return (f"{'HETATM' if het else 'ATOM  '}{serial:5d} {field}{alt}{residue:>3} {chain}{number:4d}{insertion}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{20.0:6.2f}          {element:>2}")


def radial(name, distance, x_shift=0.0):
    theta = math.tau * RING.index(name) / len(RING)
    return distance * math.cos(theta) + x_shift, distance * math.sin(theta), 0.0


class HischeckTest(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.input = self.root / "input.pdb"

    def fixture(self, nearby=(), *, missing="", alternate=False, original_h=False):
        lines = []

        def add(name, residue, chain, number, xyz, element, **kwargs):
            lines.append(atom_line(len(lines) + 1, name, residue, chain, number, xyz, element, **kwargs))

        for alt, occ, shift in [("A", 0.85, 0), ("B", 0.15, 10)] if alternate else [(" ", 1, 0)]:
            for name in RING:
                if name != missing:
                    add(name, "HIS", "A", 5, radial(name, 1.18, shift), name[0], alt=alt, occ=occ, insertion="A")
        if original_h:
            add("HE2", "HIS", "A", 5, radial("NE2", 2.19), "H", insertion="A")
            add("HA", "HIS", "A", 5, (1.0, 0, 4.0), "H", insertion="A")
        for number, (nitrogen, kind) in enumerate(nearby, 10):
            if kind == "acceptor":
                add("OD1", "ASN", "B", number, radial(nitrogen, 4.10), "O")
                # This atom falls outside the export radius but is part of the selected neighbor residue.
                add("CA", "ASN", "B", number, radial(nitrogen, 12.0), "C")
            elif kind == "donor":
                add("NZ", "LYS", "B", number, radial(nitrogen, 4.15), "N")
                add("HZ1", "LYS", "B", number, radial(nitrogen, 3.15), "H")
            elif kind in {"water", "water_with_h"}:
                ox = radial(nitrogen, 4.15)
                add("O", "HOH", "B", number, ox, "O", het=True)
                if kind == "water_with_h":
                    add("H1", "HOH", "B", number, (ox[0], ox[1], 0.97), "H", het=True)
                    add("H2", "HOH", "B", number, (*radial(nitrogen, 5.08)[:2], -0.24), "H", het=True)
            elif kind == "hydroxyl":
                add("OG1", "THR", "B", number, radial(nitrogen, 4.15), "O")
            elif kind == "amine":
                add("NZ", "LYN", "B", number, radial(nitrogen, 4.15), "N")
            elif kind == "sulfur":
                add("SD", "MET", "B", number, radial(nitrogen, 4.78), "S")
            elif kind == "thiolate":
                add("SG", "CYM", "B", number, (*radial(nitrogen, 4.18)[:2], 1.1), "S")
            elif kind == "metal":
                add("ZN", "ZN", "B", number, radial(nitrogen, 3.23), "ZN", het=True)
            elif kind in {"asp", "glu", "glu_at_limit", "glu_beyond_limit", "protonated_asp", "ash", "terminus"}:
                residue = {"asp": "ASP", "protonated_asp": "ASP", "ash": "ASH", "terminus": "ALA"}.get(kind, "GLU")
                first, second, carbon = {"GLU": ("OE1", "OE2", "CD"), "ALA": ("O", "OXT", "C")}.get(
                    residue, ("OD1", "OD2", "CG"))
                if kind in {"glu_at_limit", "glu_beyond_limit"}:
                    x, y, _ = radial(nitrogen, 1.18)
                    z = 5.0 if kind == "glu_at_limit" else 5.01
                    add(first, residue, "B", number, (x, y, z), "O")
                    add(second, residue, "B", number, (x, y, z + 2.23), "O")
                    add(carbon, residue, "B", number, (x + 0.58, y, z + 1.115), "C")
                else:
                    radius = 5.45 if kind == "glu" else 4.10
                    x, y, _ = radial(nitrogen, radius)
                    add(first, residue, "B", number, (x, y, 0), "O")
                    add(second, residue, "B", number, (x, y, 2.23), "O")
                    add(carbon, residue, "B", number, (*radial(nitrogen, radius + 0.58)[:2], 1.115), "C")
                if kind == "protonated_asp":
                    add("HD2", residue, "B", number, (x, y, 3.20), "H")
        self.input.write_text("\n".join([*lines, "END"]) + "\n", encoding="ascii")

    def run_cli(self, *options, name="out"):
        output = self.root / name
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.input), "-o", str(output), "--chain", "A", *options]), 0)
        folder = output / "A_5A_HIS"
        with (output / "summary.csv").open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["site"], "A:5A")
        return json.loads((folder / "evidence.json").read_text(encoding="utf-8")), folder

    def test_directional_contacts_and_hydrogen_export(self):
        for state, protonated, unprotonated in [("HID", "ND1", "NE2"), ("HIE", "NE2", "ND1")]:
            with self.subTest(state=state):
                self.fixture([(protonated, "acceptor"), (unprotonated, "donor")], original_h=True)
                data, folder = self.run_cli("--radius", "5", name=state)
                self.assertEqual(data["recommendation"], state)
                model = gemmi.read_pdb(str(folder / "site.pdb"))[0]
                his = model["A"][0]
                atoms = {atom.name: atom for atom in his}
                self.assertEqual(len(atoms), len(his))
                self.assertEqual(his.name, state)
                self.assertIn("HA", atoms)
                self.assertNotIn(H_NAMES[unprotonated], atoms)
                for parent, length in [(protonated, 1.01), ("CE1", 1.08), ("CD2", 1.08)]:
                    self.assertAlmostEqual(atoms[parent].pos.dist(atoms[H_NAMES[parent]].pos), length, places=2)
                    self.assertEqual(atoms[H_NAMES[parent]].pos.z, 0.0)
                self.assertIn("CA", [atom.name for atom in model["B"][0]])

    def test_two_acceptors_and_explicit_pH(self):
        self.fixture([("ND1", "acceptor"), ("NE2", "acceptor")])
        for name, options in [("default", []), ("acidic", ["--ph", "5"])]:
            with self.subTest(name=name):
                data, folder = self.run_cli(*options, name=name)
                self.assertEqual(data["recommendation"], "HIP")
                self.assertEqual(data["reasons"], ["two_primary_acceptors"])
                self.assertTrue(data["hip_evidence"]["two_primary_acceptors"])
                self.assertTrue((folder / "site.pdb").exists())
                hip = next(candidate for candidate in data["candidates"] if candidate["state"] == "HIP")
                self.assertEqual(hip["pH_prior"], 0.45 if options else 0.0)
                self.assertEqual({(c["nitrogen"], c["neighbor"]) for c in hip["supports"]},
                                 {("ND1", "B:10 ASN OD1"), ("NE2", "B:11 ASN OD1")})

    def test_carboxylate_hip_recommends_and_exports_both_protons(self):
        for kind, oxygen in [("asp", "ASP OD1"), ("terminus", "ALA O")]:
            with self.subTest(kind=kind):
                self.fixture([("ND1", kind)])
                data, folder = self.run_cli("--render", name=kind)
                self.assertEqual(data["recommendation"], "HIP")
                self.assertEqual(data["top_candidate"], "HIP")
                self.assertEqual(data["reasons"], ["nearby_carboxylate"])
                contact = data["carboxylate_contacts"][0]
                self.assertEqual(contact["atom"], f"B:10 {oxygen}")
                self.assertEqual(contact["nitrogen"], "ND1")
                self.assertAlmostEqual(contact["distance_A"], 2.92, places=2)
                his = gemmi.read_pdb(str(folder / "site.pdb"))[0]["A"][0]
                self.assertEqual(his.name, "HIP")
                self.assertTrue({"HD1", "HE2"} <= {atom.name for atom in his})
                with Image.open(folder / "site.png") as image:
                    colors = image.getcolors(image.width * image.height)
                    self.assertIn(ImageColor.getrgb(CARBOXYLATE), {color for _, color in colors})

    def test_carboxylate_distance_respects_search_and_does_not_treat_acids_as_charged(self):
        self.fixture([("NE2", "glu")])
        for name, options, expected in [("default", [], True), ("smaller", ["--carboxylate-max", "4"], False),
                                        ("capped", ["--carboxylate-max", "5", "--search-radius", "4.2"], False)]:
            with self.subTest(name=name):
                data, _ = self.run_cli(*options, name=name)
                self.assertEqual("nearby_carboxylate" in data["reasons"], expected)
                self.assertEqual(bool(data["carboxylate_contacts"]), expected)
        for kind, expected in [("glu_at_limit", True), ("glu_beyond_limit", False)]:
            with self.subTest(kind=kind):
                self.fixture([("NE2", kind)])
                data, _ = self.run_cli("--search-radius", "6", name=kind)
                self.assertEqual(bool(data["carboxylate_contacts"]), expected)
                self.assertEqual(data["recommendation"], "HIP" if expected else None)
                if expected:
                    self.assertEqual(data["carboxylate_contacts"][0]["distance_A"], 5.0)
        for kind in ("protonated_asp", "ash"):
            with self.subTest(kind=kind):
                self.fixture([("ND1", kind)])
                data, _ = self.run_cli(name=kind)
                self.assertNotIn("nearby_carboxylate", data["reasons"])
                self.assertFalse(data["carboxylate_contacts"])

    def test_metal_coordination_and_severe_clashes_take_precedence_over_carboxylate(self):
        for partner in ("metal", "donor"):
            with self.subTest(partner=partner):
                self.fixture([("ND1", "asp"), ("NE2", partner)])
                data, _ = self.run_cli(name=partner)
                hip = next(candidate for candidate in data["candidates"] if candidate["state"] == "HIP")
                self.assertTrue(data["carboxylate_contacts"])
                self.assertNotIn("nearby_carboxylate", data["reasons"])
                self.assertEqual(data["recommendation"], "HID")
                if partner == "metal":
                    self.assertEqual(data["reasons"], ["metal_coordination"])
                    self.assertFalse(hip["allowed"])
                else:
                    self.assertTrue(any(clash["severe"] for clash in hip["clashes"]))

    def test_water_without_known_hydrogens_is_weak_evidence(self):
        self.fixture([("ND1", "water")])
        data, folder = self.run_cli()
        self.assertEqual(data["status"], "review")
        self.assertIsNone(data["recommendation"])
        self.assertTrue((folder / "site_HID.pdb").exists())
        self.assertTrue((folder / "site_HIE.pdb").exists())

    def test_secondary_contacts_do_not_admit_hip_even_at_low_pH(self):
        for first, second in [("water", "water_with_h"), ("acceptor", "water"), ("acceptor", "hydroxyl"),
                              ("acceptor", "thiolate"), ("acceptor", "amine")]:
            with self.subTest(first=first, second=second):
                self.fixture([("ND1", first), ("NE2", second)])
                data, folder = self.run_cli("--ph", "3", name=f"{first}_{second}")
                hip = next(c for c in data["candidates"] if c["state"] == "HIP")
                self.assertNotEqual(data["recommendation"], "HIP")
                self.assertFalse(hip["allowed"])
                self.assertEqual(hip["exclusion_reason"], "insufficient_HIP_acceptors")
                self.assertFalse(any(e["state"] == "HIP" for e in data["exports"]))
                self.assertFalse((folder / "site_HIP.pdb").exists())

    def test_oxygen_takes_precedence_over_sulfur_and_amine_prefers_donation(self):
        self.fixture([("ND1", "thiolate"), ("NE2", "amine")])
        prior, _ = self.run_cli(name="no_oxygen")
        sulfur = next(c for c in prior["contacts"]["ND1"] if c["element"] == "S")
        self.assertGreater(sulfur["supports_His_as_donor"], 0)
        amine = next(c for c in prior["contacts"]["NE2"] if c["role"] == "amine")
        self.assertGreater(amine["supports_His_as_acceptor"], amine["supports_His_as_donor"])
        self.assertEqual(amine["priority"], "secondary")

        self.fixture([("ND1", "thiolate"), ("ND1", "acceptor"), ("NE2", "amine")])
        data, _ = self.run_cli(name="with_oxygen")
        sulfur = next(c for c in data["contacts"]["ND1"] if c["element"] == "S")
        self.assertTrue(sulfur["deferred_to_oxygen"])
        self.assertEqual(sulfur["supports_His_as_donor"], 0)
        self.assertEqual(data["recommendation"], "HID")
        hid = next(c for c in data["candidates"] if c["state"] == "HID")
        self.assertIn("B:11 ASN OD1", {c["neighbor"] for c in hid["supports"]})

    def test_ordinary_and_double_metal_coordination(self):
        self.fixture([("NE2", "metal")])
        data, _ = self.run_cli()
        self.assertEqual(data["recommendation"], "HID")
        self.assertTrue(data["metal_contacts"][0]["confirmed"])
        self.fixture([("NE2", "metal"), ("ND1", "metal")])
        data, folder = self.run_cli(name="double")
        self.assertEqual(data["status"], "unresolved")
        self.assertIn("both_nitrogens_coordinate", data["reasons"])
        self.assertTrue((folder / "site_input.pdb").exists())

    def test_contact_limits_change_evidence_and_flag_controls_pngs(self):
        self.fixture([("ND1", "acceptor"), ("NE2", "sulfur")])
        baseline, default_output = self.run_cli()
        self.assertGreater(baseline["contacts"]["ND1"][0]["supports_His_as_donor"], 0)
        self.assertGreater(baseline["contacts"]["NE2"][0]["supports_His_as_donor"], 0)
        self.assertFalse(list(default_output.glob("*.png")))

        adjusted, folder = self.run_cli("--render", "--hbond-max", "2.8", "--s-bond-max", "3.5", name="adjusted")
        self.assertEqual(adjusted["contacts"]["ND1"][0]["supports_His_as_donor"], 0)
        self.assertEqual(adjusted["contacts"]["NE2"][0]["supports_His_as_donor"], 0)
        self.assertEqual(adjusted["status"], "review")
        for exported in adjusted["exports"]:
            self.assertEqual(Path(exported["image"]).stem, Path(exported["pdb"]).stem)
            with Image.open(folder / exported["image"]) as image:
                self.assertEqual(image.format, "PNG")
                self.assertEqual(image.size, (960, 660))
                image.verify()

    def test_all_states_preserves_metal_decision_and_builds_each_variant(self):
        self.fixture([("NE2", "metal")], original_h=True)
        default, _ = self.run_cli(name="automatic")
        data, folder = self.run_cli("--all-states", "--render", name="reviewable")
        self.assertEqual(data["automatic_states"], ["HID"])
        self.assertEqual(data["exports"], [{"pdb": f"site_{state}.pdb", "state": state, "image": f"site_{state}.png"}
                                          for state in ("HID", "HIE", "HIP")])
        for field in ("status", "recommendation", "top_candidate", "candidates", "reasons"):
            self.assertEqual(data[field], default[field])
        for export, protons in zip(data["exports"], ({"HD1"}, {"HE2"}, {"HD1", "HE2"})):
            target = gemmi.read_pdb(str(folder / export["pdb"]))[0]["A"][0]
            self.assertEqual(target.name, export["state"])
            self.assertEqual({atom.name for atom in target if atom.element.is_hydrogen}, protons | {"HE1", "HD2", "HA"})
            self.assertTrue((folder / export["image"]).is_file())
        self.fixture([("NE2", "metal"), ("ND1", "metal")])
        data, _ = self.run_cli("--all-states", name="double_all")
        self.assertEqual(data["status"], "unresolved")
        self.assertEqual([export["state"] for export in data["exports"]], ["HID", "HIE", "HIP", None])
        self.fixture(missing="NE2")
        data, _ = self.run_cli("--all-states", name="missing_all")
        self.assertEqual(data["exports"], [{"pdb": "site_input.pdb", "state": None}])

    def test_lazy_variant_and_manual_workbook(self):
        from openpyxl import load_workbook

        self.fixture([("NE2", "metal")])
        evidence, folder = self.run_cli()
        original_evidence = (folder / "evidence.json").read_bytes()
        original_pdb = (folder / "site.pdb").read_bytes()
        variant = next(export for export in review_exports(evidence) if export["state"] == "HIE")
        self.assertFalse((folder / variant["pdb"]).exists())
        pdb = ensure_variant(folder, evidence, variant)
        self.assertEqual((folder / "evidence.json").read_bytes(), original_evidence)
        self.assertEqual((folder / "site.pdb").read_bytes(), original_pdb)
        self.assertEqual(gemmi.read_pdb(str(pdb))[0]["A"][0].name, "HIE")
        with (folder.parent / "summary.csv").open(encoding="utf-8-sig", newline="") as file:
            rows = list(csv.DictReader(file))
        source = rows[0]
        rows.extend([dict(source, site="A:6", directory="A_6_HIS", status="review", recommendation="",
                          top_candidate="HID"),
                     dict(source, site="B:5", chain="B", directory="B_5_HIS", status="unresolved",
                          recommendation="", top_candidate="")])
        confirm_review(folder / "site.pdb", "HID")
        confirm_review(pdb, "HIE")
        confirmed = read_reviews(folder.parent, rows)
        self.assertEqual(confirmed, {folder.name: "HIE"})
        self.assertEqual(json.loads((folder / "review.json").read_text()), {"state": "HIE", "pdb": "site_HIE.pdb"})
        workbook = self.root / "判断.xlsx"
        export_workbook(workbook, rows, confirmed)
        result = load_workbook(workbook)
        self.addCleanup(result.close)
        self.assertEqual(result.sheetnames, ["链 A", "链 B"])
        header = ("位点", "残基", "状态", "判断")
        self.assertEqual(list(result["链 A"].values), [header, ("A:5A", "HIS", "已核对", "HIE"),
                                                      ("A:6", "HIS", "待复核", "HID")])
        self.assertEqual(list(result["链 B"].values), [header, ("B:5", "HIS", "未解析", "—")])
        self.assertEqual(result["链 A"].auto_filter.ref, "A1:D3")
        self.assertEqual(result["链 A"].freeze_panes, "A2")
        self.assertEqual(manual_labels(workbook, "链 A", chain="A"), [
            {"residue_number": 5, "insertion_code": "A", "manual_state": "HIE", "manual_note": "", "xlsx_row": 2},
        ])

    def test_aspartate_default_uses_side_chain_distance_and_removes_only_acid_proton(self):
        atoms = [
            ("CG", "ASH", "A", 5, (0, 0, 0), "C"),
            ("OD1", "ASH", "A", 5, (0.6, 1.1, 0), "O"),
            ("OD2", "ASH", "A", 5, (0.6, -1.1, 0), "O"),
            ("DD2", "ASH", "A", 5, (1.6, -1.1, 0), "D"),
            ("CA", "ASH", "A", 5, (-3, 0, 0), "C"),
            ("HA", "ASH", "A", 5, (-3, 0, 1.09), "H"),
            ("O", "HOH", "B", 6, (0.6, 1.1, 2.0), "O"),
            ("H1", "HOH", "B", 6, (0.6, 1.1, 2.95), "H"),
            ("O", "HOH", "B", 7, (-3, 0, 1.5), "O"),
        ]
        lines = [atom_line(i, *args, het=args[1] == "HOH", insertion="A" if args[2] == "A" else " ")
                 for i, args in enumerate(atoms, 1)]
        self.input.write_text("\n".join([*lines, "CONECT    3    4", "CONECT    7    8", "END"]) + "\n")
        output = self.root / "aspartate"
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.input), "--residue", "asp", "--chain", "A", "--radius", "2",
                                   "--render", "-o", str(output)]), 0)
        folder = output / "A_5A_ASP"
        local = gemmi.read_pdb(str(folder / "site.pdb"))
        target = local[0]["A"][0]
        self.assertEqual(target.name, "ASP")
        self.assertEqual({atom.name for atom in target}, {"CG", "OD1", "OD2", "CA", "HA"})
        self.assertEqual([res.seqid.num for res in local[0]["B"]], [6])
        self.assertEqual({atom.name for atom in local[0]["B"][0]}, {"O", "H1"})
        serials = {atom.serial: atom.name for chain in local[0] for res in chain for atom in res}
        self.assertEqual({frozenset((serials[a], serials[b])) for a, partners in local.conect_map.items()
                          for b in partners}, {frozenset(("O", "H1"))})
        evidence = json.loads((folder / "evidence.json").read_text())
        self.assertEqual((evidence["status"], evidence["recommendation"]), ("default", "deprotonated"))
        self.assertEqual(evidence["removed_sidechain_hydrogens"], ["DD2"])
        with Image.open(folder / "site.png") as image:
            image.verify()

    def test_mixed_review_keeps_amide_coordinates_and_exports_pending_rows(self):
        from openpyxl import load_workbook

        self.fixture()
        lines = self.input.read_text().splitlines()[:-1]
        for number, name, atoms in ((8, "ASN", ("CG", "OD1", "ND2")), (9, "GLN", ("CD", "OE1", "NE2"))):
            for atom, point in zip(atoms, ((5, 0, 0), (5.6, 1.1, 0), (5.6, -1.1, 0))):
                lines.append(atom_line(len(lines) + 1, atom, name, "A", number, point, atom[0]))
        lines.extend([atom_line(len(lines) + 1, "CD", "GLU", "A", 10, (0, 6, 0), "C"), "END"])
        self.input.write_text("\n".join(lines) + "\n")
        output = self.root / "mixed"
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.input), "--all-states", "--residue", "HIS", "--residue", "ASN",
                                   "--residue", "GLN", "--residue", "GLU", "-o", str(output)]), 0)
        with (output / "summary.csv").open(encoding="utf-8-sig", newline="") as file:
            rows = list(csv.DictReader(file))
        self.assertEqual([row["input_residue"] for row in rows], ["HIS", "ASN", "GLN", "GLU"])
        self.assertEqual([row["recommendation"] for row in rows[1:]], ["no_flip", "no_flip", "deprotonated"])
        for row in rows[1:3]:
            folder = output / row["directory"]
            structure = gemmi.read_pdb(str(folder / "site.pdb"))
            target = next(residue for residue in structure[0]["A"] if residue.seqid.num == int(row["residue_number"]))
            coordinates = {atom.name: (atom.pos.x, atom.pos.y, atom.pos.z) for atom in target}
            polar_n = "ND2" if target.name == "ASN" else "NE2"
            self.assertEqual(coordinates[polar_n], (5.6, -1.1, 0.0))
            self.assertEqual(len(coordinates), 5)
            self.assertEqual({name for name in coordinates if name.startswith("H")},
                             {"HD21", "HD22"} if target.name == "ASN" else {"HE21", "HE22"})
        asp_missing = json.loads((output / "A_10_GLU" / "evidence.json").read_text())
        self.assertEqual(asp_missing["missing_functional_group_atoms"], ["OE1", "OE2"])
        self.assertEqual([entry["state"] for entry in review_exports(asp_missing)], ["deprotonated"])
        confirm_review(output / "A_8_ASN" / "site.pdb", "no_flip", kind="ASN")
        with self.assertRaises(ValueError):
            confirm_review(output / "A_8_ASN" / "site.pdb", "HIE", kind="ASN")
        workbook = self.root / "混合判断.xlsx"
        export_workbook(workbook, rows, read_reviews(output, rows))
        result = load_workbook(workbook, read_only=True)
        try:
            self.assertEqual(list(result.active.values), [
                ("位点", "残基", "状态", "判断"), ("A:5A", "HIS", "待复核", "—"),
                ("A:8", "ASN", "已核对", "不翻转"), ("A:9", "GLN", "推荐", "不翻转"),
                ("A:10", "GLU", "推荐", "去质子化"),
            ])
        finally:
            result.close()
        folder = output / "A_8_ASN"
        original_pdb = (folder / "site.pdb").read_bytes()
        cancel_review(folder)
        cancel_review(folder)
        self.assertEqual(read_reviews(output, rows), {})
        export_workbook(workbook, rows, read_reviews(output, rows))
        result = load_workbook(workbook, read_only=True)
        try:
            self.assertEqual((result.active["C3"].value, result.active["D3"].value), ("推荐", "不翻转"))
        finally:
            result.close()
        self.assertEqual((folder / "site.pdb").read_bytes(), original_pdb)

    def sidechain_run(self, kind, partners=(), *, buried=True, name, first_oxygen_long=False, conect=(), render=False):
        names = {"ASP": ("CG", "OD1", "OD2"), "GLU": ("CD", "OE1", "OE2"),
                 "ASN": ("CG", "OD1", "ND2"), "GLN": ("CD", "OE1", "NE2")}[kind]
        first = (0.675, 1.169, 0) if first_oxygen_long else (0.620, 1.074, 0)
        other = (0.620, -1.074, 0) if first_oxygen_long else (0.665, -1.152, 0)
        atoms = [(atom, kind, "A", 5, point, atom[0]) for atom, point in zip(names, ((0, 0, 0), first, other))]
        atoms.extend(partners)
        if buried:
            edges = range(-6, 7, 2)
            shell = [(x, y, z) for x in edges for y in edges for z in edges if max(abs(x), abs(y), abs(z)) == 6]
            atoms.extend(("CB", "ALA", "C", number, point, "C") for number, point in enumerate(shell, 1))
        lines = [atom_line(i, *args, het=args[1] in {"HOH", "ZN"}) for i, args in enumerate(atoms, 1)]
        self.input.write_text("\n".join([*lines, *conect, "END"]) + "\n")
        output = self.root / name
        options = ["--render"] if render else []
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.input), "--residue", kind, "--chain", "A", "-o", str(output),
                                   "--all-states", *options]), 0)
        folder = output / f"A_5_{kind}"
        return json.loads((folder / "evidence.json").read_text(encoding="utf-8")), folder

    def test_buried_acid_exports_both_states_and_places_a_cis_proton_on_single_bond(self):
        from openpyxl import load_workbook

        for kind, oxygen, proton, residue in (("ASP", "OD2", "HD2", "ASH"), ("GLU", "OE2", "HE2", "GLH")):
            with self.subTest(kind=kind):
                data, folder = self.sidechain_run(kind, name=kind, first_oxygen_long=True, render=True,
                                                 conect=("CONECT    1    2    3    3",) if kind == "ASP" else ())
                self.assertEqual(data["recommendation"], "protonated")
                with (folder.parent / "summary.csv").open(encoding="utf-8-sig", newline="") as file:
                    rows = list(csv.DictReader(file))
                workbook = self.root / f"{kind}.xlsx"
                export_workbook(workbook, rows, {})
                result = load_workbook(workbook, read_only=True)
                try:
                    self.assertEqual(list(result.active.values), [
                        ("位点", "残基", "状态", "判断"), ("A:5", kind, "待复核", "质子化"),
                    ])
                finally:
                    result.close()
                self.assertFalse(data["surface_exposure"]["exposed"])
                self.assertEqual(data["protonation_oxygen"]["input_atom"], oxygen[:-1] + "1")
                self.assertEqual(data["protonation_oxygen"]["basis"],
                                 "conect_double_bond" if kind == "ASP" else "longer_C_O_bond")
                self.assertEqual([item["state"] for item in data["exports"]], ["deprotonated", "protonated"])
                for item in data["exports"]:
                    with Image.open(folder / item["image"]) as image:
                        image.verify()
                base = gemmi.read_pdb(str(folder / "site.pdb"))[0]["A"][0]
                model = gemmi.read_pdb(str(folder / "site_protonated.pdb"))
                acid = model[0]["A"][0]
                self.assertEqual((base.name, acid.name), (kind, residue))
                donor, acceptor, hydrogen = (acid[name][0] for name in (oxygen, oxygen[:-1] + "1", proton))
                self.assertLess(donor.pos.dist(base[oxygen[:-1] + "1"][0].pos), 1e-8)
                self.assertLess(acceptor.pos.dist(base[oxygen][0].pos), 1e-8)
                carbon = acid["CG" if kind == "ASP" else "CD"][0]
                self.assertAlmostEqual(donor.pos.dist(hydrogen.pos), 0.96, places=3)
                self.assertAlmostEqual(math.degrees(gemmi.calculate_angle(carbon.pos, donor.pos, hydrogen.pos)),
                                       108.5, delta=0.1)
                self.assertAlmostEqual(gemmi.calculate_dihedral(acceptor.pos, carbon.pos, donor.pos, hydrogen.pos),
                                       0, delta=0.002)
                if kind == "ASP":
                    self.assertEqual(list(model.conect_map[carbon.serial]).count(acceptor.serial), 2)

    def test_default_acid_keeps_water_donor_cations_and_surface_and_generates_lazily(self):
        conditions = [
            ("water", ("O", "HOH", "B", 7, (3.820, 1.074, 0), "O"), "nearby_crystal_water"),
            ("nh", ("N", "ALA", "B", 7, (3.820, 1.074, 0), "N"), "nearby_hbond_donor"),
            ("oh", ("OG", "SER", "B", 7, (3.820, 1.074, 0), "O"), "nearby_hbond_donor"),
            ("lys", ("NZ", "LYS", "B", 7, (4.420, 1.074, 0), "N"), "nearby_positive_lysine"),
            ("arg", ("NH1", "ARG", "B", 7, (4.420, 1.074, 0), "N"), "nearby_positive_arginine"),
            ("metal", ("ZN", "ZN", "B", 7, (2.720, 1.074, 0), "ZN"), "sidechain_metal_coordination"),
            ("surface", None, "sidechain_exposed"),
        ]
        for case, partner, reason in conditions:
            with self.subTest(case=case):
                data, folder = self.sidechain_run("GLU", [partner] if partner else [], buried=partner is not None,
                                                 name=case, first_oxygen_long=True)
                self.assertEqual(data["recommendation"], "deprotonated")
                self.assertIn(reason, data["reasons"])
                self.assertEqual([p.name for p in folder.glob("*.pdb")], ["site.pdb"])
                if case == "surface":
                    variant = next(item for item in review_exports(data) if item["state"] == "protonated")
                    evidence_bytes = (folder / "evidence.json").read_bytes()
                    base_bytes = (folder / "site.pdb").read_bytes()
                    path = ensure_variant(folder, data, variant)
                    self.assertEqual(ensure_variant(folder, data, variant), path)
                    molecule = gemmi.read_pdb(str(path))[0]["A"][0]
                    self.assertEqual(molecule.name, "GLH")
                    self.assertLess(molecule["OE2"][0].pos.dist(gemmi.Position(0.675, 1.169, 0)), 1e-8)
                    self.assertEqual((folder / "evidence.json").read_bytes(), evidence_bytes)
                    self.assertEqual((folder / "site.pdb").read_bytes(), base_bytes)
                    confirm_review(path, "protonated", kind="GLU")
                    self.assertEqual(read_reviews(folder.parent, [{"input_residue": "GLU", "directory": folder.name}]),
                                     {folder.name: "protonated"})

    def test_short_crystal_water_contact_does_not_block_water_probe(self):
        partner = ("O", "HOH", "B", 7, (0.665, -3.462, 0), "O")
        dry, _ = self.sidechain_run("ASP", buried=False, name="dry_surface")
        wet, _ = self.sidechain_run("ASP", [partner], buried=False, name="wet_surface")
        self.assertTrue(dry["surface_exposure"]["exposed"])
        self.assertEqual(dry["surface_exposure"], wet["surface_exposure"])
        data, folder = self.sidechain_run("ASP", [partner], name="wet_buried")
        self.assertFalse(data["surface_exposure"]["exposed"])
        self.assertEqual(data["recommendation"], "deprotonated")
        self.assertEqual(data["reasons"], ["nearby_crystal_water"])
        support, = data["candidates"][0]["supports"]
        self.assertEqual((support["center"], support["interaction"], support["distance_A"]), ("OD2", "water", 2.31))
        self.assertFalse((folder / "site_protonated.pdb").exists())

    def test_display_bonds_keep_peptides_and_use_known_bond_orders(self):
        def orders(residue: str, separator: str):
            specs = [
                ("C", "LEU", "A", 4, (0, 0, 0), "C"), ("O", "LEU", "A", 4, (0, 1.23, 0), "O"),
                ("N", residue, "A", 5, (1.33, 0, 0), "N"), ("C", residue, "A", 5, (7, 0, 0), "C"),
                ("CG", residue, "A", 5, (4, 0, 0), "C"), ("OD1", residue, "A", 5, (4.6, 1.04, 0), "O"),
                ("OD2", residue, "A", 5, (4.65, -1.13, 0), "O"), ("N", "HIE", "A", 6, (8.33, 0, 0), "N"),
            ]
            specs.extend((name, "HIE", "A", 6, radial(name, 1.13, 11), name[0]) for name in RING)
            specs.extend([("C1", "LIG", "B", 1, (0, 3, 0), "C"), ("C2", "LIG", "B", 1, (1.3, 3, 0), "C")])
            lines = [atom_line(serial, *args, het=args[1] == "LIG") for serial, args in enumerate(specs, 1)]
            text = "\n".join([*lines[:2], *([separator] if separator else []), *lines[2:], "CONECT   14   15   15", "END"])
            self.input.write_text(text)
            return {(first, second): order for first, second, order in display_bonds(Structure(self.input), text)}

        acid = orders("ASH", "")
        for pair, order in { (1, 2): 2, (1, 3): 1, (4, 8): 1, (5, 6): 2, (5, 7): 1,
                             (9, 13): 2, (10, 11): 2, (14, 15): 2 }.items():
            self.assertEqual(acid[pair], order)
        self.assertNotIn((1, 3), orders("ASH", "TER"))
        charged = orders("ASP", "")
        self.assertEqual((charged[5, 6], charged[5, 7]), (1.5, 1.5))

    def test_amide_uses_in_plane_contacts_and_swaps_coordinates_for_the_alternative(self):
        acceptor = ("O", "ALA", "B", 6, (3.620, 1.074, 0), "O")
        contrary = ("O", "HOH", "B", 7, (2.169, 3.759, 0), "O")
        perpendicular = ("O", "ALA", "B", 6, (0.620, 1.074, 3.0), "O")
        amine = ("ND2", "ASN", "B", 6, (2.165, -3.750, 0), "N")
        cases = [("flip", [acceptor], True, "flip"), ("swap_acceptor", [amine], True, "flip"),
                 ("keep", [acceptor, contrary], True, "no_flip"), ("out_of_plane", [perpendicular], True, "no_flip"),
                 ("surface_amide", [acceptor], False, "no_flip"), ("empty", [], True, "no_flip")]
        for case, partners, buried, recommendation in cases:
            with self.subTest(case=case):
                data, folder = self.sidechain_run("ASN", partners, buried=buried, name=case)
                self.assertEqual(data["recommendation"], recommendation)
                alternative = folder / "site_flip.pdb"
                self.assertEqual(alternative.is_file(), recommendation == "flip")
                if case == "keep":
                    self.assertIn("original_amide_hbond", data["reasons"])
                if case == "out_of_plane":
                    self.assertTrue(any(item["plane_angle_deg"] == 90 and not item["supports"]
                                        for item in data["contacts"]["flip"]))
                if case == "surface_amide":
                    variant = next(item for item in review_exports(data) if item["state"] == "flip")
                    self.assertEqual(ensure_variant(folder, data, variant), alternative)
                    confirm_review(alternative, "flip", kind="ASN")
                if alternative.is_file():
                    base = gemmi.read_pdb(str(folder / "site.pdb"))[0]["A"][0]
                    flipped = gemmi.read_pdb(str(alternative))[0]["A"][0]
                    self.assertLess(flipped["OD1"][0].pos.dist(base["ND2"][0].pos), 1e-8)
                    self.assertLess(flipped["ND2"][0].pos.dist(base["OD1"][0].pos), 1e-8)
                    carbon, oxygen, nitrogen, h1, h2 = (flipped[name][0].pos for name in ("CG", "OD1", "ND2", "HD21", "HD22"))
                    for first, last in ((carbon, h1), (carbon, h2), (h1, h2)):
                        self.assertAlmostEqual(math.degrees(gemmi.calculate_angle(first, nitrogen, last)), 120, delta=0.1)
                    for h in (h1, h2):
                        self.assertAlmostEqual(nitrogen.dist(h), 1.0, places=3)
                        torsion = gemmi.calculate_dihedral(oxygen, carbon, nitrogen, h)
                        self.assertLess(abs(math.sin(torsion)), 0.002)

    def test_altloc_and_missing_atoms(self):
        self.fixture([("ND1", "acceptor")], alternate=True)
        data, folder = self.run_cli()
        self.assertIn("target_alternate_conformation", data["reasons"])
        self.assertEqual(data["target_altloc"], "A")
        model = gemmi.read_pdb(str(folder / "site_HID.pdb"))[0]
        self.assertTrue(all(atom.pos.x < 3 for atom in model["A"][0]))
        self.fixture(missing="NE2")
        data, folder = self.run_cli(name="missing")
        self.assertEqual(data["missing_ring_atoms"], ["NE2"])
        self.assertTrue((folder / "site_input.pdb").exists())


if __name__ == "__main__":
    unittest.main()
