"""Add protein bond orders and peptide links for the local PDB viewer."""

from itertools import pairwise

import gemmi
import numpy as np

from .chemistry import HIS_ALIAS, PROTEIN
from .structure import CARBOXYLATE_OXYGENS, RING, SiteResidue, Structure

DELOCALIZED = 1.5
PHENYL = ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")
AROMATIC_RINGS = {"PHE": (PHENYL,), "TYR": (PHENYL,),
                  "TRP": (("CG", "CD1", "NE1", "CE2", "CD2"), ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")),
                  "HIS": (RING,), "HIP": (RING,)}
DOUBLE_BONDS = {
    "ASN": (("CG", "OD1"),), "GLN": (("CD", "OE1"),),
    "ASH": (("CG", "OD1"),), "GLH": (("CD", "OE1"),),
    "HID": (("CG", "CD2"), ("CE1", "NE2")),
    "HIE": (("CG", "CD2"), ("ND1", "CE1")),
}


def residue_orders(residue: SiteResidue):
    name = HIS_ALIAS.get(residue.name, residue.name)
    if "O" in residue.atoms and "OXT" in residue.atoms:
        protonated = "HXT" in residue.atoms
        yield "C", "O", 2 if protonated else DELOCALIZED
        yield "C", "OXT", 1 if protonated else DELOCALIZED
    else:
        yield "C", "O", 2
    if name == "ASH":
        yield "CG", "OD2", 1
    elif name == "GLH":
        yield "CD", "OE2", 1
    for first, second in DOUBLE_BONDS.get(name, ()):
        yield first, second, 2
    if name in CARBOXYLATE_OXYGENS:
        carbon = "CG" if name == "ASP" else "CD"
        for oxygen in CARBOXYLATE_OXYGENS[name]:
            yield carbon, oxygen, DELOCALIZED
    for ring in AROMATIC_RINGS.get(name, ()):
        for first, second in pairwise((*ring, ring[0])):
            yield first, second, DELOCALIZED


def peptide_links(pdb_text: str):
    # A second, TER-aware read keeps chain segments distinct in the display topology.
    model = gemmi.read_pdb_string(pdb_text, split_chain_on_ter=True)[0]
    for chain in model:
        for first, second in pairwise(chain):
            if first.name not in PROTEIN or second.name not in PROTEIN:
                continue
            carbon = first.find_atom("C", "*")
            nitrogen = second.find_atom("N", "*")
            if carbon and nitrogen and 1.15 <= carbon.pos.dist(nitrogen.pos) <= 1.65:
                yield tuple(sorted((carbon.serial, nitrogen.serial)))


def display_bonds(structure: Structure, pdb_text: str) -> list[tuple[int, int, float]]:
    """Return only the bond additions and orders the PDB viewer needs to apply."""
    serials = {atom.pdb.serial for atom in structure.atoms}
    bonds = {pair: order for pair, order in structure.bond_orders.items()
             if 1 <= order <= 3 and all(serial in serials for serial in pair)}
    for pair in peptide_links(pdb_text):
        if all(serial in serials for serial in pair):
            bonds.setdefault(pair, 1)
    for residue in structure.residues_by_name.values():
        if residue.name not in PROTEIN:
            continue
        for first, second, order in residue_orders(residue):
            left, right = residue.atoms.get(first), residue.atoms.get(second)
            if left is None or right is None or not 1.05 <= np.linalg.norm(left.position - right.position) <= 1.70:
                continue
            pair = tuple(sorted((left.pdb.serial, right.pdb.serial)))
            bonds[pair] = order
    return [(*pair, order) for pair, order in sorted(bonds.items())]
