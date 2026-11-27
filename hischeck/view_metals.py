"""Find coordination links for display in the local structure viewer."""

from collections.abc import Iterator

from .analysis import STRONG_METAL_LIMIT
from .chemistry import HIS_ALIAS, HIS_PROTONS, PROTEIN
from .structure import SiteAtom, Structure, attached_hydrogens

METAL_LIMITS = {**STRONG_METAL_LIMIT, "MG": 2.45, "CA": 2.85, "NA": 2.90, "K": 3.35,
                "CD": 2.85, "HG": 2.85}
DONOR_ELEMENTS = {"N", "O", "S", "SE"}


def can_coordinate(atom: SiteAtom) -> bool:
    if atom.element not in DONOR_ELEMENTS:
        return False
    if atom.element != "N":
        return True
    residue = atom.residue
    if residue.name in PROTEIN:
        if atom.name == "N" or (residue.kind, atom.name) in {("ASN", "ND2"), ("GLN", "NE2")}:
            return False
        if residue.kind == "HIS":
            state = HIS_ALIAS.get(residue.name, residue.name)
            return atom.name in {"ND1", "NE2"} and atom.name not in HIS_PROTONS.get(state, ()) and not attached_hydrogens(atom)
        return (residue.name, atom.name) == ("LYN", "NZ")
    # A ligand amino N may carry one or two H and still coordinate the metal.
    return len(attached_hydrogens(atom)) < 3


def coordination_pairs(structure: Structure, metals: list[SiteAtom]) -> Iterator[tuple[SiteAtom, SiteAtom, float]]:
    for metal in metals:
        declared = structure.declared_partners.get(metal, {})
        cutoff = METAL_LIMITS.get(metal.element, 2.7)
        for atom, distance in structure.nearest(metal.position, max(cutoff, 3.35)):
            if atom is metal or atom.pdb.element.is_metal or not 1.6 <= distance:
                continue
            if atom in declared or (distance <= cutoff and can_coordinate(atom)):
                yield metal, atom, distance
