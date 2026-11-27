"""Build local structures for the two states of each acidic or amide side chain."""

from math import cos, radians, sin, sqrt
from pathlib import Path

import gemmi
import numpy as np

from .i18n import tr

from .analysis import angle, unit
from .structure import CARBOXYLATE_OXYGENS, FUNCTIONAL_GROUPS, SiteAtom, SiteResidue, Structure, write_local_pdb

SIDECHAIN_STATES = {
    "ASP": ("deprotonated", "protonated"), "GLU": ("deprotonated", "protonated"),
    "ASN": ("no_flip", "flip"), "GLN": ("no_flip", "flip"),
}
ACID_NAMES = {"ASP": "ASH", "GLU": "GLH"}
AMIDE_H_NAMES = {"ASN": ("HD21", "HD22"), "GLN": ("HE21", "HE22")}
ACID_O_H = 0.96
ACID_C_O_H = 108.5
ACID_LENGTH_DIFFERENCE = 0.04


def functional_group(target: SiteResidue) -> tuple[SiteAtom, SiteAtom, SiteAtom]:
    names = FUNCTIONAL_GROUPS[target.kind]
    if missing := [name for name in names if name not in target.atoms]:
        raise ValueError(tr("error.missing_sidechain_atoms", atoms=", ".join(missing)))
    carbon, oxygen, other = (target.atoms[name] for name in names)
    expected = "O" if target.kind in ACID_NAMES else "N"
    first, second = oxygen.position - carbon.position, other.position - carbon.position
    if (carbon.element, oxygen.element, other.element) != ("C", "O", expected):
        raise ValueError(tr("error.bad_sidechain_element"))
    if (not all(1.08 <= np.linalg.norm(vector) <= 1.65 for vector in (first, second))
            or not 90 <= angle(first, second) <= 150):
        raise ValueError(tr("error.bad_sidechain_geometry"))
    return carbon, oxygen, other


def polar_hydrogens(target: SiteResidue) -> dict[SiteAtom, SiteAtom]:
    """Map side-chain H/D to their parent using names and the nearest covalent heavy atom."""
    names = FUNCTIONAL_GROUPS[target.kind][1:]
    heavy = [atom for atom in target.atoms.values() if not atom.pdb.is_hydrogen()]
    if not heavy:
        return {}
    if target.kind in ACID_NAMES:
        named_parents = {f"H{name[1:]}": name for name in names}
        parent_names = set(names)
    else:
        named_parents = dict.fromkeys(AMIDE_H_NAMES[target.kind], names[1])
        parent_names = {names[1]}
    hydrogens = {}
    for hydrogen in target.atoms.values():
        if not hydrogen.pdb.is_hydrogen():
            continue
        known = named_parents.get("H" + hydrogen.name[1:])
        if known in target.atoms and hydrogen.name[0] in "HD":
            hydrogens[hydrogen] = target.atoms[known]
            continue
        parent = min(heavy, key=lambda atom: np.linalg.norm(atom.position - hydrogen.position))
        if parent.name in parent_names and 0.60 <= np.linalg.norm(parent.position - hydrogen.position) <= 1.30:
            hydrogens[hydrogen] = parent
    return hydrogens


def protonation_site(source: Structure, target: SiteResidue) -> dict:
    carbon, first, second = functional_group(target)
    oxygens = first, second
    lengths = {atom.name: float(np.linalg.norm(atom.position - carbon.position)) for atom in oxygens}
    orders = {atom.name: source.bond_orders.get(tuple(sorted((carbon.pdb.serial, atom.pdb.serial))), 0)
              for atom in oxygens}
    heavy_carbonyl = [atom for atom in oxygens if orders[atom.name] == 2]
    bound = set(polar_hydrogens(target).values())
    if len(heavy_carbonyl) == 1:
        oxygen = second if heavy_carbonyl[0] is first else first
        basis = "conect_double_bond"
    elif abs(lengths[first.name] - lengths[second.name]) >= ACID_LENGTH_DIFFERENCE:
        oxygen = max(oxygens, key=lambda atom: lengths[atom.name])
        basis = "longer_C_O_bond"
    elif len(bound) == 1:
        oxygen = next(iter(bound))
        basis = "input_hydrogen"
    else:
        oxygen = second
        basis = "amber_oxygen_2"
    return {"input_atom": oxygen.name, "output_atom": second.name, "basis": basis,
            "C_O_lengths_A": {name: round(length, 3) for name, length in lengths.items()}}


def acid_hydrogen(carbon: np.ndarray, oxygen: np.ndarray, carbonyl: np.ndarray) -> np.ndarray:
    axis = unit(carbon - oxygen)
    toward_carbonyl = carbonyl - oxygen
    in_plane = unit(toward_carbonyl - np.dot(toward_carbonyl, axis) * axis)
    theta = radians(ACID_C_O_H)
    return oxygen + ACID_O_H * (cos(theta) * axis + sin(theta) * in_plane)


def amide_hydrogens(carbon: np.ndarray, oxygen: np.ndarray, nitrogen: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    axis = unit(carbon - nitrogen)
    toward_oxygen = oxygen - nitrogen
    in_plane = unit(toward_oxygen - np.dot(toward_oxygen, axis) * axis)
    # Unit N-H vectors at 120 degrees to each other and to the N-C vector.
    return nitrogen - 0.5 * axis + sqrt(3) / 2 * in_plane, nitrogen - 0.5 * axis - sqrt(3) / 2 * in_plane


def attached_amide_positions(target: SiteResidue, nitrogen: SiteAtom) -> tuple[np.ndarray, ...]:
    return tuple(atom.position for atom, parent in polar_hydrogens(target).items() if parent is nitrogen
                 and 0.75 <= np.linalg.norm(atom.position - parent.position) <= 1.25)


def nonstandard_group_links(source: Structure, group: tuple[SiteAtom, ...]) -> list[str]:
    return [f"{atom.name} to {partner.label} ({kind})" for atom in group
            for partner, kind in source.declared_partners.get(atom, {}).items()
            if partner.residue is not atom.residue and not partner.pdb.is_hydrogen()
            and (atom.residue.kind not in ACID_NAMES or not partner.pdb.element.is_metal)]


def write_sidechain_pdb(source: Structure, neighbors: list[SiteResidue], target: SiteResidue, path: Path,
                        state: str, protonation_oxygen: str | None = None) -> None:
    default, alternative = SIDECHAIN_STATES[target.kind]
    if state not in (default, alternative):
        raise ValueError(tr("error.invalid_state", residue=target.kind, state=state))
    is_acid = target.kind in ACID_NAMES
    target_name = target.kind
    omitted = set(polar_hydrogens(target)) if is_acid else set()
    hydrogens = {}
    edits = {}
    try:
        carbon, oxygen, other = functional_group(target)
    except ValueError:
        if state == alternative:
            raise
    else:
        group = (carbon, oxygen, other)
        linked = nonstandard_group_links(source, group)
        if state == alternative and linked:
            raise ValueError(tr("error.nonstandard_bond"))
        if state == "protonated":
            input_name = protonation_oxygen or protonation_site(source, target)["input_atom"]
            if input_name not in (oxygen.name, other.name):
                raise ValueError(tr("error.bad_carboxyl_oxygen", atom=input_name))
            hydroxyl = target.atoms[input_name]
            carbonyl = other if hydroxyl is oxygen else oxygen
            first_name, second_name = CARBOXYLATE_OXYGENS[target.kind]
            target_name = ACID_NAMES[target.kind]
            for atom, name in ((carbonyl, first_name), (hydroxyl, second_name)):
                replacement = atom.pdb.clone()
                replacement.name = name
                replacement.charge = 0
                edits[atom] = replacement
            hydrogens["H" + second_name[1:]] = (
                hydroxyl, acid_hydrogen(carbon.position, hydroxyl.position, carbonyl.position))
        elif not is_acid and not linked:
            flipped = state == "flip"
            polar_o, polar_n = (other.position, oxygen.position) if flipped else (oxygen.position, other.position)
            if flipped:
                for atom, point in ((oxygen, polar_o), (other, polar_n)):
                    replacement = atom.pdb.clone()
                    replacement.pos = gemmi.Position(*point)
                    edits[atom] = replacement
            existing = polar_hydrogens(target)
            if flipped or len(existing) != 2 or len(attached_amide_positions(target, other)) != 2:
                omitted = set(existing)
                positions = amide_hydrogens(carbon.position, polar_o, polar_n)
                hydrogens = {name: (other, point) for name, point in zip(AMIDE_H_NAMES[target.kind], positions)}
    write_local_pdb(source, neighbors, target, path, target_name=target_name, omitted=omitted,
                    hydrogens=hydrogens, atom_edits=edits)
