"""Donor and acceptor roles for common PDB protein and solvent atoms."""

import csv
from dataclasses import dataclass
from pathlib import Path

from .i18n import tr

from .structure import CARBOXYLATE_OXYGENS, HIS_NAMES, NITROGENS, SiteAtom, SiteResidue, attached_hydrogens


@dataclass(frozen=True)
class Role:
    label: str
    donor: float = 0.0
    acceptor: float = 0.0


DONOR = Role("donor", donor=1.0)
ACCEPTOR = Role("acceptor", acceptor=1.0)
BOTH = Role("donor/acceptor", donor=0.55, acceptor=0.55)
AMINE = Role("amine", donor=0.68, acceptor=0.30)
NONE = Role("neither")
UNKNOWN = Role("unknown")
WATER = Role("water", donor=0.32, acceptor=0.32)
HYDROXYL = Role("hydroxyl", donor=0.58, acceptor=0.44)
SULFUR = Role("weak sulfur acceptor", acceptor=0.28)
HIS_UNKNOWN = Role("unassigned His", donor=0.36, acceptor=0.36)
PROTEIN = set("ALA ARG ASN ASP ASH CYS CYM CYX GLN GLU GLH GLY ILE LEU LYS LYN MET MSE PHE PRO SER THR TRP TYR VAL".split()) | HIS_NAMES
SIDE_CHAIN = {
    "ARG": {"NE": DONOR, "NH1": DONOR, "NH2": DONOR},
    "ASN": {"OD1": ACCEPTOR, "ND2": DONOR},
    "ASP": {"OD1": ACCEPTOR, "OD2": ACCEPTOR},
    "CYS": {"SG": Role("thiol", donor=0.30, acceptor=0.23)},
    "CYM": {"SG": Role("thiolate", acceptor=0.52)},
    "CYX": {"SG": SULFUR},
    "GLN": {"OE1": ACCEPTOR, "NE2": DONOR},
    "GLU": {"OE1": ACCEPTOR, "OE2": ACCEPTOR},
    "LYS": {"NZ": DONOR},
    "LYN": {"NZ": AMINE},
    "MET": {"SD": SULFUR},
    "MSE": {"SE": SULFUR},
    "SER": {"OG": HYDROXYL},
    "THR": {"OG1": HYDROXYL},
    "TRP": {"NE1": DONOR},
    "TYR": {"OH": Role("phenol", donor=0.55, acceptor=0.30)},
}
HIS_PROTONS = {
    "HID": frozenset({"ND1"}),
    "HIE": frozenset({"NE2"}),
    "HIP": frozenset({"ND1", "NE2"}),
}
HIS_ALIAS = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
OVERRIDE_ROLES = {"donor": DONOR, "acceptor": ACCEPTOR, "both": BOTH, "amine": AMINE,
                  "neither": NONE, "unknown": UNKNOWN}
WATER_NAMES = {"WAT", "TIP", "SOL", "H2O"}
WEAK_ACCEPTOR_ELEMENTS = {"S", "SE"}


def load_overrides(path: Path | None) -> dict[tuple[str, str], Role]:
    overrides = {}
    if path is None:
        return overrides
    with path.open(encoding="utf-8-sig", newline="") as file:
        for line_number, row in enumerate(csv.DictReader(file), 2):
            try:
                residue = row["residue"].strip().upper()
                atom = row["atom"].strip().upper()
                role = OVERRIDE_ROLES[row["role"].strip().lower()]
                if not residue or not atom or (residue, atom) in overrides:
                    raise ValueError
            except (AttributeError, KeyError, ValueError) as error:
                raise ValueError(tr("error.bad_role", path=path, line=line_number)) from error
            overrides[residue, atom] = role
    if not overrides:
        raise ValueError(tr("error.no_roles", path=path))
    return overrides


def is_water(residue: SiteResidue) -> bool:
    return residue.pdb.is_water() or residue.name in WATER_NAMES


def is_primary_role(atom: SiteAtom, role: Role) -> bool:
    return role in (DONOR, ACCEPTOR) and atom.element not in WEAK_ACCEPTOR_ELEMENTS and not is_water(atom.residue)


def role_for(atom: SiteAtom, overrides: dict[tuple[str, str], Role]) -> Role:
    name, residue = atom.name, atom.residue
    if override := overrides.get((residue.name, name)):
        return override
    if is_water(residue):
        return WATER if atom.element == "O" else NONE
    if residue.name not in PROTEIN:
        return UNKNOWN if atom.element in {"N", "O", "S", "SE", "F", "CL"} else NONE
    if name in {"O", "OXT", "OT1", "OT2"}:
        return ACCEPTOR
    if name == "N":
        return DONOR if residue.name != "PRO" or attached_hydrogens(atom) else NONE
    if residue.name in HIS_NAMES and name in NITROGENS:
        state = HIS_ALIAS.get(residue.name, residue.name)
        if attached_hydrogens(atom):
            return DONOR
        other_name = "NE2" if name == "ND1" else "ND1"
        if state in HIS_PROTONS:
            return DONOR if name in HIS_PROTONS[state] else ACCEPTOR
        if (other := residue.atoms.get(other_name)) and attached_hydrogens(other):
            return ACCEPTOR
        return HIS_UNKNOWN
    if (residue.name, name) in {("ASH", "OD1"), ("ASH", "OD2"), ("GLH", "OE1"), ("GLH", "OE2")}:
        if attached_hydrogens(atom):
            return DONOR
        other_name = name[:-1] + ("2" if name[-1] == "1" else "1")
        if (other := residue.atoms.get(other_name)) and attached_hydrogens(other):
            return ACCEPTOR
        return UNKNOWN
    return SIDE_CHAIN.get(residue.name, {}).get(name, NONE)


def is_carboxylate_oxygen(atom: SiteAtom, overrides: dict[tuple[str, str], Role]) -> bool:
    residue = atom.residue
    names = CARBOXYLATE_OXYGENS.get(residue.name, ())
    if atom.name not in names:
        if residue.name not in PROTEIN or atom.name not in {"O", "OXT"} or "OXT" not in residue.atoms:
            return False
        names = ("O", "OXT")
    for name in names:
        oxygen = residue.atoms.get(name)
        if oxygen is None or oxygen.element != "O" or oxygen.pdb.charge > 0:
            return False
        if role_for(oxygen, overrides) != ACCEPTOR or attached_hydrogens(oxygen):
            return False
    return True
