"""PDB input, consistent conformer selection and local structure export."""

from collections import Counter, defaultdict
from collections.abc import Collection
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import gemmi
import numpy as np

from .i18n import tr

HIS_NAMES = frozenset({"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"})
RING = ("CG", "ND1", "CE1", "NE2", "CD2")
NITROGENS = ("ND1", "NE2")
CARBOXYLATE_OXYGENS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
FUNCTIONAL_GROUPS = {
    "HIS": RING,
    "ASP": ("CG", *CARBOXYLATE_OXYGENS["ASP"]),
    "GLU": ("CD", *CARBOXYLATE_OXYGENS["GLU"]),
    "ASN": ("CG", "OD1", "ND2"),
    "GLN": ("CD", "OE1", "NE2"),
}
RESIDUE_TYPES = tuple(FUNCTIONAL_GROUPS)
RESIDUE_ALIASES = {**dict.fromkeys(HIS_NAMES, "HIS"), "ASH": "ASP", "GLH": "GLU"}


def residue_type(name: str) -> str:
    return RESIDUE_ALIASES.get(name, name)


def clean_char(value: str) -> str:
    return value.strip(" \0")


def path_token(text: str) -> str:
    return "".join(c if c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789" else f"~{ord(c):02X}" for c in text)


@dataclass(frozen=True)
class ResidueID:
    chain: str
    number: int
    insertion: str = ""

    @property
    def label(self) -> str:
        return f"{self.chain or '_'}:{self.number}{self.insertion}"

    def directory_name(self, kind: str) -> str:
        return f"{path_token(self.chain) or 'blank'}_{self.number}{path_token(self.insertion)}_{kind}"

    @property
    def directory(self) -> str:
        return self.directory_name("HIS")


@dataclass(eq=False)
class SiteAtom:
    residue: "SiteResidue"
    pdb: gemmi.Atom
    position: np.ndarray

    @property
    def name(self) -> str:
        return self.pdb.name

    @property
    def element(self) -> str:
        return self.pdb.element.name.upper()

    @property
    def label(self) -> str:
        return f"{self.residue.id.label} {self.residue.name} {self.name}"


@dataclass(eq=False)
class SiteResidue:
    id: ResidueID
    pdb: gemmi.Residue
    atoms: dict[str, SiteAtom]
    altloc: str
    options: tuple[str, ...]

    @property
    def name(self) -> str:
        return self.pdb.name

    @property
    def kind(self) -> str:
        return residue_type(self.name)

    @property
    def functional_atoms(self) -> list[SiteAtom]:
        return [self.atoms[name] for name in FUNCTIONAL_GROUPS.get(self.kind, ()) if name in self.atoms]

    def site_metadata(self) -> dict:
        return {"site": self.id.label, "chain": self.id.chain, "residue_number": self.id.number,
                "insertion_code": self.id.insertion, "input_residue": self.name,
                "target_altloc": self.altloc, "target_altloc_options": list(self.options)}


def select_residue(chain: str, residue: gemmi.Residue) -> SiteResidue:
    occupancies = defaultdict(list)
    for atom in residue:
        if altloc := clean_char(atom.altloc):
            occupancies[altloc].append(atom.occ)
    options = tuple(sorted(occupancies))
    preferred = min(options, key=lambda label: (-float(np.mean(occupancies[label])), label)) if options else ""
    key = ResidueID(chain, residue.seqid.num, clean_char(residue.seqid.icode))
    chosen = SiteResidue(key, residue, {}, preferred, options)
    for atom in residue:
        altloc = clean_char(atom.altloc)
        if altloc and altloc != preferred:
            continue
        existing = chosen.atoms.get(atom.name)
        if existing is not None:
            if existing.pdb.altloc == atom.altloc:
                raise ValueError(tr("error.duplicate_atom", site=key.label, residue=residue.name, atom=atom.name))
            if not altloc:
                continue
        position = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
        if not np.all(np.isfinite(position)):
            raise ValueError(tr("error.bad_coordinates", site=key.label, residue=residue.name, atom=atom.name))
        chosen.atoms[atom.name] = SiteAtom(chosen, atom, position)
    return chosen


class Structure:
    def __init__(self, path: Path, model_number: int | None = None):
        pdb_bytes = path.read_bytes()
        self.sha256 = sha256(pdb_bytes).hexdigest()
        self.pdb = gemmi.read_pdb_string(pdb_bytes)
        models = [model for model in self.pdb if model_number is None or model.num == model_number]
        if not models:
            raise ValueError(tr("error.no_model", path=path))
        model = models[0]
        self.model_number = model.num
        self.model_count = len(self.pdb)
        self.residues = {}
        for chain in model:
            for residue in chain:
                selected = select_residue(chain.name.strip(), residue)
                if selected.id in self.residues:
                    raise ValueError(tr("error.duplicate_residue", site=selected.id.label))
                self.residues[selected.id] = selected
        if not self.residues:
            raise ValueError(tr("error.empty_model", model=self.model_number))
        self.atoms = [atom for residue in self.residues.values() for atom in residue.atoms.values()]
        self.positions = np.asarray([atom.position for atom in self.atoms])
        self.heavy_atoms = [atom for atom in self.atoms if not atom.pdb.is_hydrogen()]
        self.heavy_positions = np.asarray([atom.position for atom in self.heavy_atoms]).reshape(-1, 3)
        self.bond_orders = {}
        for serial, neighbors in self.pdb.conect_map.items():
            for other, count in Counter(neighbors).items():
                if serial != other:
                    pair = tuple(sorted((serial, other)))
                    self.bond_orders[pair] = max(count, self.bond_orders.get(pair, 0))
        if self.bond_orders and len({atom.pdb.serial for atom in self.atoms}) != len(self.atoms):
            raise ValueError(tr("error.duplicate_serials"))
        self.links = []
        self.declared_partners: dict[SiteAtom, dict[SiteAtom, str]] = defaultdict(dict)
        serial_atoms = {atom.pdb.serial: atom for atom in self.atoms}
        for first_serial, second_serial in self.bond_orders:
            first, second = serial_atoms.get(first_serial), serial_atoms.get(second_serial)
            if first is not None and second is not None:
                self.declared_partners[first][second] = "CONECT"
                self.declared_partners[second][first] = "CONECT"
        for link in self.pdb.connections:
            if link.asu == gemmi.Asu.Different:
                continue
            first, second = self.resolve(link.partner1), self.resolve(link.partner2)
            if first is not None and second is not None:
                self.links.append((link, first, second))
                self.declared_partners[first][second] = "LINK"
                self.declared_partners[second][first] = "LINK"

    def resolve(self, address: gemmi.AtomAddress) -> SiteAtom | None:
        key = ResidueID(address.chain_name.strip(), address.res_id.seqid.num, clean_char(address.res_id.seqid.icode))
        residue = self.residues.get(key)
        if residue is None or (address.res_id.name and residue.name != address.res_id.name):
            return None
        atom = residue.atoms.get(address.atom_name)
        altloc = clean_char(address.altloc)
        return atom if atom is not None and (not altloc or clean_char(atom.pdb.altloc) == altloc) else None

    def nearest(self, center: np.ndarray, radius: float, *, heavy_only=True) -> list[tuple[SiteAtom, float]]:
        atoms, positions = (self.heavy_atoms, self.heavy_positions) if heavy_only else (self.atoms, self.positions)
        offsets = positions - center
        squared_distances = np.einsum("ij,ij->i", offsets, offsets)
        return [(atoms[i], float(np.sqrt(squared_distances[i]))) for i in np.flatnonzero(squared_distances <= radius**2)]

    def neighborhood(self, target: SiteResidue, radius: float) -> list[SiteResidue]:
        centers = target.functional_atoms
        if not centers:
            centers = [target.atoms["CA"]] if "CA" in target.atoms else list(target.atoms.values())
        keys = {target.id}
        for center in centers:
            keys.update(atom.residue.id for atom, _ in self.nearest(center.position, radius))
        return [residue for key, residue in self.residues.items() if key in keys]


def attached_hydrogens(atom: SiteAtom) -> list[SiteAtom]:
    return [other for other in atom.residue.atoms.values() if other.pdb.is_hydrogen()
            and 0.60 <= np.linalg.norm(other.position - atom.position) <= 1.30]


def ring_hydrogens(residue: SiteResidue) -> set[SiteAtom]:
    heavy = [atom for atom in residue.atoms.values() if not atom.pdb.is_hydrogen()]
    ring_names = {"HD1", "HE2", "HE1", "HD2", "DD1", "DE2", "DE1", "DD2"}
    found = set()
    for hydrogen in residue.atoms.values():
        if not hydrogen.pdb.is_hydrogen():
            continue
        parent = min(heavy, key=lambda atom: np.linalg.norm(atom.position - hydrogen.position), default=None)
        if hydrogen.name in ring_names or (parent is not None and parent.name in RING
                                           and np.linalg.norm(parent.position - hydrogen.position) <= 1.30):
            found.add(hydrogen)
    return found


def export_pdb(source: Structure, residues: list[SiteResidue], target: SiteResidue, path: Path,
               state: str | None, hydrogens: dict[str, tuple[SiteAtom, np.ndarray]]) -> None:
    write_local_pdb(source, residues, target, path, target_name=state,
                    omitted=ring_hydrogens(target) if state else (), hydrogens=hydrogens if state else None)


def write_local_pdb(source: Structure, residues: list[SiteResidue], target: SiteResidue, path: Path, *,
                    target_name: str | None = None, omitted: Collection[SiteAtom] = (),
                    hydrogens: dict[str, tuple[SiteAtom, np.ndarray]] | None = None,
                    atom_edits: dict[SiteAtom, gemmi.Atom] | None = None) -> None:
    local = gemmi.Structure()
    atom_edits = atom_edits or {}
    model = gemmi.Model(1)
    emitted = {}
    chain = None
    serial = 0
    for residue in residues:
        if chain is None or chain.name != residue.id.chain:
            chain = model.add_chain(residue.id.chain)
        new_residue = gemmi.Residue()
        new_residue.name = target_name if target_name and residue is target else residue.name
        new_residue.seqid = residue.pdb.seqid
        new_residue.het_flag = residue.pdb.het_flag
        new_residue.entity_type = residue.pdb.entity_type
        for atom in residue.atoms.values():
            if atom in omitted:
                continue
            copy = atom_edits.get(atom, atom.pdb).clone()
            copy.altloc = "\0"
            serial += 1
            copy.serial = serial
            new_residue.add_atom(copy)
            emitted[atom] = serial
        if residue is target and hydrogens:
            for name, (parent, point) in hydrogens.items():
                hydrogen = gemmi.Atom()
                hydrogen.name = name
                hydrogen.element = gemmi.Element("H")
                hydrogen.pos = gemmi.Position(*point)
                hydrogen.occ = parent.pdb.occ
                hydrogen.b_iso = parent.pdb.b_iso
                serial += 1
                hydrogen.serial = serial
                new_residue.add_atom(hydrogen)
        chain.add_residue(new_residue)
    local.add_model(model)
    serial_map = {atom.pdb.serial: number for atom, number in emitted.items()}
    for (first, second), order in source.bond_orders.items():
        if first in serial_map and second in serial_map:
            local.add_conect(serial_map[first], serial_map[second], order)

    def address(atom: SiteAtom) -> gemmi.AtomAddress:
        residue = atom.residue
        name = target_name if target_name and residue is target else residue.name
        return gemmi.AtomAddress(residue.id.chain, residue.pdb.seqid, name, atom_edits.get(atom, atom.pdb).name)

    for original, first, second in source.links:
        if first in emitted and second in emitted:
            link = gemmi.Connection()
            link.name, link.type, link.asu = original.name, original.type, gemmi.Asu.Same
            link.partner1, link.partner2 = address(first), address(second)
            link.reported_distance = original.reported_distance
            local.connections.append(link)

    options = gemmi.PdbWriteOptions()
    options.cryst1_record = False
    options.seqres_records = False
    options.preserve_serial = True
    options.numbered_ter = False
    options.conect_records = True
    path.write_text(local.make_pdb_string(options), encoding="ascii")
