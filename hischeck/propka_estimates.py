"""Estimate side-chain pKa once for the selected full PDB model."""

from importlib.metadata import version
from io import StringIO
from math import isfinite

import gemmi

from .structure import ResidueID, Structure

TITRATABLE = frozenset({"HIS", "ASP", "GLU"})


def selected_pdb(source: Structure) -> str:
    # Keep the original model's chain/TER topology, heteroatoms and coordinates.
    structure = gemmi.Structure()
    model = structure.add_model(next(model for model in source.pdb if model.num == source.model_number).clone())
    for chain in model:
        for residue in chain:
            key = ResidueID(chain.name.strip(), residue.seqid.num, residue.seqid.icode.strip(" \0"))
            selected = source.residues_by_name[key, residue.name]
            for index in reversed(range(len(residue))):
                atom = residue[index]
                chosen = selected.atoms.get(atom.name)
                if chosen is None or atom.altloc != chosen.pdb.altloc or atom.is_hydrogen():
                    del residue[index]
                else:
                    atom.altloc = "\0"
            if selected.kind in TITRATABLE and selected.kind != residue.name:
                residue.name = selected.kind
    return structure.make_minimal_pdb()


def estimate_pkas(source: Structure) -> tuple[dict[ResidueID, float], dict]:
    """Return estimates and provenance; unsupported residues have no entry."""
    try:
        from propka.run import single
    except ImportError:
        return {}, {"status": "unavailable"}

    info = {"method": "PROPKA", "version": version("propka")}
    try:
        molecule = single("selected_model.pdb", stream=StringIO(selected_pdb(source)), write_pka=False)
        groups = molecule.conformations[molecule.conformation_names[0]].get_titratable_groups()
        values = {}
        for group in groups:
            if group.residue_type not in TITRATABLE or not isfinite(group.pka_value):
                continue
            atom = group.atom
            key = ResidueID("" if atom.chain_id == "_" else atom.chain_id, atom.res_num, atom.icode.strip())
            residue = source.residues.get(key)
            if residue is not None and residue.kind == group.residue_type:
                values[key] = round(group.pka_value, 4)
    except (Exception, SystemExit) as error:
        return {}, {**info, "status": "failed", "detail": str(error)}
    return values, {**info, "status": "completed"}
