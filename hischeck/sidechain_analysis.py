"""Estimate solvent exposure and directional contact support for acidic and amide side chains."""

from math import pi
from pathlib import Path

import numpy as np

from .analysis import HEAVY_MIN, STRONG_METAL_LIMIT, Settings, angle, backbone_hydrogen, describe, hydrogen_bond, unit
from .chemistry import Role, WEAK_ACCEPTOR_ELEMENTS, is_water, role_for
from .sidechain import (SIDECHAIN_STATES, acid_hydrogen, amide_hydrogens, attached_amide_positions, functional_group,
                        nonstandard_group_links, polar_hydrogens, protonation_site, write_sidechain_pdb)
from .structure import FUNCTIONAL_GROUPS, SiteAtom, SiteResidue, Structure, attached_hydrogens

PROBE_RADIUS = 1.4
PROBE_REACH = 8.0
OPEN_FRACTION = 0.15
SAMPLES = 96
CATION_LIMIT = 4.0
CATION_CONTACTS = {"LYS": ({"NZ"}, "lysine"), "ARG": ({"NE", "NH1", "NH2"}, "arginine")}
MIN_WATER_O_DISTANCE = 2.20
MIN_AMIDE_SUPPORT = 0.24
MIN_ACID_SUPPORT = 0.20
METAL_LIMITS = {**STRONG_METAL_LIMIT, "CA": 2.85, "MG": 2.60, "NA": 2.75, "K": 3.05, "CD": 2.75}


def sphere_directions() -> np.ndarray:
    index = np.arange(SAMPLES)
    z = 1 - (index + 0.5) * 2 / SAMPLES
    theta = pi * (3 - np.sqrt(5)) * index
    radius = np.sqrt(1 - z**2)
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))


DIRECTIONS = sphere_directions()


def solvent_exposure(source: Structure, group: tuple[SiteAtom, SiteAtom, SiteAtom]) -> dict:
    carbon, first, second = group
    axis = unit(first.position - carbon.position)
    normal = unit(np.cross(axis, second.position - carbon.position))
    # The sample rotates with the side chain, independently of the input coordinate frame.
    rays = DIRECTIONS @ np.stack((axis, np.cross(normal, axis), normal))
    fractions = {}
    for center in (first, second):
        # Deposited water is contact evidence, so it never blocks the water probe.
        blockers = [atom for atom, _ in source.nearest(center.position, PROBE_REACH + 5.0)
                    if atom is not center and not is_water(atom.residue) and atom.pdb.occ > 0]
        offsets = np.asarray([atom.position - center.position for atom in blockers]).reshape(-1, 3)
        radii = np.array([atom.pdb.element.vdw_r + PROBE_RADIUS for atom in blockers])
        projection = offsets @ rays.T
        along = np.clip(projection, center.pdb.element.vdw_r + PROBE_RADIUS, PROBE_REACH)
        square_distances = np.sum(offsets**2, axis=1)[:, None] - 2 * projection * along + along**2
        open_rays = ~np.any(square_distances < radii[:, None]**2, axis=0)
        fractions[center.name] = float(np.mean(open_rays))
    return {"exposed": max(fractions.values()) >= OPEN_FRACTION, "probe_radius_A": PROBE_RADIUS,
            "reach_A": PROBE_REACH, "directions_per_atom": SAMPLES, "threshold": OPEN_FRACTION,
            "open_fraction": {name: round(value, 4) for name, value in fractions.items()}}


def contact(atom: SiteAtom, center: str, distance: float, interaction: str, **details) -> dict:
    return {**describe(atom), "center": center, "distance_A": round(distance, 3), "interaction": interaction, **details}


def possible_donor(source: Structure, donor: SiteAtom, acceptor: np.ndarray) -> float:
    hydrogens = [atom.position for atom in attached_hydrogens(donor)]
    if not hydrogens and (estimated := backbone_hydrogen(source, donor)) is not None:
        hydrogens = [estimated]
    if (not hydrogens and donor.residue.kind in {"ASN", "GLN"}
            and donor.name == FUNCTIONAL_GROUPS[donor.residue.kind][-1]):
        try:
            carbon, oxygen, _ = functional_group(donor.residue)
        except ValueError:
            pass
        else:
            hydrogens = amide_hydrogens(carbon.position, oxygen.position, donor.position)
    if hydrogens:
        return max(hydrogen_bond(donor.position, point, acceptor)[0] for point in hydrogens)
    return 0.8  # H coordinates are absent, so this is a possible heavy-atom contact.


def acid_partners(source: Structure, group: tuple[SiteAtom, SiteAtom, SiteAtom], neighbors: dict,
                  settings: Settings, overrides: dict[tuple[str, str], Role]) -> list[dict]:
    records = []
    for center in group[1:]:
        for atom, distance in neighbors[center]:
            role = role_for(atom, overrides)
            cation_names, cation_kind = CATION_CONTACTS.get(atom.residue.name, ((), ""))
            if atom.pdb.element.is_metal:
                declared = source.declared_partners.get(center, {}).get(atom)
                limit = 3.05 if declared else METAL_LIMITS.get(atom.element, 0)
                if 1.65 <= distance <= limit and atom.pdb.occ >= 0.50:
                    records.append(contact(atom, center.name, distance, "metal", declaration=declared))
            elif atom.name in cation_names and role.donor and HEAVY_MIN <= distance <= CATION_LIMIT:
                records.append(contact(atom, center.name, distance, cation_kind))
            elif atom.element == "O" and is_water(atom.residue):
                if MIN_WATER_O_DISTANCE <= distance <= settings.hbond_max:
                    records.append(contact(atom, center.name, distance, "water"))
            elif atom.element in {"N", "O"} and role.donor >= 0.40 and HEAVY_MIN <= distance <= settings.hbond_max:
                strength = role.donor * min(1.0, atom.pdb.occ) * possible_donor(source, atom, center.position)
                if strength >= MIN_ACID_SUPPORT:
                    records.append(contact(atom, center.name, distance, "hydrogen", strength=round(strength, 3)))
    return sorted(records, key=lambda item: item["distance_A"])


def amide_partners(source: Structure, group: tuple[SiteAtom, SiteAtom, SiteAtom], neighbors: dict, flipped: bool,
                   settings: Settings, overrides: dict[tuple[str, str], Role]) -> list[dict]:
    carbon, oxygen, nitrogen = group
    acceptor, donor = (nitrogen, oxygen) if flipped else (oxygen, nitrogen)
    normal = unit(np.cross(acceptor.position - carbon.position, donor.position - carbon.position))
    points = () if flipped else attached_amide_positions(nitrogen.residue, nitrogen)
    if len(points) != 2:
        points = amide_hydrogens(carbon.position, acceptor.position, donor.position)
    records = []
    for center, name in ((acceptor, oxygen.name), (donor, nitrogen.name)):
        for atom, distance in neighbors[center]:
            role = role_for(atom, overrides)
            is_sulfur = atom.element in WEAK_ACCEPTOR_ELEMENTS
            limit = settings.s_bond_max if is_sulfur else settings.hbond_max
            weight = role.donor if center is acceptor else role.acceptor
            if not weight or not HEAVY_MIN <= distance <= limit:
                continue
            direction = unit(atom.position - center.position)
            plane_angle = abs(90 - angle(normal, direction))
            plane_weight = float(np.clip((55 - plane_angle) / 35, 0, 1))
            water = is_water(atom.residue)
            if center is acceptor:
                out_angle = angle(center.position - carbon.position, direction)
                geometry = possible_donor(source, atom, center.position) * float(np.clip((90 - out_angle) / 40, 0, 1))
                weight = max(weight, 0.55) if water else weight
                flip_support = not water and weight >= 0.60
            else:
                geometry = max(hydrogen_bond(center.position, point, atom.position)[0] for point in points)
                flip_support = not role.donor
            strength = geometry * plane_weight * weight * min(1.0, atom.pdb.occ)
            records.append(contact(atom, name, distance, "hydrogen", role=role.label, plane_angle_deg=round(plane_angle, 1),
                                   strength=round(strength, 3), supports=strength >= MIN_AMIDE_SUPPORT,
                                   supports_flip=flip_support and strength >= MIN_AMIDE_SUPPORT))
    return sorted(records, key=lambda item: (-item["strength"], item["distance_A"]))


def trial_clashes(source: Structure, group: tuple[SiteAtom, SiteAtom, SiteAtom], points) -> list[dict]:
    excluded = set(group) | polar_hydrogens(group[0].residue).keys()
    limits = {"C": 1.65, "N": 1.37, "O": 1.37, "H": 1.18, "D": 1.18}
    return [{**describe(atom), "distance_A": round(distance, 3), "severe": True}
            for point in points for atom, distance in source.nearest(point, 1.65, heavy_only=False)
            if atom not in excluded and atom.pdb.occ > 0 and distance < limits.get(atom.element, 1.57)]


def sidechain_evidence(source: Structure, target: SiteResidue, settings: Settings,
                       overrides: dict[tuple[str, str], Role]) -> dict:
    default, alternative = SIDECHAIN_STATES[target.kind]
    data = {**target.site_metadata(), "status": "default", "recommendation": default, "top_candidate": default,
            "reasons": [], "candidates": [], "contacts": {}, "metal_contacts": [], "carboxylate_contacts": [],
            "available_states": [default]}
    if missing := [name for name in FUNCTIONAL_GROUPS[target.kind] if name not in target.atoms]:
        data.update(reasons=["missing_functional_group_atoms"], missing_functional_group_atoms=sorted(missing))
        return data
    try:
        group = functional_group(target)
    except ValueError:
        data["reasons"] = ["invalid_functional_group_geometry"]
        return data
    if links := nonstandard_group_links(source, group):
        data.update(reasons=["nonstandard_sidechain_link"], nonstandard_links=links)
        return data

    carbon, oxygen, other = group
    exposure = solvent_exposure(source, group)
    data.update(surface_exposure=exposure, available_states=[default, alternative])
    neighbors = {center: [(atom, distance) for atom, distance in source.nearest(center.position, settings.search_radius)
                          if atom not in group and atom.pdb.occ > 0] for center in group[1:]}
    if default == "deprotonated":
        data["protonation_oxygen"] = protonation_site(source, target)
        hydroxyl = target.atoms[data["protonation_oxygen"]["input_atom"]]
        carbonyl = other if hydroxyl is oxygen else oxygen
        points = (acid_hydrogen(carbon.position, hydroxyl.position, carbonyl.position),)
        default_supports = acid_partners(source, group, neighbors, settings, overrides)
        alternative_supports = []
        acid_reasons = {"metal": "sidechain_metal_coordination", "water": "nearby_crystal_water",
                        "lysine": "nearby_positive_lysine", "arginine": "nearby_positive_arginine",
                        "hydrogen": "nearby_hbond_donor"}
        why = list(dict.fromkeys(acid_reasons[item["interaction"]] for item in default_supports))
        changed = not why
        change_reason = "buried_without_donors"
    else:
        points = amide_hydrogens(carbon.position, other.position, oxygen.position)
        no_flip = amide_partners(source, group, neighbors, False, settings, overrides)
        flip = amide_partners(source, group, neighbors, True, settings, overrides)
        data["contacts"] = {default: no_flip, alternative: flip}
        default_supports = [item for item in no_flip if item["supports"]]
        alternative_supports = [item for item in flip if item["supports"]]
        changed = not default_supports and any(item["supports_flip"] for item in alternative_supports)
        why = ["original_amide_hbond" if default_supports else "no_directional_flip_support"]
        change_reason = "flipped_amide_hbond"

    clashes = trial_clashes(source, group, points)
    data["candidates"] = [{"state": default, "supports": default_supports, "clashes": []},
                          {"state": alternative, "supports": alternative_supports, "clashes": clashes}]
    if exposure["exposed"]:
        why = ["sidechain_exposed", *why] if default_supports else ["sidechain_exposed"]
    elif changed:
        if clashes:
            why = ["sidechain_hydrogen_clash"]
        elif min(atom.pdb.occ for atom in group) < 0.75 or target.options:
            why = ["uncertain_sidechain_conformation"]
        else:
            data.update(status="recommended", recommendation=alternative, top_candidate=alternative)
            why = [change_reason]
    data["reasons"] = why
    return data


def export_sidechain(source: Structure, neighbors: list[SiteResidue], target: SiteResidue, folder: Path,
                     settings: Settings, overrides: dict[tuple[str, str], Role]) -> dict:
    data = sidechain_evidence(source, target, settings, overrides)
    default, alternative = SIDECHAIN_STATES[target.kind]
    data["exports"] = [{"pdb": "site.pdb", "state": default}]
    data["automatic_states"] = [data["recommendation"]]
    if data["recommendation"] == alternative:
        data["exports"].append({"pdb": f"site_{alternative}.pdb", "state": alternative})
    proton = data.get("protonation_oxygen", {}).get("input_atom")
    for item in data["exports"]:
        write_sidechain_pdb(source, neighbors, target, folder / item["pdb"], item["state"], proton)
    if default == "deprotonated" and (removed := polar_hydrogens(target)):
        data["removed_sidechain_hydrogens"] = sorted(atom.name for atom in removed)
    return data
