"""Rank His states using directional contacts, carboxylates, coordination, and clashes."""

from dataclasses import dataclass, replace
from itertools import product
from math import acos, degrees

import numpy as np

from .chemistry import HIS_PROTONS, WEAK_ACCEPTOR_ELEMENTS, Role, is_carboxylate_oxygen, is_primary_role, role_for
from .structure import NITROGENS, RING, SiteAtom, SiteResidue, Structure, attached_hydrogens, ring_hydrogens

N_H_NAMES = {"ND1": "HD1", "NE2": "HE2"}
H_NAMES = {**N_H_NAMES, "CD2": "HD2", "CE1": "HE1"}
STRONG_METAL_LIMIT = {"ZN": 2.55, "FE": 2.55, "CU": 2.60, "NI": 2.55, "CO": 2.55, "MN": 2.65}
METAL_SEARCH = 3.05
HEAVY_MIN = 2.35
H_BOND_MAX = 3.55
S_BOND_MAX = 3.75
CARBOXYLATE_MAX = 5.0
MIN_PRIMARY_SUPPORT = 0.35
PRIMARY_CONTACT_BONUS = 0.35
MIN_SUPPORT = 0.55
MIN_MARGIN = 0.26
RING_MIN_OCC = 0.75
TIE = 1e-8


@dataclass(frozen=True)
class Settings:
    search_radius: float = 5.0
    ph: float | None = None
    pka: float = 6.0
    hbond_max: float = H_BOND_MAX
    s_bond_max: float = S_BOND_MAX
    carboxylate_max: float = CARBOXYLATE_MAX


@dataclass(frozen=True)
class Contact:
    atom: SiteAtom
    distance: float
    role: str
    axis_angle: float
    his_donor: float
    his_acceptor: float
    donor_geometry: dict | None
    acceptor_geometry: dict | None
    primary: bool
    deferred_to_oxygen: bool = False

    @property
    def priority(self) -> str:
        return "primary" if self.primary else "secondary"

    @property
    def qualifying_primary(self) -> bool:
        return self.primary and max(self.his_donor, self.his_acceptor) >= MIN_PRIMARY_SUPPORT

    def evidence(self) -> dict:
        return {
            **describe(self.atom),
            "distance_A": round(self.distance, 3), "role": self.role,
            "priority": self.priority, "deferred_to_oxygen": self.deferred_to_oxygen,
            "angle_from_His_outward_axis_deg": round(self.axis_angle, 1),
            "supports_His_as_donor": round(self.his_donor, 3),
            "supports_His_as_acceptor": round(self.his_acceptor, 3),
            "His_donor_geometry": self.donor_geometry,
            "His_acceptor_geometry": self.acceptor_geometry,
        }


def describe(atom: SiteAtom) -> dict:
    return {"atom": atom.label, "element": atom.element, "serial": atom.pdb.serial,
            "altloc": atom.residue.altloc, "occupancy": round(atom.pdb.occ, 3)}


def unit(vector: np.ndarray) -> np.ndarray:
    length = np.linalg.norm(vector)
    if length < 1e-8:
        raise ValueError("Degenerate vector")
    return vector / length


def angle(first: np.ndarray, second: np.ndarray) -> float:
    a, b = np.linalg.norm(first), np.linalg.norm(second)
    if a < 1e-8 or b < 1e-8:
        return 180.0
    return degrees(acos(float(np.clip(np.dot(first, second) / (a * b), -1, 1))))


def ring_geometry(target: SiteResidue) -> tuple[dict[str, tuple[SiteAtom, np.ndarray]], float]:
    atoms = [target.atoms[name] for name in RING]
    positions = np.array([atom.position for atom in atoms])
    lengths = np.linalg.norm(positions - np.roll(positions, 1, axis=0), axis=1)
    if np.any((lengths < 1.10) | (lengths > 1.70)):
        raise ValueError("Invalid His ring bond length")
    _, singular, axes = np.linalg.svd(positions - positions.mean(axis=0))
    normal = axes[-1]
    deviation = float(np.max(np.abs((positions - positions.mean(axis=0)) @ normal)))
    if singular[1] < 0.3 or deviation > 0.22:
        raise ValueError(f"His ring plane is unreliable (maximum deviation {deviation:.3f} A)")
    hydrogen_sites = {}
    for index, atom in enumerate(atoms):
        if atom.name not in H_NAMES:
            continue
        if atom.element != ("N" if atom.name in NITROGENS else "C"):
            raise ValueError(f"Invalid element for ring atom {atom.name}")
        incoming = positions[index - 1] - atom.position
        outgoing = positions[(index + 1) % len(atoms)] - atom.position
        direction = -unit(unit(incoming - np.dot(incoming, normal) * normal)
                          + unit(outgoing - np.dot(outgoing, normal) * normal))
        hydrogen_sites[H_NAMES[atom.name]] = (atom, atom.position + (1.01 if atom.name in NITROGENS else 1.08) * direction)
    return hydrogen_sites, deviation


def hydrogens_by_state(sites: dict) -> dict:
    return {state: {name: pair for name, pair in sites.items()
                    if pair[0].name not in NITROGENS or pair[0].name in protons}
            for state, protons in HIS_PROTONS.items()}


def hydrogen_bond(donor: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray) -> tuple[float, dict]:
    separation = float(np.linalg.norm(acceptor - hydrogen))
    dha = angle(donor - hydrogen, acceptor - hydrogen) if separation > 1e-8 else 0.0
    factor = 0.0
    if 1.35 <= separation <= 2.70 and dha >= 120.0:
        factor = min(1.0, (dha - 110.0) / 45.0) * min(1.0, (3.0 - separation) / 0.70)
    return factor, {"H_to_acceptor_A": round(separation, 3), "D_H_A_deg": round(dha, 1)}


def backbone_hydrogen(structure: Structure, atom: SiteAtom) -> np.ndarray | None:
    alpha = atom.residue.atoms.get("CA")
    if atom.name != "N" or atom.residue.pdb.het_flag == "H" or alpha is None:
        return None
    alpha_bond = alpha.position - atom.position
    if not 1.15 <= np.linalg.norm(alpha_bond) <= 1.65:
        return None
    previous = [candidate for candidate, distance in structure.nearest(atom.position, 1.60)
                if candidate.name == "C" and candidate.residue is not atom.residue
                and candidate.residue.id.chain == atom.residue.id.chain and 1.15 <= distance]
    if len(previous) != 1:
        return None
    direction = unit(alpha_bond) + unit(previous[0].position - atom.position)
    if np.linalg.norm(direction) < 0.2:
        return None
    return atom.position - 1.01 * unit(direction)


def evaluate_contact(structure: Structure, acceptor: SiteAtom, trial_h: np.ndarray, atom: SiteAtom, distance: float,
                     role: Role, settings: Settings) -> Contact:
    axis = trial_h - acceptor.position
    axis_angle = angle(axis, atom.position - acceptor.position)
    limit = settings.s_bond_max if atom.element in WEAK_ACCEPTOR_ELEMENTS else settings.hbond_max
    his_donor = his_acceptor = 0.0
    donating, accepting = None, None
    if HEAVY_MIN <= distance <= limit and axis_angle <= 75 and (role.donor or role.acceptor):
        distance_factor = min(1.0, (limit + 0.43 - distance) / 0.95)
        occupancy = min(1.0, max(0.0, atom.pdb.occ))
        if role.acceptor:
            factor, donating = hydrogen_bond(acceptor.position, trial_h, atom.position)
            if (atom.element in WEAK_ACCEPTOR_ELEMENTS and donating["D_H_A_deg"] >= 120
                    and 2.70 < donating["H_to_acceptor_A"] <= 2.90):
                factor = 0.25
            his_donor = factor * distance_factor * occupancy * role.acceptor
        if role.donor:
            hydrogens = attached_hydrogens(atom)
            if hydrogens:
                options = [(h.position, f"input {h.name}", 1.0) for h in hydrogens]
            elif (estimated := backbone_hydrogen(structure, atom)) is not None:
                options = [(estimated, "estimated backbone N-H", 0.92)]
            else:
                options = []
            if options:
                tested = []
                for point, source, weight in options:
                    factor, geometry = hydrogen_bond(atom.position, point, acceptor.position)
                    alignment = angle(axis, point - acceptor.position)
                    tested.append((factor * weight * max(0.0, min(1.0, (70 - alignment) / 45)),
                                   {**geometry, "neighbor_H": source, "acceptor_axis_angle_deg": round(alignment, 1)}))
                factor, accepting = max(tested, key=lambda item: item[0])
                weight = 0.48 if role.label == "water" and hydrogens else role.donor
            else:
                factor = 0.78 * max(0.0, min(1.0, (65 - axis_angle) / 38))
                weight = role.donor
                accepting = {"neighbor_H": "unavailable; heavy-atom geometry only"}
            his_acceptor = factor * weight * distance_factor * occupancy
    return Contact(atom, distance, role.label, axis_angle, his_donor, his_acceptor, donating, accepting,
                   is_primary_role(atom, role))


def prefer_oxygen(entries: list[Contact]) -> list[Contact]:
    if not any(contact.atom.element == "O" and contact.his_donor > 0 for contact in entries):
        return entries
    return [replace(contact, his_donor=0.0, deferred_to_oxygen=True)
            if contact.his_donor > 0 and contact.atom.element in WEAK_ACCEPTOR_ELEMENTS else contact
            for contact in entries]


def coordination(structure: Structure, target: SiteResidue, neighbors: dict) -> tuple[set[str], list[dict], list[str]]:
    confident = set()
    records = []
    extra_links = []
    for name in RING:
        center = target.atoms[name]
        for other, source in structure.declared_partners.get(center, {}).items():
            if other.residue is target or other.pdb.is_hydrogen():
                continue
            if name not in NITROGENS or not other.pdb.element.is_metal:
                extra_links.append(f"{name} to {other.label} ({source})")
            elif not any(other is nearby for nearby, _ in neighbors[name]):
                records.append({"nitrogen": name, **describe(other),
                                "distance_A": round(float(np.linalg.norm(other.position - center.position)), 3),
                                "source": source, "confirmed": False})
    for name, nearby in neighbors.items():
        center = target.atoms[name]
        for atom, distance in nearby:
            if not atom.pdb.element.is_metal:
                continue
            source = structure.declared_partners.get(center, {}).get(atom)
            if distance > METAL_SEARCH and not source:
                continue
            reliable = 1.65 <= distance <= (2.85 if source else STRONG_METAL_LIMIT.get(atom.element, 0))
            reliable = reliable and atom.pdb.occ >= RING_MIN_OCC and center.pdb.occ >= RING_MIN_OCC
            if reliable:
                confident.add(name)
            records.append({"nitrogen": name, **describe(atom), "distance_A": round(distance, 3),
                            "source": source or "distance", "confirmed": reliable})
    return confident, records, extra_links


def clash_contacts(structure: Structure, target: SiteResidue, sites: dict) -> dict[str, list[dict]]:
    removed = ring_hydrogens(target)
    excluded = removed | {target.atoms[name] for name in RING}
    hard_cutoffs = {"H": 1.18, "D": 1.18, "C": 1.65, "N": 1.37, "O": 1.37, "S": 1.57, "SE": 1.60}
    clashes = {}
    for name, (_, point) in sites.items():
        findings = []
        for atom, distance in structure.nearest(point, 1.92, heavy_only=False):
            if atom in excluded:
                continue
            cutoff = hard_cutoffs.get(atom.element, 1.60)
            if distance < cutoff + 0.25:
                findings.append({**describe(atom), "hydrogen": name, "distance_A": round(distance, 3),
                                 "severe": distance < cutoff, "penalty": round(min(2.5, 3 * (cutoff + 0.25 - distance)), 3)})
        clashes[name] = sorted(findings, key=lambda item: item["distance_A"])
    return clashes


def best_contacts(contacts: dict[str, list[Contact]], protons: frozenset[str]) -> tuple[float, list[dict]]:
    choices = []
    for name in NITROGENS:
        donor = name in protons
        entries = contacts[name]
        primary = [contact for contact in entries if contact.qualifying_primary]
        # Only directional primary evidence takes precedence at this nitrogen.
        pairs = [(contact.his_donor if donor else contact.his_acceptor, contact) for contact in primary or entries]
        shortlist = sorted((pair for pair in pairs if pair[0] > 0), key=lambda pair: -pair[0])[:2]
        choices.append([*shortlist, (0.0, None)])
    best_value = -1.0
    best_pair = None
    for first, second in product(*choices):
        # Do not let one water oxygen, or one acceptor atom, establish both sites independently.
        if first[1] is not None and second[1] is not None and first[1].atom is second[1].atom:
            continue
        bonus = PRIMARY_CONTACT_BONUS * sum(contact.qualifying_primary for _, contact in (first, second)
                                            if contact is not None)
        rank_value = first[0] + second[0] + bonus
        if rank_value > best_value:
            best_value, best_pair = rank_value, (first, second)
    evidence = []
    for name, (value, contact) in zip(NITROGENS, best_pair):
        if contact is not None:
            evidence.append({"nitrogen": name, "His_role": "donor" if name in protons else "acceptor",
                             "score": round(value, 3), "neighbor": contact.atom.label,
                             "distance_A": round(contact.distance, 3), "priority": contact.priority,
                             "qualifying_primary": contact.qualifying_primary})
    return sum(value for value, _ in best_pair), evidence


def analyze(structure: Structure, target: SiteResidue, settings: Settings,
            overrides: dict[tuple[str, str], Role]) -> tuple[dict, dict]:
    data = {**target.site_metadata(),
            "status": "unresolved", "recommendation": None, "top_candidate": None,
            "reasons": [], "candidates": [], "contacts": {}, "metal_contacts": [], "carboxylate_contacts": []}
    missing = sorted(set(RING) - target.atoms.keys())
    if missing:
        data.update(reasons=["missing_ring_atoms"], missing_ring_atoms=missing)
        return data, {}
    try:
        sites, deviation = ring_geometry(target)
    except (ValueError, np.linalg.LinAlgError) as error:
        data.update(reasons=["invalid_ring_geometry"], geometry_error=str(error))
        return data, {}
    data["max_ring_plane_deviation_A"] = round(deviation, 3)
    hydrogen_sets = hydrogens_by_state(sites)
    nearby = {name: sorted((item for item in structure.nearest(target.atoms[name].position, settings.search_radius)
                            if item[0].residue is not target), key=lambda item: item[1]) for name in NITROGENS}
    blocked, metal_records, nonstandard_links = coordination(structure, target, nearby)
    data["metal_contacts"] = metal_records
    if nonstandard_links:
        data.update(reasons=["nonstandard_ring_link"], nonstandard_ring_links=nonstandard_links)
        return data, hydrogen_sets
    if len(blocked) > 1:
        data["reasons"] = ["both_nitrogens_coordinate"]
        return data, hydrogen_sets
    # Include ring-to-backbone side-chain contacts within the target residue in scoring.
    for name in NITROGENS:
        for atom in target.atoms.values():
            if atom.name in RING or atom.pdb.is_hydrogen():
                continue
            distance = float(np.linalg.norm(atom.position - target.atoms[name].position))
            if distance <= settings.search_radius:
                nearby[name].append((atom, distance))
        nearby[name].sort(key=lambda item: item[1])
    carboxylates = []
    for name, entries in nearby.items():
        for atom, distance in entries:
            if not HEAVY_MIN <= distance <= settings.carboxylate_max or atom.pdb.occ < RING_MIN_OCC:
                continue
            if is_carboxylate_oxygen(atom, overrides):
                carboxylates.append({"nitrogen": name, **describe(atom), "distance_A": round(distance, 3),
                                     "alternate_conformation": len(atom.residue.options) > 1})
    carboxylates.sort(key=lambda item: (item["alternate_conformation"], item["distance_A"]))
    data["carboxylate_contacts"] = carboxylates
    contacts = {}
    for name in NITROGENS:
        relevant = []
        for atom, distance in nearby[name]:
            role = role_for(atom, overrides)
            if role.label != "neither":
                relevant.append(evaluate_contact(structure, target.atoms[name], sites[N_H_NAMES[name]][1],
                                                 atom, distance, role, settings))
        contacts[name] = prefer_oxygen(relevant)
    data["contacts"] = {name: [contact.evidence() for contact in entries] for name, entries in contacts.items()}
    clashes = clash_contacts(structure, target, sites)
    raw_scores = {}
    usable = {name: [] if name in blocked else entries for name, entries in contacts.items()}
    for state, protons in HIS_PROTONS.items():
        contact_score, support = best_contacts(usable, protons)
        exclusion = "metal_coordination" if blocked & protons else None
        if state == "HIP":
            two_acceptors = len(support) == 2 and all(item["qualifying_primary"] for item in support)
            data["hip_evidence"] = {"two_primary_acceptors": two_acceptors,
                                    "nearby_carboxylate": bool(carboxylates)}
            if exclusion is None and not (two_acceptors or carboxylates):
                exclusion = "insufficient_HIP_acceptors"
        candidate_clashes = [item for h in hydrogen_sets[state] for item in clashes[h]]
        penalty = sum(max((entry["penalty"] for entry in clashes[h]), default=0) for h in hydrogen_sets[state])
        ph_prior = 0.0
        if state == "HIP" and exclusion is None and settings.ph is not None:
            ph_prior = float(np.clip(0.45 * (settings.pka - settings.ph), -1.5, 1.5))
        primary_bonus = PRIMARY_CONTACT_BONUS * sum(item["qualifying_primary"] for item in support)
        geometry_score = contact_score + primary_bonus - penalty
        raw_scores[state] = geometry_score + ph_prior
        data["candidates"].append({"state": state, "allowed": exclusion is None, "exclusion_reason": exclusion,
                                   "contact_score": round(contact_score, 3), "primary_contact_bonus": round(primary_bonus, 3),
                                   "clash_penalty": round(penalty, 3), "geometry_score": round(geometry_score, 3),
                                   "pH_prior": round(ph_prior, 3), "score": round(raw_scores[state], 3),
                                   "supports": support, "clashes": candidate_clashes,
                                   "hydrogens": list(hydrogen_sets[state])})
    data["candidates"].sort(key=lambda candidate: (not candidate["allowed"], -raw_scores[candidate["state"]]))
    compatible = [candidate for candidate in data["candidates"] if candidate["allowed"]]
    ranked = [candidate for candidate in compatible if not any(clash["severe"] for clash in candidate["clashes"])]
    reasons = []
    if target.options and len(target.options) > 1:
        reasons.append("target_alternate_conformation")
    if min(target.atoms[name].pdb.occ for name in RING) < RING_MIN_OCC:
        reasons.append("low_ring_occupancy")
    if any(not record["confirmed"] for record in metal_records):
        reasons.append("possible_metal_contact")
    if not ranked:
        reasons.append("all_candidates_have_clashes")
        if not blocked and not any(data["hip_evidence"].values()):
            reasons.append("insufficient_HIP_acceptors")
    elif blocked:
        top = ranked[0]
        data["top_candidate"] = top["state"]
        if any(record["confirmed"] and record["altloc"] for record in metal_records):
            reasons.append("metal_alternate_conformation")
        if not reasons:
            data.update(status="recommended", recommendation=top["state"], reasons=["metal_coordination"])
            return data, hydrogen_sets
    else:
        hip = next((candidate for candidate in ranked if candidate["state"] == "HIP"), None)
        if carboxylates and hip is not None:
            if carboxylates[0]["alternate_conformation"]:
                reasons.append("carboxylate_alternate_conformation")
            elif not reasons:
                data.update(status="recommended", recommendation="HIP", top_candidate="HIP",
                            reasons=["nearby_carboxylate"])
                return data, hydrogen_sets
        top = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        gap = raw_scores[top["state"]] - raw_scores[second["state"]] if second else None
        data["decision_margin"] = round(gap, 3) if gap is not None else None
        data["top_candidate"] = top["state"] if gap is None or gap > TIE else None
        if gap is not None and gap < MIN_MARGIN:
            reasons.append("close_candidate_scores")
        if top["state"] != "HIP" and max((item["score"] for item in top["supports"]), default=0) < MIN_SUPPORT:
            reasons.append("weak_directional_evidence")
        if any(contact.role == "unknown" and HEAVY_MIN <= contact.distance <= settings.hbond_max and contact.axis_angle <= 55
               for entries in contacts.values() for contact in entries):
            reasons.append("untyped_polar_neighbor")
        used = {item["neighbor"] for item in top["supports"]}
        if any(cont.atom.label in used and len(cont.atom.residue.options) > 1
               for entries in contacts.values() for cont in entries):
            reasons.append("supporting_neighbor_alternate_conformation")
        if not reasons:
            decision = "two_primary_acceptors" if top["state"] == "HIP" else "directional_contacts"
            data.update(status="recommended", recommendation=top["state"], reasons=[decision])
            return data, hydrogen_sets
    data.update(status="review", reasons=reasons)
    return data, hydrogen_sets
