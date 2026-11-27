"""Draw fixed-view PNGs from the exported coordinates and recorded evidence."""

from functools import lru_cache
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from .analysis import H_NAMES
from .i18n import current_language, tr
from .review import state_text
from .sidechain import polar_hydrogens
from .structure import FUNCTIONAL_GROUPS, RING, ResidueID, SiteResidue, Structure, ring_hydrogens

SIZE = (960, 660)
PLOT = (36, 108, 924, 590)
INK = "#193149"
MUTED = "#52677a"
FIELD = "#f3f7fa"
H_BOND = "#007e76"
CARBOXYLATE = "#aa571a"
METAL = "#8854a3"
WATER_CONTACT = "#4a73a7"
ION_CONTACT = "#985861"
CONTACT_COLORS = {"hydrogen": H_BOND, "metal": METAL, "water": WATER_CONTACT,
                  "lysine": ION_CONTACT, "arginine": ION_CONTACT}
ELEMENT = {"C": "#3d5265", "N": "#3275c5", "O": "#d25b64", "H": "#fffefe", "D": "#fffefe",
           "S": "#d1a136", "SE": "#d1a136", "P": "#cf853d"}


@lru_cache
def _font(size: int, bold: bool, language: str):
    latin = ("segoeuib.ttf", "DejaVuSans-Bold.ttf") if bold else ("segoeui.ttf", "DejaVuSans.ttf")
    chinese = ("msyhbd.ttc" if bold else "msyh.ttc", "NotoSansCJK-Regular.ttc") if language.startswith("zh") else ()
    for name in (*chinese, *latin):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default(size=size)


def font(size: int, bold: bool = False):
    return _font(size, bold, current_language())


def site_view(target: SiteResidue, evidence: dict | None = None, state: str | None = None) -> tuple[np.ndarray, np.ndarray]:
    atoms = target.functional_atoms
    points = np.array([atom.position for atom in atoms or list(target.atoms.values())])
    center = points.mean(axis=0)
    if len(atoms) >= 3:
        _, singular, axes = np.linalg.svd(points - center, full_matrices=False)
        if singular[1] > 1e-6:
            normal = axes[-1]
            if target.kind == "HIS":
                left = target.atoms.get("ND1", atoms[0]).position
                right = target.atoms.get("NE2", atoms[-1]).position
                apex = target.atoms.get("CE1", atoms[1]).position
            else:
                apex, left, right = points
                proton = (evidence or {}).get("protonation_oxygen", {})
                if state == "flip" or (state == "protonated" and proton.get("input_atom") != proton.get("output_atom")):
                    left, right = right, left
            horizontal = right - left
            horizontal -= np.dot(horizontal, normal) * normal
            length = np.linalg.norm(horizontal)
            horizontal = horizontal / length if length > 1e-6 else axes[0]
            vertical = np.cross(normal, horizontal)
            if np.dot(vertical, apex - (left + right) / 2) < 0:
                vertical = -vertical
            return center, np.stack((horizontal, vertical))
    return center, np.eye(3)[:2]


def recorded_links(evidence: dict, state: str | None):
    acid_link = None
    if state == "HIP" and "nearby_carboxylate" in evidence.get("reasons", []):
        contact = evidence["carboxylate_contacts"][0]
        acid_link = (contact["nitrogen"], contact["atom"])
        yield (*acid_link, contact["distance_A"], CARBOXYLATE)
    for candidate in evidence["candidates"]:
        if candidate["state"] == state:
            seen = set()
            for support in candidate["supports"]:
                if "center" in support:
                    key = (support["center"], support["interaction"])
                    if key not in seen:
                        seen.add(key)
                        yield support["center"], support["atom"], support["distance_A"], CONTACT_COLORS[support["interaction"]]
                else:
                    link = (support["nitrogen"], support["neighbor"])
                    if link != acid_link:
                        yield (*link, support["distance_A"], H_BOND)
    for record in evidence["metal_contacts"]:
        if record["confirmed"]:
            yield record["nitrogen"], record["atom"], record["distance_A"], METAL


def dashed(draw, start, end, color, width=3):
    delta = np.asarray(end) - start
    length = np.linalg.norm(delta)
    for step in range(0, int(length), 14):
        first = start + step / length * delta
        last = start + min(step + 8, length) / length * delta
        draw.line((*first, *last), fill=color, width=width)


def dot(draw, point, radius, color, outline=None, width=1):
    x, y = point
    draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=color, outline=outline, width=width)


def label(draw, text, point, occupied, color, *, size=19):
    typeface = font(size, True)
    left, top, right, bottom = draw.textbbox((0, 0), text, font=typeface)
    width, height = right - left + 14, bottom - top + 10
    x, y = point
    # Try the side away from the ring first so labels do not cover the bonds.
    side = 1 if x >= SIZE[0] / 2 else -1
    options = [(side, -0.5), (side, -1.15), (side, 0.15), (-side, -0.5),
               (-side, -1.15), (-side, 0.15), (0, -1.65), (0, 0.65)]

    def overlap(box):
        a, b, c, d = box
        return sum(max(0, min(c, right) - max(a, left)) * max(0, min(d, bottom) - max(b, top))
                   for left, top, right, bottom in occupied)

    boxes = []
    for dx, dy in options:
        x0 = x + 15 if dx > 0 else x - 15 - width if dx < 0 else x - width / 2
        x0 = min(max(x0, PLOT[0] + 6), PLOT[2] - width - 6)
        y0 = min(max(y + dy * height, PLOT[1] + 6), PLOT[3] - height - 6)
        boxes.append((x0, y0, x0 + width, y0 + height))
    chosen = min(boxes, key=overlap)
    occupied.append(chosen)
    draw.rounded_rectangle(chosen, radius=6, fill="white", outline="#dae4ec")
    draw.text((chosen[0] + 7 - left, chosen[1] + 5 - top), text, fill=color, font=typeface)


def render_image(pdb_path: Path, evidence: dict, state: str | None) -> Image.Image:
    structure = Structure(pdb_path)
    target = structure.residues[ResidueID(evidence["chain"], evidence["residue_number"], evidence["insertion_code"])]
    center, axes = site_view(target, evidence, state)
    atoms = target.atoms
    group = target.functional_atoms
    is_his = target.kind == "HIS"
    if is_his:
        stem = "CB"
        foreground = set(group) | ring_hydrogens(target)
        bonds = [*zip(RING, RING[1:] + RING[:1]), ("CG", stem), *H_NAMES.items()]
        labels = ("ND1", "NE2", *H_NAMES.values())
    else:
        stem = "CG" if target.kind in {"GLU", "GLN"} else "CB"
        carbon, *polar = FUNCTIONAL_GROUPS[target.kind]
        hydrogens = polar_hydrogens(target)
        foreground = set(group) | hydrogens.keys()
        bonds = [(carbon, atom) for atom in (*polar, stem)]
        bonds.extend((parent.name, hydrogen.name) for hydrogen, parent in hydrogens.items())
        labels = [*polar, *(hydrogen.name for hydrogen in hydrogens)]
    if stem in atoms:
        foreground.add(atoms[stem])
    if not foreground:
        foreground.update(atoms.values())

    by_label = {atom.label: atom for atom in structure.atoms}
    by_label.update({f"{target.id.label} {evidence['input_residue']} {atom.name}": atom
                     for atom in atoms.values()})
    scale_atoms = set(group) or foreground
    preview_states = evidence.get("available_states") or [item["state"] for item in evidence["exports"]]
    scale_atoms = scale_atoms | {by_label[neighbor] for candidate in preview_states
                                 for _, neighbor, _, _ in recorded_links(evidence, candidate) if neighbor in by_label}
    max_extent = np.max(np.abs([(atom.position - center) @ axes.T for atom in scale_atoms]), axis=0)
    scale = min(58.0, 300 / (max_extent[0] + 0.6), 200 / (max_extent[1] + 0.6))
    middle = np.array([(PLOT[0] + PLOT[2]) / 2, (PLOT[1] + PLOT[3]) / 2])
    screen = {atom: (atom.position - center) @ axes.T * (scale, -scale) + middle for atom in structure.atoms}
    links = [(atoms[nitrogen], by_label[neighbor], distance, color)
             for nitrogen, neighbor, distance, color in recorded_links(evidence, state)
             if nitrogen in atoms and neighbor in by_label]
    partners = {atom for _, atom, _, _ in links}

    image = Image.new("RGB", SIZE, "white")
    draw = ImageDraw.Draw(image)
    label_state = state_text(state) if state else tr("png.input")
    draw.text((36, 22), f"{evidence['site']}   /   {label_state}", fill=INK, font=font(33, True))
    description = tr("png.his_plane" if is_his else "png.sidechain")
    draw.text((38, 67), tr("png.local", description=description), fill=MUTED, font=font(20))
    draw.rounded_rectangle(PLOT, radius=17, fill=FIELD, outline="#dce6ee", width=2)
    for atom in structure.heavy_atoms:
        if atom in foreground or atom in partners:
            continue
        x, y = screen[atom]
        if not (PLOT[0] + 12 < x < PLOT[2] - 12 and PLOT[1] + 12 < y < PLOT[3] - 12):
            continue
        if np.linalg.norm(atom.position - center) <= 5.0:
            dot(draw, (x, y), 3, "#b6c6d3")
    for atom in sorted(partners, key=lambda item: item.label):
        if atom.pdb.element.is_metal:
            continue
        for other in atom.residue.atoms.values():
            if other is atom or other in foreground or other.pdb.is_hydrogen():
                continue
            reach = atom.pdb.element.covalent_r + other.pdb.element.covalent_r + 0.25
            x, y = screen[other]
            if (PLOT[0] + 12 < x < PLOT[2] - 12 and PLOT[1] + 12 < y < PLOT[3] - 12
                    and 0.8 <= np.linalg.norm(atom.position - other.position) <= reach):
                draw.line((*screen[atom], x, y), fill="#a0b3c2", width=5)
                dot(draw, (x, y), 5, ELEMENT.get(other.element, "#a0b3c2"))
    for first, second, _, color in links:
        dashed(draw, screen[first], screen[second], color)
    for a, b in bonds:
        if a in atoms and b in atoms:
            first, second = atoms[a], atoms[b]
            if np.linalg.norm(first.position - second.position) <= 1.9:
                draw.line((*screen[first], *screen[second]), fill=INK, width=7)

    occupied = []
    for atom in sorted(foreground | partners, key=lambda item: item.label):
        radius = 8 if atom.pdb.is_hydrogen() else 12 if atom in foreground else 10
        color = METAL if atom.pdb.element.is_metal else ELEMENT.get(atom.element, "#738595")
        dot(draw, screen[atom], radius, color, "#304558" if atom.pdb.is_hydrogen() else "white", width=2)
        x, y = screen[atom]
        occupied.append((x - radius - 2, y - radius - 2, x + radius + 2, y + radius + 2))
    for name in labels:
        if name in atoms:
            label(draw, name, screen[atoms[name]], occupied, INK, size=18)
    for atom in sorted(partners, key=lambda item: item.label):
        distances = [f"{distance:.2f}" for _, other, distance, _ in links if atom is other]
        label(draw, f"{atom.label}  {' / '.join(distances)} Å", screen[atom], occupied,
              next(color for _, other, _, color in links if atom is other))
    legend = [(36, tr("png.hbond"), H_BOND), (255, tr("png.metal"), METAL)] if is_his else []
    if any(color == CARBOXYLATE for _, _, _, color in links):
        legend.append((500, tr("png.carboxylate"), CARBOXYLATE))
    if not is_his:
        captions = ((tr("png.hbond"), H_BOND), (tr("png.water"), WATER_CONTACT), (tr("png.cation"), ION_CONTACT), (tr("png.metal"), METAL))
        legend = [(36 + 165 * i, caption, color) for i, (caption, color) in enumerate(captions)
                  if any(color == link[3] for link in links)]
    for x, caption, color in legend:
        dashed(draw, np.array([x, 624]), np.array([x + 39, 624]), color)
        draw.text((x + 50, 611), caption, fill=MUTED, font=font(19))
    dot(draw, (710, 624), 4, "#b6c6d3")
    draw.text((723, 611), tr("png.nearby"), fill=MUTED, font=font(19))
    return image


def render_site(pdb_path: Path, evidence: dict, state: str | None, png_path: Path) -> None:
    render_image(pdb_path, evidence, state).save(png_path)
