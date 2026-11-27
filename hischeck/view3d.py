"""Qt pane for the offline 3D viewer and its PNG fallback."""

import json
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QImage, QPixmap
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QLabel, QStackedWidget, QVBoxLayout, QWidget

from .chemistry import is_water
from .i18n import current_language, messages, tr
from .preview import CARBOXYLATE, H_BOND, ION_CONTACT, METAL, WATER_CONTACT, recorded_links, render_image, site_view
from .sidechain import polar_hydrogens
from .structure import NITROGENS, ResidueID, Structure, ring_hydrogens
from .view_bonds import display_bonds
from .view_metals import coordination_pairs

PAGE = Path(__file__).with_name("web") / "viewer.html"
LINK_TYPES = {CARBOXYLATE: "carboxylate", H_BOND: "hydrogen", METAL: "metal", WATER_CONTACT: "water", ION_CONTACT: "cation"}
ATOM_MARK_RADIUS = 5.0


def build_scene(pdb: Path, evidence: dict, state: str | None) -> dict:
    structure = Structure(pdb)
    target = structure.residues[ResidueID(evidence["chain"], evidence["residue_number"], evidence["insertion_code"])]
    ids = {atom.label: atom for atom in structure.atoms}
    # The evidence identifies upstream atoms by label, while PDB export renumbers serials.
    ids.update({f"{target.id.label} {evidence['input_residue']} {atom.name}": atom for atom in target.atoms.values()})
    links = []
    for nitrogen, neighbor, distance, color in recorded_links(evidence, state):
        if nitrogen in target.atoms and neighbor in ids:
            first, second = target.atoms[nitrogen], ids[neighbor]
            links.append({"from": first.pdb.serial, "to": second.pdb.serial, "label": neighbor,
                          "distance": distance, "kind": LINK_TYPES[color], "color": color})
    metals = [atom for atom in structure.heavy_atoms if atom.pdb.element.is_metal]
    recorded = {frozenset((link["from"], link["to"])) for link in links if link["kind"] == "metal"}
    coordination = [{"from": metal.pdb.serial, "to": atom.pdb.serial, "distance": distance}
                    for metal, atom, distance in coordination_pairs(structure, metals)
                    if frozenset((metal.pdb.serial, atom.pdb.serial)) not in recorded]
    center, axes = site_view(target, evidence, state)
    focus = target.functional_atoms
    if target.kind == "HIS":
        centers = [target.atoms[name] for name in NITROGENS if name in target.atoms]
        if "CB" in target.atoms:
            focus.append(target.atoms["CB"])
        hydrogens = sorted(atom.pdb.serial for atom in ring_hydrogens(target))
    else:
        centers = [atom for atom in focus if atom.element in {"O", "N"}]
        bound = polar_hydrogens(target)
        hydrogens = sorted(atom.pdb.serial for atom in bound)
    nearby_atoms = sorted({atom.pdb.serial for center in centers
                           for atom, _ in structure.nearest(center.position, ATOM_MARK_RADIUS)
                           if atom.element != "C" and atom not in centers})
    pdb_text = pdb.read_text(encoding="utf-8", errors="replace")
    return {
        "site": str(pdb.parent.resolve()),
        "pdb": pdb_text,
        "bonds": display_bonds(structure, pdb_text),
        "target": [atom.pdb.serial for atom in target.atoms.values()],
        "focusAtoms": [atom.pdb.serial for atom in focus],
        "hydrogens": hydrogens,
        "labels": [atom.pdb.serial for atom in centers],
        "waters": [atom.pdb.serial for atom in structure.heavy_atoms if is_water(atom.residue)],
        "metals": {atom.pdb.serial: atom.pdb.element.name for atom in metals},
        "coordination": coordination,
        "nearbyAtoms": nearby_atoms,
        "links": links,
        "center": center.tolist(),
        "axes": axes.tolist(),
    }


class StructureView(QStackedWidget):
    def __init__(self):
        super().__init__()
        self.setMinimumHeight(260)
        self.web = QWebEngineView()
        self.web.settings().setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, False)
        self.web.loadFinished.connect(self._page_loaded)
        self.web.renderProcessTerminated.connect(lambda *_: self._renderer_failed(tr("viewer.process_exited")))
        self.addWidget(self.web)
        fallback = QWidget()
        fallback.setStyleSheet("background: #f2f6f9; color: #486275;")
        layout = QVBoxLayout(fallback)
        self.notice = QLabel()
        self.notice.setWordWrap(True)
        self.notice.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.notice)
        self.image_label = QLabel()
        self.image_label.setMinimumSize(1, 1)
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.image_label, 1)
        self.addWidget(fallback)
        self.fallback = fallback
        self.pixmap = QPixmap()
        self.source = None
        self.ready = False
        self.renderer_error = ""
        self.fallback_message = ""
        self.request = 0
        self.setCurrentWidget(self.web)
        self.web.load(QUrl.fromLocalFile(str(PAGE)))

    @staticmethod
    def locale_payload() -> str:
        return json.dumps({"language": current_language(), "text": messages("viewer.")})

    def retranslate(self):
        self.web.page().runJavaScript(f"window.hisViewer?.setLocale({self.locale_payload()})")
        if self.currentWidget() is self.fallback:
            self._show_fallback(self.fallback_message)

    def _page_loaded(self, ok: bool):
        if ok:
            locale = self.locale_payload()
            missing = json.dumps(tr("viewer.script_missing"))
            self.web.page().runJavaScript(
                f"window.hisViewer ? window.hisViewer.initialize({locale}) : {missing}", self._initialized)
        else:
            self._renderer_failed(tr("viewer.load_failed"))

    def _initialized(self, error):
        if error != "":
            self._renderer_failed(error or tr("viewer.no_webgl"))
            return
        self.ready = True
        self.renderer_error = ""
        self.retranslate()
        self._show_scene()

    def _renderer_failed(self, error: str):
        self.ready = False
        self.renderer_error = error
        self._show_fallback(error)

    def clear(self):
        self.request += 1
        self.source = None
        if self.ready:
            self.setCurrentWidget(self.web)
            self.web.page().runJavaScript("window.hisViewer.clear()")
        elif self.renderer_error:
            self._show_fallback(self.renderer_error)

    def load(self, pdb: Path, png: Path, evidence: dict, state: str | None):
        self.source = (pdb, png, evidence, state)
        self._show_scene()

    def _show_scene(self):
        if self.renderer_error:
            self._show_fallback(self.renderer_error)
        elif self.ready and self.source is not None:
            pdb, _, evidence, state = self.source
            self.request += 1
            request = self.request
            try:
                data = json.dumps(build_scene(pdb, evidence, state), allow_nan=False)
            except (OSError, KeyError, TypeError, ValueError, RuntimeError) as error:
                self._show_fallback(str(error))
                return
            self.web.page().runJavaScript(f"window.hisViewer.load({data})", lambda result: self._loaded(request, result))

    def _loaded(self, request: int, error):
        if request != self.request or not self.ready:
            return
        if error == "":
            self.setCurrentWidget(self.web)
        else:
            self._show_fallback(error or tr("viewer.render_failed"))

    def _show_fallback(self, message: str):
        self.fallback_message = message
        self.setCurrentWidget(self.fallback)
        self.pixmap = QPixmap()
        if self.source is not None:
            pdb, png, evidence, state = self.source
            self.pixmap = QPixmap(str(png))
            if self.pixmap.isNull():
                try:
                    picture = render_image(pdb, evidence, state).convert("RGB")
                    image = QImage(picture.tobytes(), picture.width, picture.height, 3 * picture.width,
                                   QImage.Format.Format_RGB888)
                    self.pixmap = QPixmap.fromImage(image.copy())
                except (OSError, KeyError, TypeError, ValueError, RuntimeError):
                    pass
        self.notice.setText((tr("viewer.fallback") if not self.pixmap.isNull() else "") + message)
        self._fit_image()

    def _fit_image(self):
        self.image_label.setPixmap(self.pixmap.scaled(self.image_label.size(), Qt.AspectRatioMode.KeepAspectRatio,
                                                    Qt.TransformationMode.SmoothTransformation)
                                   if not self.pixmap.isNull() else QPixmap())

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.currentWidget() is self.fallback:
            self._fit_image()
