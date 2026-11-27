"""Desktop front end for the hischeck command-line program."""

import argparse
import csv
import json
from math import isfinite
import os
from pathlib import Path
import re
import subprocess
import sys

from PySide6.QtCore import QEvent, QPoint, QSettings, QSignalBlocker, QStandardPaths, QThread, QTimer, Qt, QUrl
from PySide6.QtGui import QBrush, QColor, QDesktopServices, QDragEnterEvent, QDropEvent, QFont, QIcon
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QFileDialog, QFormLayout, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QMessageBox, QProgressBar, QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from .__main__ import DEFAULT_OUTPUT_RADIUS, DEFAULTS, METAL_SEARCH
from .i18n import LANGUAGE_ENV, current_language, language_directory, languages, set_language, tr
from .pymol import find_pymol, launch_pymol, pymol_path
from .review import (
    cancel_review, confirm_review, ensure_variant, export_workbook, read_reviews,
    result_columns, result_status, result_values, review_exports, state_text,
)
from .structure import RESIDUE_TYPES, residue_type
from .table_filters import FilterHeader
from .view3d import StructureView

PROJECT = Path(__file__).resolve().parents[1]
FROZEN = getattr(sys, "frozen", False)
PYMOL_KEY = "pymol/executable"
RESIDUES_KEY = "processing/residues"
LANGUAGE_KEY = "interface/language"
PROPKA_KEY = "processing/propka"
RUN_NOTICE_KEY = "interface/hide_run_notice"
APP_ICON = Path(__file__).with_name("assets") / "app.ico"
PDB_SUFFIXES = {".pdb", ".ent"}
STATUS_COLOR = {"recommended": "#157362", "review": "#965e26", "unresolved": "#9c435b", "checked": "#157362"}
DISTANCE_OPTIONS = (
    ("search-radius", "settings.search_radius", DEFAULTS.search_radius),
    ("radius", "settings.output_radius", DEFAULT_OUTPUT_RADIUS),
    ("hbond-max", "settings.hbond_max", DEFAULTS.hbond_max),
    ("s-bond-max", "settings.s_bond_max", DEFAULTS.s_bond_max),
    ("carboxylate-max", "settings.carboxylate_max", DEFAULTS.carboxylate_max),
)
STYLE = """
HischeckApp, QDialog { background: #f7f9fb; }
QLabel { color: #1d384c; }
QLabel#muted { color: #53697a; }
QLineEdit, QComboBox { background: white; color: #19374a; border: 1px solid #d1dde6; border-radius: 5px; padding: 5px 7px; }
QLineEdit:focus, QComboBox:focus { border-color: #417a90; }
QLineEdit:read-only { background: #fbfcfd; }
QPushButton { background: white; color: #22455a; border: 1px solid #cad8e3; border-radius: 5px; padding: 5px 11px; }
QPushButton:hover { background: #eaf2f7; }
QPushButton:disabled { color: #8093a2; }
QPushButton#run, QPushButton#confirm, QPushButton#confirm-options, QPushButton#cancel-review { background: #17677b; color: white; border-color: #17677b; padding: 6px 18px; }
QPushButton#confirm { border-top-right-radius: 0; border-bottom-right-radius: 0; }
QPushButton#confirm-options { border-top-left-radius: 0; border-bottom-left-radius: 0; padding: 6px 8px; border-left-color: #7ba8b5; }
QPushButton#run:hover, QPushButton#confirm:hover, QPushButton#confirm-options:hover, QPushButton#confirm-options:checked, QPushButton#cancel-review:hover { background: #1a788d; }
QPushButton#run:disabled, QPushButton#confirm:disabled, QPushButton#confirm-options:disabled { background: #8ea6b1; border-color: #8ea6b1; }
QLabel#variant-warning { color: #965e26; }
QTableWidget { background: white; alternate-background-color: #f5f8fb; border: 1px solid #d1dde6; color: #203c50; }
QTableWidget::item:selected { background: #dcedf4; color: #194357; }
QHeaderView::section { background: #ecf2f7; color: #365568; padding: 7px 23px 7px 7px; border: 0; border-bottom: 1px solid #d1dde6; border-right: 1px solid #d1dde6; }
QToolButton#column-filter { color: #365568; background: #f7fafc; border: 1px solid #c6d5e0; border-radius: 3px; padding: 0; }
QToolButton#column-filter:hover, QToolButton#column-filter:focus { border-color: #417a90; background: #e2eef5; }
QToolButton#column-filter[filtered="true"] { color: white; background: #17677b; border-color: #17677b; }
QToolButton#column-filter:disabled { color: #8b9daa; }
QSplitter::handle { background: #e7eff4; }
"""


class AnalysisProcess(QThread):
    """Run the CLI off the GUI thread and capture its exit status and diagnostic output."""

    def __init__(self, command: list[str], parent: QWidget):
        super().__init__(parent)
        self.command = command
        self.error: str | None = None

    def run(self):
        try:
            result = subprocess.run(
                self.command,
                cwd=Path(sys.executable).parent if FROZEN else PROJECT,
                env={**os.environ, LANGUAGE_ENV: current_language()},
                stdin=subprocess.DEVNULL,
                capture_output=True,
                encoding="utf-8",
                errors="replace",
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            if result.returncode:
                self.error = (result.stderr or result.stdout).strip()[-2500:] or tr("gui.worker_exit", code=result.returncode)
        except (OSError, ValueError) as error:
            self.error = str(error)


class LanguagePicker(QComboBox):
    def __init__(self):
        super().__init__()
        self.refresh_languages(current_language())

    def refresh_languages(self, selection: str):
        available = languages()
        if list(available.items()) == [(self.itemData(i), self.itemText(i)) for i in range(self.count())]:
            return
        with QSignalBlocker(self):
            self.clear()
            for code, name in available.items():
                self.addItem(name, code)
            self.setCurrentIndex(max(0, self.findData(selection)))

    def showPopup(self):
        self.refresh_languages(self.currentData() or current_language())
        super().showPopup()


class HischeckApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(tr("gui.title"))
        self.setWindowIcon(QIcon(str(APP_ICON)))
        self.resize(1240, 800)
        self.setMinimumSize(940, 650)
        self.setAcceptDrops(True)
        self.settings = {option: f"{default:g}" for option, _, default in DISTANCE_OPTIONS}
        self.settings.update(chain="", model="", ph="", pka="", roles="")
        self.result_root: Path | None = None
        self.site_dir: Path | None = None
        self.current_pdb: Path | None = None
        self.evidence: dict | None = None
        self.exports: list[dict] = []
        self.rows: list[dict] = []
        self.confirmed: dict[str, str] = {}
        self.column_filters: dict[int, set[str]] = {}
        self.suggested_output: Path | None = None
        self.destination: Path | None = None
        self.process: AnalysisProcess | None = None
        self.preferences = QSettings("HisCheck", "HisCheck")
        self._translations = []
        self._result_message = ("gui.start", None, {})
        self._site_error = None
        self._variant_error = None
        self.with_propka = self.preferences.value(PROPKA_KEY, True, type=bool)
        saved_types = self.preferences.value(RESIDUES_KEY, "", type=str).split(",")
        self.target_types = set(RESIDUE_TYPES).intersection(saved_types) or set(RESIDUE_TYPES)
        documents = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.DocumentsLocation)
        self.default_results = (Path(documents or Path.home()) / "HisCheck") if FROZEN else PROJECT / "results"
        self._build()
        self.clear_site()

    def bind_text(self, setter, key: str):
        setter(tr(key))
        self._translations.append((setter, key))

    def text_label(self, key: str) -> QLabel:
        label = QLabel()
        self.bind_text(label.setText, key)
        return label

    def text_button(self, key: str, callback, parent: QWidget | None = None, **options) -> QPushButton:
        button = QPushButton(parent, clicked=callback, **options)
        self.bind_text(button.setText, key)
        return button

    def _build(self):
        body = QVBoxLayout(self)
        body.setContentsMargins(18, 14, 18, 15)
        body.setSpacing(9)
        header = QHBoxLayout()
        title = self.text_label("gui.title")
        title.setFont(QFont(self.font().family(), 20, QFont.Weight.DemiBold))
        header.addWidget(title)
        header.addStretch()
        self.open_button = self.text_button("gui.open_results", self.browse_results)
        self.settings_button = self.text_button("gui.settings", self.open_settings)
        header.addWidget(self.open_button)
        header.addWidget(self.settings_button)
        body.addLayout(header)
        files = QFormLayout()
        self.input_path = QLineEdit()
        self.bind_text(self.input_path.setPlaceholderText, "gui.drop_pdb")
        self.output_path = QLineEdit()
        for field in (self.input_path, self.output_path):
            field.setAcceptDrops(False)  # The window handles file drops without inserting URL text.
        self.file_controls = [self.input_path, self.output_path]
        for label, edit, text, callback in (
            ("gui.input", self.input_path, "gui.select_file", self.browse_pdb),
            ("gui.output", self.output_path, "gui.select_location", self.browse_output),
        ):
            row = QHBoxLayout()
            row.addWidget(edit, 1)
            button = self.text_button(text, callback)
            row.addWidget(button)
            self.file_controls.append(button)
            files.addRow(self.text_label(label), row)
        body.addLayout(files)
        actions = QHBoxLayout()
        self.run_button = self.text_button("gui.run", self.run, objectName="run")
        self.result_message = QLabel(tr("gui.start"), objectName="muted")
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setFixedWidth(110)
        self.progress.hide()
        actions.addWidget(self.run_button)
        actions.addWidget(self.result_message, 1)
        actions.addWidget(self.progress)
        body.addLayout(actions)

        panes = QSplitter(Qt.Orientation.Horizontal)
        panes.setChildrenCollapsible(False)
        left, right = QWidget(), QWidget()
        left.setMinimumWidth(360)
        panes.addWidget(left)
        panes.addWidget(right)
        panes.setSizes([390, 805])
        panes.setStretchFactor(0, 0)
        panes.setStretchFactor(1, 1)
        body.addWidget(panes, 1)
        results = QVBoxLayout(left)
        results.setContentsMargins(0, 3, 8, 0)
        self.summary = QLabel(tr("gui.results"))
        self.clear_filters_button = self.text_button("gui.clear_filters", self.reset_filters)
        self.clear_filters_button.setEnabled(False)
        result_header = QHBoxLayout()
        result_header.addWidget(self.summary, 1)
        result_header.addWidget(self.clear_filters_button)
        results.addLayout(result_header)
        columns = (*result_columns(), tr("column.pka"))
        self.table = QTableWidget(0, len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.filter_header = FilterHeader(self.table, columns)
        self.table.setHorizontalHeader(self.filter_header)
        self.filter_header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.filter_header.filter_requested.connect(self.filter_column)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(27)
        self.table.setAlternatingRowColors(True)
        self.table.setMouseTracking(True)
        self.table.setShowGrid(False)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.currentCellChanged.connect(self._row_changed)
        results.addWidget(self.table)
        self.root_folder = self.text_button("gui.open_output", lambda: self.show_directory(self.result_root))
        self.root_folder.setEnabled(False)
        self.export_button = self.text_button("gui.export_button", self.save_workbook)
        self.export_button.setEnabled(False)
        result_actions = QHBoxLayout()
        result_actions.addWidget(self.root_folder)
        result_actions.addWidget(self.export_button)
        result_actions.addStretch()
        results.addLayout(result_actions)

        details = QVBoxLayout(right)
        details.setContentsMargins(10, 0, 0, 0)
        candidate_row = QGridLayout()
        candidate_row.setColumnStretch(0, 1)
        candidate_row.setVerticalSpacing(4)
        self.site_title = QLabel()
        self.site_title.setFont(QFont(self.font().family(), 13, QFont.Weight.DemiBold))
        self.candidate = QComboBox()
        self.candidate.setMinimumWidth(105)
        self.candidate.currentIndexChanged.connect(self.load_export)
        self.confirm_button = self.text_button("gui.confirm_next", self.confirm_and_advance, objectName="confirm")
        self.confirm_options = QPushButton("▾", objectName="confirm-options", checkable=True)
        self.bind_text(self.confirm_options.setAccessibleName, "gui.review_actions")
        self.bind_text(self.confirm_options.setToolTip, "gui.review_actions")
        self.cancel_button = self.text_button("gui.cancel_review", self.cancel_current_review,
                                              parent=self, objectName="cancel-review")
        self.cancel_button.hide()
        self.confirm_options.toggled.connect(self.toggle_review_overlay)
        review_buttons = QHBoxLayout()
        review_buttons.setSpacing(0)
        review_buttons.addWidget(self.confirm_button)
        review_buttons.addWidget(self.confirm_options)
        candidate_row.addWidget(self.site_title, 0, 0)
        candidate_row.addWidget(self.text_label("gui.review_choice"), 0, 1)
        candidate_row.addWidget(self.candidate, 0, 2)
        candidate_row.addLayout(review_buttons, 0, 3)
        details.addLayout(candidate_row)
        self.reason = QLabel(objectName="muted")
        self.reason.setWordWrap(True)
        details.addWidget(self.reason)
        self.variant_warning = QLabel(objectName="variant-warning")
        self.variant_warning.setWordWrap(True)
        details.addWidget(self.variant_warning)
        self.viewer = StructureView()
        details.addWidget(self.viewer, 1)
        paths = QFormLayout()
        self.pdb_path = QLineEdit(readOnly=True)
        self.png_path = QLineEdit(readOnly=True)
        self.pdb_path.setAcceptDrops(False)
        self.png_path.setAcceptDrops(False)
        paths.addRow(self.text_label("gui.pdb_path"), self.pdb_path)
        paths.addRow(self.text_label("gui.png_path"), self.png_path)
        details.addLayout(paths)
        external = QHBoxLayout()
        self.pymol_button = self.text_button("gui.open_pymol", self.open_pymol)
        self.site_folder = self.text_button("gui.open_site_directory", lambda: self.show_directory(self.site_dir))
        external.addWidget(self.pymol_button)
        external.addWidget(self.site_folder)
        external.addStretch()
        details.addLayout(external)

    def suggest_output(self, parent: Path):
        stem = Path(self.input_path.text().strip()).stem or "hischeck"
        base = parent / f"{stem}_HisCheck"
        candidate, suffix = base, 2
        while candidate.exists():
            candidate = parent / f"{base.name}_{suffix:02d}"
            suffix += 1
        self.suggested_output = candidate.resolve()
        self.output_path.setText(str(self.suggested_output))

    def set_input_pdb(self, path: Path):
        self.input_path.setText(str(path))
        parent = Path(self.output_path.text()).parent if self.output_path.text() else self.default_results
        self.suggest_output(parent)

    def browse_pdb(self):
        selected, _ = QFileDialog.getOpenFileName(self, tr("gui.select_pdb"), self.input_path.text(), tr("format.pdb"))
        if selected:
            self.set_input_pdb(Path(selected))

    def _dropped_pdb(self, event: QDropEvent) -> Path | None:
        urls = event.mimeData().urls()
        if self.destination is not None or len(urls) != 1 or not urls[0].isLocalFile():
            return None
        path = Path(urls[0].toLocalFile())
        return path if path.suffix.lower() in PDB_SUFFIXES and path.is_file() else None

    def dragEnterEvent(self, event: QDragEnterEvent):
        if self._dropped_pdb(event) is not None:
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent):
        path = self._dropped_pdb(event)
        if path is None:
            event.ignore()
            return
        self.set_input_pdb(path)
        event.acceptProposedAction()

    def browse_output(self):
        selected = QFileDialog.getExistingDirectory(self, tr("gui.choose_output"))
        if selected:
            self.suggest_output(Path(selected))

    def open_settings(self):
        dialog = QDialog(self)
        dialog.setWindowTitle(tr("gui.settings"))
        dialog.setMinimumWidth(620)
        layout = QVBoxLayout(dialog)
        form = QFormLayout()
        residue_row = QHBoxLayout()
        residue_checks = {}
        for kind in RESIDUE_TYPES:
            check = QCheckBox(kind)
            check.setChecked(kind in self.target_types)
            residue_row.addWidget(check)
            residue_checks[kind] = check
        residue_row.addStretch()
        form.addRow(tr("settings.residues"), residue_row)
        language = LanguagePicker()
        language_row = QHBoxLayout()
        language_row.addWidget(language, 1)
        language_row.addWidget(QPushButton(tr("settings.open_languages"),
                                          clicked=lambda: self.open_language_directory(dialog)))
        form.addRow(tr("settings.language"), language_row)
        propka_check = QCheckBox(tr("settings.propka"), checked=self.with_propka)
        form.addRow("PROPKA", propka_check)
        options = [(key, tr(title)) for key, title, _ in DISTANCE_OPTIONS]
        distance_keys = {key for key, _ in options}
        options += [("chain", tr("settings.chains")), ("model", tr("settings.model")), ("ph", tr("settings.ph")),
                    ("pka", tr("settings.reference_pka")), ("roles", tr("settings.roles"))]
        edits = {key: QLineEdit(self.settings[key]) for key, _ in options}

        def browse_roles():
            selected, _ = QFileDialog.getOpenFileName(dialog, tr("settings.select_roles"), edits["roles"].text(), "CSV (*.csv)")
            if selected:
                edits["roles"].setText(selected)

        def reset_distances():
            for key, _, default in DISTANCE_OPTIONS:
                edits[key].setText(f"{default:g}")

        advanced = QWidget()
        advanced_form = QFormLayout(advanced)
        advanced_form.setContentsMargins(0, 0, 0, 0)
        for key, title in options:
            destination = form if key in distance_keys else advanced_form
            if key == "roles":
                field = QHBoxLayout()
                field.addWidget(edits[key])
                field.addWidget(QPushButton(tr("common.browse"), clicked=browse_roles))
                destination.addRow(title, field)
            else:
                destination.addRow(title, edits[key])
        advanced.setVisible(False)
        saved_pymol = self.preferences.value(PYMOL_KEY, "", type=str)
        pymol_edit = QLineEdit(str(pymol_path(saved_pymol) or ""))

        def show_detected_path():
            detected = find_pymol()
            pymol_edit.setPlaceholderText(str(detected) if detected else tr("settings.pymol_placeholder"))
            pymol_edit.setToolTip(str(detected) if detected else "")

        def browse_pymol():
            selected = self.choose_pymol(dialog, pymol_edit.text())
            if selected:
                pymol_edit.setText(selected)

        def use_auto_pymol():
            pymol_edit.clear()
            show_detected_path()

        def apply_settings():
            chosen = [kind for kind, check in residue_checks.items() if check.isChecked()]
            if not chosen:
                QMessageBox.warning(dialog, tr("settings.residue_required"), tr("settings.residue_required_detail"))
                return
            value = pymol_edit.text().strip()
            executable = pymol_path(value)
            if value and executable is None:
                QMessageBox.warning(dialog, tr("settings.check_pymol"), tr("settings.pymol_path_detail"))
                pymol_edit.setFocus()
                return
            if executable:
                self.preferences.setValue(PYMOL_KEY, str(executable))
            else:
                self.preferences.remove(PYMOL_KEY)
            self.target_types = set(chosen)
            self.preferences.setValue(RESIDUES_KEY, ",".join(chosen))
            self.with_propka = propka_check.isChecked()
            self.preferences.setValue(PROPKA_KEY, self.with_propka)
            self.preferences.setValue(LANGUAGE_KEY, language.currentData())
            self.preferences.sync()
            dialog.accept()

        show_detected_path()
        pymol_row = QHBoxLayout()
        pymol_row.addWidget(pymol_edit, 1)
        pymol_row.addWidget(QPushButton(tr("common.browse"), clicked=browse_pymol))
        pymol_row.addWidget(QPushButton(tr("settings.auto_find"), clicked=use_auto_pymol))
        form.addRow(tr("settings.pymol_path"), pymol_row)
        layout.addLayout(form)
        advanced_toggle = QPushButton(checkable=True, objectName="advanced-settings")

        def show_advanced(expanded: bool):
            advanced_toggle.setText(("▾  " if expanded else "▸  ") + tr("settings.advanced"))
            advanced.setVisible(expanded)
            dialog.adjustSize()

        advanced_toggle.toggled.connect(show_advanced)
        layout.addWidget(advanced_toggle, alignment=Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(advanced)
        show_advanced(False)
        help_text = QLabel(tr("settings.help", metal_radius=METAL_SEARCH), objectName="muted", wordWrap=True)
        layout.addWidget(help_text)
        buttons = QHBoxLayout()
        buttons.addWidget(QPushButton(tr("settings.reset_distances"), clicked=reset_distances))
        buttons.addStretch()
        buttons.addWidget(QPushButton(tr("common.cancel"), clicked=dialog.reject))
        buttons.addWidget(QPushButton(tr("common.apply"), clicked=apply_settings))
        layout.addLayout(buttons)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.settings.update({key: edit.text().strip() for key, edit in edits.items()})
            for key, _, default in DISTANCE_OPTIONS:
                self.settings[key] = self.settings[key] or f"{default:g}"
            if language.currentData() != current_language():
                self.apply_language(language.currentData())

    def apply_language(self, language: str):
        set_language(language)
        self.setWindowTitle(tr("gui.title"))
        QApplication.setApplicationDisplayName(tr("gui.title"))
        for setter, key in self._translations:
            setter(tr(key))
        columns = (*result_columns(), tr("column.pka"))
        self.table.setHorizontalHeaderLabels(columns)
        self.filter_header.labels = columns
        with QSignalBlocker(self.table):
            for index in range(len(self.rows)):
                self.update_result_row(index)
        self.refresh_filter_visibility()
        self.refresh_result_message()
        index = self.table.currentRow()
        if index >= 0:
            self.show_review_status(self.rows[index])
            self.show_site_reason()
            for index, entry in enumerate(self.exports):
                state = entry["state"]
                self.candidate.setItemText(index, state_text(state) if state is not None else tr("gui.input_structure"))
            self.show_variant_details()
        else:
            self.site_title.setText(tr("gui.site_details"))
        self.viewer.retranslate()

    def set_result_message(self, key: str, *, suffix: str | None = None, **values):
        self._result_message = (key, suffix, values)
        self.refresh_result_message()

    def refresh_result_message(self):
        key, suffix, values = self._result_message
        if key == "gui.confirmed":
            values = {**values, "state": state_text(values["state"])}
        self.result_message.setText(tr(key, **values) + (tr(suffix) if suffix else ""))

    def toggle_review_overlay(self, visible: bool):
        self.cancel_button.setVisible(visible)
        app = QApplication.instance()
        if visible:
            self.position_review_overlay()
            self.cancel_button.raise_()
            app.installEventFilter(self)
        else:
            app.removeEventFilter(self)

    def position_review_overlay(self):
        if not self.confirm_options.isChecked():
            return
        left = self.confirm_button.mapTo(self, QPoint()).x()
        right = self.confirm_options.mapTo(self, QPoint(self.confirm_options.width(), 0)).x()
        below = self.confirm_options.mapTo(self, QPoint(0, self.confirm_options.height())).y() + 4
        width = max(right - left, self.cancel_button.sizeHint().width())
        self.cancel_button.setGeometry(max(0, right - width), below, width, self.cancel_button.sizeHint().height())

    def eventFilter(self, watched, event):
        if self.confirm_options.isChecked():
            kind = event.type()
            if kind == QEvent.Type.MouseButtonPress:
                position = event.globalPosition().toPoint()
                if not any(widget.rect().contains(widget.mapFromGlobal(position))
                           for widget in (self.cancel_button, self.confirm_options)):
                    self.confirm_options.setChecked(False)
            elif kind == QEvent.Type.KeyPress and event.key() == Qt.Key.Key_Escape:
                self.confirm_options.setChecked(False)
                self.confirm_options.setFocus()
                return True
            elif watched is self and kind == QEvent.Type.WindowDeactivate:
                self.confirm_options.setChecked(False)
            elif watched in (self.confirm_button, self.confirm_options) and kind in (QEvent.Type.Move, QEvent.Type.Resize):
                QTimer.singleShot(0, self.position_review_overlay)
        return super().eventFilter(watched, event)

    def build_command(self, source: Path, destination: Path) -> list[str]:
        if FROZEN:
            command = [str(Path(sys.executable).with_name("hischeck-cli.exe"))]
        else:
            # A source launch through pythonw still uses console Python for the piped worker.
            interpreter = Path(sys.executable)
            if interpreter.name.lower() == "pythonw.exe":
                interpreter = interpreter.with_name("python.exe")
            command = [str(interpreter), "-X", "utf8", "-m", "hischeck"]
        command.extend([str(source), "-o", str(destination), "--render", "--all-states"])
        if self.with_propka:
            command.append("--propka")
        for kind in RESIDUE_TYPES:
            if kind in self.target_types:
                command.extend(["--residue", kind])
        for chain in filter(None, re.split(r"[\s,，]+", self.settings["chain"])):
            command.extend(["--chain", chain])
        for key, value in self.settings.items():
            if key != "chain" and value:
                if key == "roles":
                    value = str(Path(value).expanduser().resolve())
                command.extend([f"--{key}", value])
        return command

    def show_run_notice(self) -> bool:
        if self.preferences.value(RUN_NOTICE_KEY, False, type=bool):
            return True
        notice = QMessageBox(QMessageBox.Icon.Information, tr("notice.title"), tr("notice.body"), parent=self)
        notice.setTextFormat(Qt.TextFormat.PlainText)
        skip = QCheckBox(tr("notice.dont_show"))
        notice.setCheckBox(skip)
        proceed = notice.addButton(tr("gui.run"), QMessageBox.ButtonRole.AcceptRole)
        cancel = notice.addButton(tr("common.cancel"), QMessageBox.ButtonRole.RejectRole)
        notice.setDefaultButton(proceed)
        notice.setEscapeButton(cancel)
        notice.exec()
        if notice.clickedButton() is not proceed:
            return False
        if skip.isChecked():
            self.preferences.setValue(RUN_NOTICE_KEY, True)
            self.preferences.sync()
        return True

    def run(self):
        source = Path(self.input_path.text().strip()).expanduser()
        if not source.is_file():
            QMessageBox.warning(self, tr("gui.pdb_required"), tr("gui.pdb_required_detail"))
            return
        if not self.output_path.text().strip():
            self.suggest_output(self.default_results)
        destination = Path(self.output_path.text().strip()).expanduser().resolve()
        if destination.exists() and destination == self.suggested_output:
            self.suggest_output(destination.parent)
            destination = self.suggested_output
        if destination.exists():
            QMessageBox.warning(self, tr("gui.output_exists"), tr("gui.output_exists_detail"))
            return
        if not self.show_run_notice():
            return
        self.destination = destination
        self.output_path.setText(str(destination))
        self.set_result_message("gui.running")
        self._busy(True)
        process = AnalysisProcess(self.build_command(source.resolve(), destination), self)
        self.process = process
        process.finished.connect(self._run_finished)
        process.start()

    def _busy(self, running: bool):
        for control in (self.run_button, self.settings_button, self.open_button, *self.file_controls):
            control.setEnabled(not running)
        self.progress.setVisible(running)

    def _run_finished(self):
        process, self.process = self.process, None
        destination, self.destination = self.destination, None
        process.deleteLater()
        self._busy(False)
        if process.error is not None:
            self.set_result_message("gui.run_failed_hint")
            QMessageBox.warning(self, tr("gui.run_failed"), process.error)
        elif destination is not None:
            self.load_results(destination)

    def browse_results(self):
        start = self.result_root or (self.default_results if FROZEN else PROJECT / "validation")
        selected = QFileDialog.getExistingDirectory(self, tr("gui.select_result_directory"), str(start))
        if selected:
            self.load_results(Path(selected))

    def load_results(self, directory: Path):
        directory = directory.resolve()
        try:
            with (directory / "summary.csv").open(encoding="utf-8-sig", newline="") as file:
                reader = csv.DictReader(file)
                fields = {"site", "chain", "residue_number", "insertion_code", "status", "recommendation",
                          "top_candidate", "directory"}
                if not fields <= set(reader.fieldnames or []):
                    raise ValueError(tr("gui.missing_summary_columns"))
                rows = list(reader)
            confirmed = read_reviews(directory, rows)
        except (OSError, ValueError, csv.Error) as error:
            QMessageBox.warning(self, tr("gui.read_failed"), str(error))
            return
        self.result_root, self.rows, self.confirmed = directory, rows, confirmed
        self.clear_site()
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.setRowCount(len(rows))
        for index in range(len(rows)):
            self.update_result_row(index)
        self.table.blockSignals(False)
        self.reset_filters()
        message_key = "gui.loaded" if rows else "gui.no_targets"
        try:
            metadata = json.loads((directory / "run.json").read_text(encoding="utf-8"))
            propka_status = metadata.get("propka", {}).get("status")
        except (OSError, ValueError, AttributeError):
            propka_status = None
        suffix = f"gui.propka_{propka_status}" if propka_status in {"unavailable", "failed"} else None
        self.set_result_message(message_key, suffix=suffix, directory=directory.name)
        self.root_folder.setEnabled(True)
        self.export_button.setEnabled(bool(rows))

    def update_result_row(self, index: int):
        row = self.rows[index]
        confirmed = self.confirmed.get(row["directory"])
        status, state = result_status(row, confirmed)
        pka = row.get("propka_pka")
        try:
            estimate = float(pka) if pka else float("nan")
        except ValueError:
            estimate = float("nan")
        values = (*result_values(row, confirmed), f"{estimate:.2f}" if isfinite(estimate) else "—")
        for column, text in enumerate(values):
            item = QTableWidgetItem(text)
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            key = {2: status, 3: state or ""}.get(column, text)
            item.setData(Qt.ItemDataRole.UserRole, key)
            if column == 2 and status in STATUS_COLOR:
                item.setForeground(QBrush(QColor(STATUS_COLOR[status])))
            tooltip = text
            if column == 4:
                tooltip += "\n" + tr("gui.propka_tooltip" if isfinite(estimate) else "gui.propka_tooltip_empty")
            item.setToolTip(tooltip)
            self.table.setItem(index, column, item)

    def filter_column(self, column: int):
        rows = self.rows
        if not rows:
            return
        labels: dict[str, set[str]] = {}
        for index in range(len(rows)):
            item = self.table.item(index, column)
            labels.setdefault(item.text(), set()).add(item.data(Qt.ItemDataRole.UserRole))
        previous = self.column_filters.get(column)
        display_selection = None if previous is None else {label for label, keys in labels.items() if keys & previous}
        selected = self.filter_header.select_values(column, set(labels), display_selection)
        if self.rows is not rows or selected is display_selection:
            return  # The result set changed or the menu was dismissed.
        if selected is None:
            self.column_filters.pop(column, None)
        else:
            self.column_filters[column] = {key for label in selected for key in labels[label]}
        self.apply_filters()

    def reset_filters(self):
        self.column_filters.clear()
        self.apply_filters()

    def refresh_filter_visibility(self):
        for index in range(len(self.rows)):
            hidden = any(self.table.item(index, column).data(Qt.ItemDataRole.UserRole) not in values
                         for column, values in self.column_filters.items())
            self.table.setRowHidden(index, hidden)
        self.filter_header.update_buttons(set(self.column_filters), bool(self.rows))
        self.clear_filters_button.setEnabled(bool(self.column_filters))
        self.update_progress()

    def apply_filters(self):
        self.refresh_filter_visibility()
        current = self.table.currentRow()
        if current < 0 or self.table.isRowHidden(current):
            following = self.next_unchecked(current, default=self.first_visible_row())
            self.table.setCurrentCell(following, 0 if following >= 0 else -1)

    def first_visible_row(self) -> int:
        return next((index for index in range(len(self.rows)) if not self.table.isRowHidden(index)), -1)

    def update_progress(self):
        visible = [row for index, row in enumerate(self.rows) if not self.table.isRowHidden(index)]
        count = sum(row["directory"] in self.confirmed for row in visible)
        text = tr("gui.checked_count", count=count, total=len(visible))
        if len(visible) != len(self.rows):
            text += tr("gui.filtered_count", shown=len(visible), total=len(self.rows))
        self.summary.setText(text)

    def next_unchecked(self, current: int, *, default: int = -1) -> int:
        for indices in (range(current + 1, len(self.rows)), range(max(current, 0))):
            for index in indices:
                if not self.table.isRowHidden(index) and self.rows[index]["directory"] not in self.confirmed:
                    return index
        return default

    def confirm_and_advance(self):
        state = self.candidate.currentData()
        index = self.table.currentRow()
        if state is None or self.current_pdb is None or index < 0:
            return
        row = self.rows[index]
        try:
            confirm_review(self.current_pdb, state, kind=residue_type(row.get("input_residue") or "HIS"))
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, tr("gui.confirm_failed"), str(error))
            return
        self.confirmed[row["directory"]] = state
        self.update_result_row(index)
        self.refresh_filter_visibility()
        following = self.next_unchecked(index)
        if following >= 0:
            self.set_result_message("gui.confirmed", site=row["site"], state=state)
            self.table.setCurrentCell(following, 0)
        else:
            if self.table.isRowHidden(index):
                visible = self.first_visible_row()
                self.table.setCurrentCell(visible, 0 if visible >= 0 else -1)
            else:
                self.show_review_status(row)
            self.set_result_message("gui.all_confirmed" if len(self.confirmed) == len(self.rows) else "gui.filtered_confirmed")

    def cancel_current_review(self):
        index = self.table.currentRow()
        if index < 0 or self.result_root is None:
            return
        row = self.rows[index]
        key = row["directory"]
        if key not in self.confirmed:
            return
        try:
            cancel_review(self.result_root / key)
        except OSError as error:
            QMessageBox.warning(self, tr("gui.cancel_review_failed"), str(error))
            return
        del self.confirmed[key]
        self.update_result_row(index)
        self.apply_filters()
        if self.table.currentRow() == index:
            self.show_review_status(row)
            self.select_result_candidate()
        self.set_result_message("gui.review_cancelled", site=row["site"])

    def show_review_status(self, row: dict):
        state = self.confirmed.get(row["directory"])
        site, kind, status, _ = result_values(row, state)
        self.site_title.setText(f"{site} {kind}  ·  {status}")
        self.confirm_options.setChecked(False)
        self.confirm_options.setEnabled(state is not None)

    def save_workbook(self):
        directory, rows = self.result_root, self.rows
        if directory is None or not rows:
            return
        dialog = QFileDialog(self, tr("gui.export_workbook"), str(directory), tr("format.xlsx"))
        dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        dialog.setDefaultSuffix("xlsx")
        dialog.selectFile(tr("gui.workbook_name"))
        if not dialog.exec():
            return
        path = Path(dialog.selectedFiles()[0])
        if path.suffix.lower() != ".xlsx":
            QMessageBox.warning(self, tr("gui.xlsx_required"), tr("gui.xlsx_required_detail"))
            return
        try:
            confirmed = read_reviews(directory, rows)
            export_workbook(path, rows, confirmed)
        except (ImportError, OSError, ValueError) as error:
            QMessageBox.warning(self, tr("gui.export_failed"), str(error))
            return
        if self.result_root == directory:
            self.set_result_message("gui.workbook_saved", name=path.name, count=len(confirmed))

    def clear_site(self):
        self.site_dir = self.current_pdb = self.evidence = None
        self.exports = []
        self._site_error = self._variant_error = None
        self.site_title.setText(tr("gui.site_details"))
        self.reason.clear()
        self.variant_warning.clear()
        self.variant_warning.hide()
        self.pdb_path.clear()
        self.png_path.clear()
        self.candidate.clear()
        self.candidate.setEnabled(False)
        self.confirm_button.setEnabled(False)
        self.confirm_options.setChecked(False)
        self.confirm_options.setEnabled(False)
        self.site_folder.setEnabled(False)
        self.pymol_button.setEnabled(False)
        self.viewer.clear()

    def _row_changed(self, row: int, _column: int, previous_row: int, _previous_column: int):
        if row != previous_row:
            self.load_site()

    def load_site(self):
        self.clear_site()
        index = self.table.currentRow()
        if index < 0 or self.result_root is None:
            return
        row = self.rows[index]
        self.show_review_status(row)
        folder = self.result_root / row["directory"]
        try:
            data = json.loads((folder / "evidence.json").read_text(encoding="utf-8"))
            exports = review_exports(data)
        except (OSError, ValueError, KeyError, TypeError) as error:
            self._site_error = str(error)
            self.show_site_reason()
            return
        self.site_dir, self.evidence, self.exports = folder, data, exports
        self.show_site_reason()
        self.site_folder.setEnabled(True)
        if exports:
            with QSignalBlocker(self.candidate):
                for entry in exports:
                    option = entry["state"]
                    self.candidate.addItem(state_text(option) if option is not None else tr("gui.input_structure"), option)
            self.candidate.setEnabled(True)
            self.select_result_candidate()

    def show_site_reason(self):
        if self._site_error is not None:
            self.reason.setText(tr("gui.site_failed", error=self._site_error))
            return
        data = self.evidence
        if data is None:
            return
        reason = " · ".join(tr(f"reason.{key}", default=key) for key in data["reasons"])
        hip = next((candidate for candidate in data["candidates"] if candidate["state"] == "HIP"), None)
        if (data["status"] == "review" and hip and hip.get("exclusion_reason") == "insufficient_HIP_acceptors"
                and "insufficient_HIP_acceptors" not in data["reasons"]):
            reason += " · " + tr("reason.insufficient_HIP_acceptors")
        carboxylates = data.get("carboxylate_contacts", [])
        if "nearby_carboxylate" in data["reasons"] and carboxylates:
            contact = carboxylates[0]
            reason += tr("gui.carboxylate_contact", nitrogen=contact['nitrogen'], atom=contact['atom'], distance=contact['distance_A'])
        if exposure := data.get("surface_exposure"):
            reason += tr("gui.probe_exposure", fraction=max(exposure['open_fraction'].values()))
        if missing := data.get("missing_functional_group_atoms"):
            reason += tr("gui.missing_atoms", atoms=", ".join(missing))
        self.reason.setText(reason)

    def select_result_candidate(self):
        if not self.exports or self.evidence is None:
            return
        data = self.evidence
        states = [entry["state"] for entry in self.exports]
        row = self.rows[self.table.currentRow()]
        preferred = self.confirmed.get(row["directory"]) or data.get("recommendation") or data.get("top_candidate")
        if preferred is None and None not in states:
            preferred = (data.get("automatic_states") or [data["exports"][0]["state"]])[0]
        with QSignalBlocker(self.candidate):
            self.candidate.setCurrentIndex(states.index(preferred) if preferred in states else 0)
        self.load_export()

    def load_export(self):
        index = self.candidate.currentIndex()
        if index < 0 or self.evidence is None or self.site_dir is None:
            return
        export = self.exports[index]
        state = export["state"]
        try:
            pdb = ensure_variant(self.site_dir, self.evidence, export)
            if not pdb.is_file():
                raise FileNotFoundError(pdb)
        except (OSError, KeyError, ValueError, RuntimeError) as error:
            self.current_pdb = None
            self.pdb_path.clear()
            self.png_path.clear()
            self.pymol_button.setEnabled(False)
            self.confirm_button.setEnabled(False)
            self._variant_error = str(error)
            self.show_variant_details()
            self.viewer.clear()
            return
        png = self.site_dir / export.get("image", pdb.with_suffix(".png").name)
        self._variant_error = None
        self.current_pdb = pdb
        self.pdb_path.setText(str(pdb))
        self.pymol_button.setEnabled(True)
        self.confirm_button.setEnabled(state is not None)
        self.show_variant_details()
        self.viewer.load(pdb, png, self.evidence, state)

    def show_variant_details(self):
        if self._variant_error is not None:
            self.variant_warning.setText(tr("gui.variant_failed", error=self._variant_error))
            self.variant_warning.show()
            return
        if self.current_pdb is None or self.site_dir is None or self.evidence is None:
            return
        export = self.exports[self.candidate.currentIndex()]
        state = export["state"]
        png = self.site_dir / export.get("image", self.current_pdb.with_suffix(".png").name)
        self.png_path.setText(str(png) if png.is_file() else tr("gui.no_png"))
        selected = next((item for item in self.evidence["candidates"] if item["state"] == state), None)
        warnings = []
        if selected:
            exclusion = selected.get("exclusion_reason")
            if exclusion:
                warnings.append(tr("gui.excluded", state=state, reason=tr(f"reason.{exclusion}", default=exclusion)))
            if any(clash["severe"] for clash in selected["clashes"]):
                warnings.append(tr("gui.severe_clash"))
        elif state is not None and self.evidence["status"] != "default":
            warnings.append(tr("gui.manual_only"))
        self.variant_warning.setText(" · ".join(warnings))
        self.variant_warning.setVisible(bool(warnings))

    def open_language_directory(self, parent: QWidget):
        directory = language_directory()
        try:
            directory.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            QMessageBox.warning(parent, tr("gui.directory_failed"), str(error))
            return
        self.show_directory(directory, parent=parent)

    def show_directory(self, directory: Path | None, *, parent: QWidget | None = None):
        if directory is not None and not QDesktopServices.openUrl(QUrl.fromLocalFile(str(directory))):
            QMessageBox.warning(parent or self, tr("gui.directory_failed"), str(directory))

    def choose_pymol(self, parent: QWidget, start: str = "") -> str:
        types = tr("format.exe") if sys.platform == "win32" else tr("format.all")
        selected, _ = QFileDialog.getOpenFileName(parent, tr("gui.select_pymol"), start, types)
        return selected

    def open_pymol(self):
        pdb = self.current_pdb
        if pdb is None or not pdb.is_file():
            return
        saved = self.preferences.value(PYMOL_KEY, "", type=str)
        executable = pymol_path(saved) or find_pymol()
        selected = False
        if executable is None:
            path = self.choose_pymol(self)
            if not path:
                return
            executable = pymol_path(path)
            if executable is None:
                QMessageBox.warning(self, tr("settings.check_pymol"), tr("gui.pymol_executable_required"))
                return
            selected = True
        try:
            launch_pymol(executable, pdb)
        except OSError as error:
            QMessageBox.warning(self, tr("gui.pymol_failed"), tr("gui.pymol_failed_detail", error=error))
        else:
            if selected:
                self.preferences.setValue(PYMOL_KEY, str(executable))
                self.preferences.sync()

    def closeEvent(self, event):
        if self.destination is not None:
            self.result_message.setText(tr("gui.wait_to_close"))
            event.ignore()
        else:
            super().closeEvent(event)


def main(argv: list[str] | None = None) -> int:
    preferences = QSettings("HisCheck", "HisCheck")
    set_language(os.environ.get(LANGUAGE_ENV) or preferences.value(LANGUAGE_KEY, "", type=str))
    parser = argparse.ArgumentParser(description=tr("gui.description"))
    parser.add_argument("results", nargs="?", type=Path, help=tr("gui.open_result_help"))
    args = parser.parse_args(argv)
    app = QApplication([])
    app.setApplicationDisplayName(tr("gui.title"))
    app.setWindowIcon(QIcon(str(APP_ICON)))
    if sys.platform == "win32":
        app.setFont(QFont("Microsoft YaHei UI", 10))
    app.setStyleSheet(STYLE)
    window = HischeckApp()
    window.show()
    if args.results:
        window.load_results(args.results)
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
