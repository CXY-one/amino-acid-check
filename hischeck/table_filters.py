"""Per-column value menus for the residue review table."""

from PySide6.QtCore import QPoint, QSignalBlocker, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox, QHBoxLayout, QHeaderView, QLineEdit, QListWidget, QListWidgetItem, QMenu,
    QPushButton, QTableWidget, QToolButton, QVBoxLayout, QWidget, QWidgetAction,
)

from .i18n import tr


class FilterHeader(QHeaderView):
    filter_requested = Signal(int)

    def __init__(self, table: QTableWidget, labels: tuple[str, ...]):
        super().__init__(Qt.Orientation.Horizontal, table)
        self.labels = labels
        self.buttons = []
        for column, label in enumerate(labels):
            button = QToolButton(self.viewport(), text="▾", objectName="column-filter")
            button.setFixedSize(17, 18)
            button.setFocusPolicy(Qt.FocusPolicy.TabFocus)
            button.clicked.connect(lambda _checked=False, c=column: self.filter_requested.emit(c))
            self.buttons.append(button)
        self.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.setSectionsClickable(True)
        self.setHighlightSections(False)
        self.sectionClicked.connect(self.filter_requested.emit)
        self.sectionResized.connect(self._position_buttons)
        self.geometriesChanged.connect(self._position_buttons)
        table.horizontalScrollBar().valueChanged.connect(self._position_buttons)
        self.update_buttons(set(), False)

    def _position_buttons(self, *_):
        for column, button in enumerate(self.buttons):
            if column < self.count():
                right = self.sectionViewportPosition(column) + self.sectionSize(column) - 4
                button.move(right - button.width(), (self.viewport().height() - button.height()) // 2)

    def update_buttons(self, filtered: set[int], enabled: bool):
        for column, (label, button) in enumerate(zip(self.labels, self.buttons)):
            active = column in filtered
            title = tr("filter.title", column=label) + (tr("filter.active") if active else "")
            button.setEnabled(enabled)
            button.setAccessibleName(title)
            button.setToolTip(title)
            if button.property("filtered") != active:
                button.setProperty("filtered", active)
                button.style().unpolish(button)
                button.style().polish(button)
                button.update()

    def select_values(self, column: int, values: set[str], selected: set[str] | None) -> set[str] | None:
        menu = QMenu(self)
        panel = QWidget(menu)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(10, 9, 10, 9)
        search = QLineEdit(placeholderText=tr("filter.search"))
        layout.addWidget(search)
        select_all = QCheckBox(tr("filter.select_all"))
        layout.addWidget(select_all)
        choices = QListWidget()
        choices.setMinimumWidth(195)
        choices.setFixedHeight(min(224, max(32, 24 * len(values) + 4)))
        choices.setSpacing(2)
        layout.addWidget(choices)
        checked, unchecked = Qt.CheckState.Checked, Qt.CheckState.Unchecked
        items = []
        for value in sorted(values, key=str.casefold):
            item = QListWidgetItem(value, choices)
            item.setCheckState(checked if selected is None or value in selected else unchecked)
            items.append(item)

        def update_select_all():
            shown = [item for item in items if not item.isHidden()]
            count = sum(item.checkState() == checked for item in shown)
            with QSignalBlocker(select_all):
                select_all.setEnabled(bool(shown))
                select_all.setCheckState(unchecked if not count else
                                         checked if count == len(shown) else Qt.CheckState.PartiallyChecked)

        def search_choices(text: str):
            term = text.strip().casefold()
            for item in items:
                item.setHidden(term not in item.text().casefold())
            update_select_all()

        def check_visible():
            state = unchecked if select_all.checkState() == unchecked else checked
            with QSignalBlocker(choices):
                for item in items:
                    if not item.isHidden():
                        item.setCheckState(state)
            update_select_all()

        result = selected

        def apply_filter():
            nonlocal result
            matching = {item.text() for item in items if not item.isHidden() and item.checkState() == checked}
            result = None if matching == values else matching
            menu.close()

        def clear_filter():
            nonlocal result
            result = None
            menu.close()

        search.textChanged.connect(search_choices)
        search.returnPressed.connect(apply_filter)
        select_all.clicked.connect(check_visible)
        choices.itemChanged.connect(update_select_all)
        update_select_all()
        actions = QHBoxLayout()
        actions.addWidget(QPushButton(tr("filter.clear"), clicked=clear_filter))
        actions.addStretch()
        actions.addWidget(QPushButton(tr("common.apply"), clicked=apply_filter))
        layout.addLayout(actions)
        action = QWidgetAction(menu)
        action.setDefaultWidget(panel)
        menu.addAction(action)
        anchor = self.viewport().mapToGlobal(QPoint(self.sectionViewportPosition(column), self.viewport().height()))
        QTimer.singleShot(0, search.setFocus)
        menu.exec(anchor)
        menu.deleteLater()
        return result
