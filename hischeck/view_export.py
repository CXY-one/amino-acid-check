"""Save the interactive renderer's canvas without changing the user's camera."""

from base64 import b64decode
from collections import deque
from dataclasses import dataclass
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from PySide6.QtCore import QObject, QTimer, QUrl, Signal
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView

from .preview import render_site
from .view3d import PAGE, StructureView, build_scene


@dataclass(frozen=True)
class PreviewJob:
    pdb: Path
    png: Path
    evidence: dict
    state: str | None


class PreviewExporter(QObject):
    saved = Signal(object, str)
    finished = Signal()

    def __init__(self, parent):
        super().__init__(parent)
        self.queue = deque()
        self.pending = set()
        self.job = None
        self.web = None
        self.ready = False
        self.failed = False
        self.timeout = QTimer(self, singleShot=True, interval=15000)
        self.timeout.timeout.connect(self._disable_3d)

    def enqueue(self, jobs: list[PreviewJob]):
        for job in jobs:
            if job.png not in self.pending:
                self.pending.add(job.png)
                self.queue.append(job)
        if self.web is None and not self.failed:
            self.web = QWebEngineView(self.parent())
            self.web.hide()
            settings = self.web.settings()
            settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, False)
            self.web.loadFinished.connect(self._page_loaded)
            self.web.renderProcessTerminated.connect(self._disable_3d)
            self.timeout.start()
            self.web.load(QUrl.fromLocalFile(str(PAGE)))
        elif self.job is None:
            QTimer.singleShot(0, self._next)

    def _page_loaded(self, ok: bool):
        if not ok:
            self._disable_3d()
        elif not self.failed:
            self.web.page().runJavaScript(
                f"window.hisViewer?.initialize({StructureView.locale_payload()}, true)", self._initialized)

    def _initialized(self, result):
        if self.failed:
            return
        if result != "":
            self._disable_3d()
            return
        self.timeout.stop()
        self.ready = True
        self._next()

    def _disable_3d(self, *_):
        self.timeout.stop()
        self.failed = True
        if self.job is not None:
            self._save(self.job, None)
        else:
            self._next()

    def _next(self):
        if self.job is not None or not (self.ready or self.failed):
            return
        if not self.queue:
            self.finished.emit()
            return
        job = self.job = self.queue.popleft()
        if self.failed:
            self._save(job, None)
            return
        try:
            scene = json.dumps(build_scene(job.pdb, job.evidence, job.state, snapshot=True), allow_nan=False)
        except (OSError, KeyError, TypeError, ValueError, RuntimeError):
            self._save(job, None)
            return
        self.timeout.start()
        self.web.page().runJavaScript(
            f"JSON.stringify(window.hisViewer.snapshot({scene}))", lambda data: self._save(job, data))

    def _save(self, job: PreviewJob, data):
        if self.job is not job:
            return
        self.timeout.stop()
        kind = ""
        try:
            with TemporaryDirectory(prefix=".hischeck-image-", dir=job.png.parent) as temp:
                staged = Path(temp) / job.png.name
                if data:
                    try:
                        uri = json.loads(data)["png"]
                        prefix = "data:image/png;base64,"
                        if isinstance(uri, str) and uri.startswith(prefix):
                            picture = b64decode(uri[len(prefix):], validate=True)
                            if picture.startswith(b"\x89PNG\r\n\x1a\n"):
                                staged.write_bytes(picture)
                                kind = "3d"
                    except (KeyError, TypeError, ValueError):
                        pass
                if not kind:
                    render_site(job.pdb, job.evidence, job.state, staged)
                    kind = "2d"
                staged.replace(job.png)
        except (OSError, KeyError, TypeError, ValueError, RuntimeError):
            kind = ""
        self.pending.discard(job.png)
        self.job = None
        self.saved.emit(job, kind)
        QTimer.singleShot(0, self._next)
