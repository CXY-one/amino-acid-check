"""Asynchronous checks for the latest stable GitHub release."""

import json
import re

from PySide6.QtCore import QObject, QUrl, Signal
from PySide6.QtNetwork import QNetworkAccessManager, QNetworkReply, QNetworkRequest

from . import __version__

REPOSITORY = "CXY-one/amino-acid-check"
LATEST_RELEASE = f"https://api.github.com/repos/{REPOSITORY}/releases/latest"
VERSION_PATTERN = re.compile(r"v?([0-9]{1,9})\.([0-9]{1,9})\.([0-9]{1,9})")


def stable_version(tag: str) -> tuple[int, int, int]:
    match = VERSION_PATTERN.fullmatch(tag) if isinstance(tag, str) else None
    if not match:
        raise ValueError("Expected a stable major.minor.patch version")
    major, minor, patch = map(int, match.groups())
    return major, minor, patch


class UpdateCheck(QObject):
    changed = Signal()

    def __init__(self, parent: QObject):
        super().__init__(parent)
        self.network = QNetworkAccessManager(self)
        self.status = "idle"
        self.latest = ""
        self.release_page = QUrl()

    def check(self):
        if self.status == "checking":
            return
        self.status = "checking"
        self.changed.emit()
        request = QNetworkRequest(QUrl(LATEST_RELEASE))
        request.setTransferTimeout(8000)
        request.setAttribute(QNetworkRequest.Attribute.RedirectPolicyAttribute,
                             QNetworkRequest.RedirectPolicy.SameOriginRedirectPolicy)
        request.setRawHeader(b"Accept", b"application/vnd.github+json")
        request.setRawHeader(b"User-Agent", b"Amino-Acid-Check")
        request.setRawHeader(b"X-GitHub-Api-Version", b"2022-11-28")
        reply = self.network.get(request)
        reply.finished.connect(lambda: self._finished(reply))

    def _finished(self, reply: QNetworkReply):
        try:
            if reply.error() != QNetworkReply.NetworkError.NoError:
                self.status = "failed"
                return
            release = json.loads(bytes(reply.readAll()))
            if not isinstance(release, dict) or release.get("draft") or release.get("prerelease"):
                raise ValueError("Expected a published stable release")
            tag = release.get("tag_name")
            version = stable_version(tag)
            self.latest = ".".join(map(str, version))
            if version > stable_version(__version__):
                self.status = "available"
                self.release_page = QUrl(f"https://github.com/{REPOSITORY}/releases/tag/{tag}")
            else:
                self.status = "current"
                self.release_page = QUrl()
        except ValueError:
            self.status = "failed"
        finally:
            reply.deleteLater()
            self.changed.emit()
