"""Load plain JSON language packs shared by Qt, the CLI and the local 3D page."""

from functools import lru_cache
import json
import os
from pathlib import Path
import re
import sys

DEFAULT_LANGUAGE = "zh-CN"
LANGUAGE_ENV = "AMINO_ACID_CHECK_LANGUAGE"
BUNDLED = Path(__file__).with_name("locales")
LOCALE_ID = re.compile(r"[a-z]{2,3}(?:-[A-Za-z0-9]{2,8})*")
_active = DEFAULT_LANGUAGE


def language_directory() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).parent / "locales"
    return BUNDLED


def _folders() -> tuple[Path, ...]:
    directory = language_directory()
    return (BUNDLED,) if directory == BUNDLED else (directory, BUNDLED)


def _read_pack(language: str) -> dict:
    if not LOCALE_ID.fullmatch(language):
        return {}
    for directory in _folders():
        try:
            data = json.loads((directory / f"{language}.json").read_text(encoding="utf-8-sig"))
        except (OSError, ValueError):
            continue
        if isinstance(data, dict) and isinstance(data.get("name"), str) and isinstance(data.get("messages"), dict):
            return data
    return {}


_pack = lru_cache()(_read_pack)


def languages() -> dict[str, str]:
    candidates = {path.stem for directory in _folders() for path in directory.glob("*.json")}
    return {language: pack["name"] for language in sorted(candidates) if (pack := _read_pack(language))}


def set_language(language: str | None = None) -> str:
    global _active
    _pack.cache_clear()
    requested = language or os.environ.get(LANGUAGE_ENV, DEFAULT_LANGUAGE)
    _active = requested if _pack(requested) else DEFAULT_LANGUAGE
    return _active


def current_language() -> str:
    return _active


def messages(prefix: str, language: str | None = None) -> dict[str, str]:
    defaults = _pack(DEFAULT_LANGUAGE).get("messages", {})
    translated = _pack(language or _active).get("messages", {})
    return {key: translated[key] if isinstance(translated.get(key), str) else value
            for key, value in defaults.items() if key.startswith(prefix)}


def tr(key: str, *, language: str | None = None, default: str | None = None, **values) -> str:
    fallback = _pack(DEFAULT_LANGUAGE).get("messages", {}).get(key)
    if not isinstance(fallback, str):
        fallback = default if default is not None else key
    text = _pack(language or _active).get("messages", {}).get(key)
    if not isinstance(text, str):
        text = fallback
    if not values:
        return text
    try:
        return text.format(**values)
    except (KeyError, ValueError, IndexError, AttributeError, TypeError):
        return fallback.format(**values)
