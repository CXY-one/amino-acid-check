"""Locate and launch the optional desktop PyMOL installation."""

from collections.abc import Iterator
import ctypes
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


WINDOWS_NAMES = ("PyMOLWin.exe", "PyMOL.exe")
NO_CONSOLE = getattr(subprocess, "CREATE_NO_WINDOW", 0)


def pymol_path(value: str | Path | None) -> Path | None:
    """Resolve a full executable path or a PyMOL installation directory."""
    if not value:
        return None
    try:
        path = Path(os.path.expandvars(str(value).strip().strip('"'))).expanduser()
        if not path.is_absolute():
            return None
        if path.is_file():
            valid = path.suffix.lower() == ".exe" if sys.platform == "win32" else os.access(path, os.X_OK)
            return path.resolve() if valid else None
        if sys.platform == "win32" and path.is_dir():
            for subdir in (path, path / "bin", path / "Scripts"):
                for filename in WINDOWS_NAMES:
                    executable = subdir / filename
                    if executable.is_file():
                        return executable.resolve()
    except (OSError, ValueError, RuntimeError):
        pass
    return None


def _registry_paths() -> Iterator[str | Path]:
    import winreg

    def read(key, field):
        try:
            value, _ = winreg.QueryValueEx(key, field)
            return value if isinstance(value, str) else ""
        except OSError:
            return ""

    current_version = r"Software\Microsoft\Windows\CurrentVersion"
    for hive in (winreg.HKEY_CURRENT_USER, winreg.HKEY_LOCAL_MACHINE):
        for view in (winreg.KEY_WOW64_64KEY, winreg.KEY_WOW64_32KEY):
            access = winreg.KEY_READ | view
            for filename in WINDOWS_NAMES:
                try:
                    with winreg.OpenKey(hive, current_version + rf"\App Paths\{filename}", 0, access) as key:
                        if value := read(key, ""):
                            yield value
                except OSError:
                    pass
            try:
                with winreg.OpenKey(hive, current_version + r"\Uninstall", 0, access) as key:
                    for index in range(winreg.QueryInfoKey(key)[0]):
                        try:
                            name = winreg.EnumKey(key, index)
                            with winreg.OpenKey(key, name, 0, access) as entry:
                                if not re.search(r"\bpymol\b", read(entry, "DisplayName"), re.IGNORECASE):
                                    continue
                                if directory := read(entry, "InstallLocation"):
                                    yield directory
                                if icon := read(entry, "DisplayIcon"):
                                    # DisplayIcon can point to the uninstaller and carry an icon index.
                                    icon = re.sub(r",\s*-?\d+\s*$", "", icon.strip()).strip('"')
                                    yield Path(os.path.expandvars(icon)).parent
                        except OSError:
                            continue
            except OSError:
                pass


def _pymol_candidates() -> Iterator[str | Path]:
    if sys.platform != "win32":
        if executable := shutil.which("pymol"):
            yield executable
        return
    yield from _registry_paths()
    for name in WINDOWS_NAMES:
        if executable := shutil.which(name):
            yield executable
    locations = {
        "ProgramData": ("PyMOL",),
        "ProgramFiles": ("PyMOL", "Schrodinger/PyMOL", "Schrödinger/PyMOL"),
        "ProgramFiles(x86)": ("PyMOL", "Schrodinger/PyMOL", "Schrödinger/PyMOL"),
        "LOCALAPPDATA": ("Programs/PyMOL", "PyMOL"),
    }
    for key, folders in locations.items():
        if root := os.environ.get(key):
            for folder in folders:
                yield Path(root) / folder


def find_pymol() -> Path | None:
    for candidate in _pymol_candidates():
        if executable := pymol_path(candidate):
            return executable
    return None


def launch_pymol(executable: Path, pdb: Path) -> subprocess.Popen:
    arguments = [str(executable), str(pdb), "-d", "as sticks, all"]
    options = {"cwd": pdb.parent, "creationflags": NO_CONSOLE}
    if sys.platform != "win32" or not getattr(sys, "frozen", False):
        return subprocess.Popen(arguments, **options)

    # The frozen runtime adds its DLL folder, Qt plugins and QML paths. PyMOL
    # runs with a separate Qt installation, so only pass through other paths.
    bundled = os.path.normcase(os.path.normpath(sys._MEIPASS))
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(path for path in env.get("PATH", "").split(os.pathsep)
                                  if path and os.path.normcase(os.path.normpath(path)) != bundled)
    for key in ("QT_PLUGIN_PATH", "QML2_IMPORT_PATH"):
        value = os.path.normcase(os.path.normpath(env.get(key, "")))
        if value.startswith(bundled + os.sep):
            env.pop(key)

    set_dll_directory = ctypes.windll.kernel32.SetDllDirectoryW
    set_dll_directory.argtypes = [ctypes.c_wchar_p]
    set_dll_directory(None)
    try:
        return subprocess.Popen(arguments, env=env, **options)
    finally:
        set_dll_directory(sys._MEIPASS)
