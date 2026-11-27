# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path


from PyInstaller.utils.hooks import collect_data_files, copy_metadata

packaging_dir = Path(SPECPATH).resolve()
project = packaging_dir.parents[1]
package = project / "hischeck"
assets = package / "assets"
app_name = "Amino Acid Check"
entries = ("gui_entry", "cli_entry")
analysis = Analysis(
    [str(packaging_dir / f"{entry}.py") for entry in entries],
    pathex=[str(project)],
    binaries=[],
    datas=[(str(package / folder), f"hischeck/{folder}") for folder in ("web", "locales", "assets")]
          + collect_data_files("propka") + copy_metadata("propka"),
)
pyz = PYZ(analysis.pure)


def build_executable(entry, name, console, **options):
    scripts = [script for script in analysis.scripts if script[0] not in entries or script[0] == entry]
    assert any(script[0] == entry for script in scripts), entry
    return EXE(pyz, scripts, [], name=name, exclude_binaries=True, console=console, upx=False, **options)


app = build_executable("gui_entry", app_name, False, icon=str(assets / "app.ico"))
cli = build_executable("cli_entry", "hischeck-cli", True)
# Release Qt uses the ordinary WebEngine resource packs and snapshots.
resources = Path("PySide6") / "resources"
datas = [entry for entry in analysis.datas
         if not (Path(entry[0]).parent == resources and ".debug." in Path(entry[0]).name)]
release = COLLECT(app, cli, analysis.binaries, datas, name=app_name, upx=False)
