"""Persist human decisions and present results in the desktop table and Excel."""

from collections import defaultdict
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory

from .analysis import hydrogens_by_state, ring_geometry
from .chemistry import HIS_PROTONS
from .i18n import tr
from .sidechain import SIDECHAIN_STATES, write_sidechain_pdb
from .structure import ResidueID, Structure, export_pdb, residue_type

COLUMN_KEYS = ("column.site", "column.residue", "column.status", "column.judgment")


def result_columns(language: str | None = None) -> tuple[str, ...]:
    return tuple(tr(key, language=language) for key in COLUMN_KEYS)


def state_text(state: str | None) -> str:
    return tr(f"state.{state}", default=state) if state else "—"


def result_status(row: dict, confirmed_state: str | None) -> tuple[str, str | None]:
    kind = residue_type(row.get("input_residue") or "HIS")
    status = row["status"]
    state = row["recommendation"] or row["top_candidate"]
    if confirmed_state is not None:
        status, state = "checked", confirmed_state
    elif kind in SIDECHAIN_STATES and status in {"default", "recommended"}:
        _, alternative = SIDECHAIN_STATES[kind]
        status = "review" if state == alternative else "recommended"
    return status, state


def result_values(row: dict, confirmed_state: str | None) -> tuple[str, str, str, str]:
    status, state = result_status(row, confirmed_state)
    return (row["site"], residue_type(row.get("input_residue") or "HIS"),
            tr(f"status.{status}", default=status), state_text(state))


def review_exports(evidence: dict) -> list[dict]:
    """List saved structures and the states that can be generated on selection."""
    exports = evidence["exports"]
    by_state = {export["state"]: export for export in exports}
    kind = residue_type(evidence["input_residue"])
    if kind in SIDECHAIN_STATES:
        if "available_states" in evidence:
            states = evidence["available_states"]
        elif evidence.get("missing_functional_group_atoms"):
            return exports
        else:
            states = SIDECHAIN_STATES[kind]  # Older runs only saved the default structure.
    elif "max_ring_plane_deviation_A" in evidence:
        states = HIS_PROTONS
    else:
        return exports
    return [by_state.get(state, {"pdb": f"site_{state}.pdb", "state": state}) for state in states] + [
        export for export in exports if export["state"] is None
    ]


def ensure_variant(folder: Path, evidence: dict, export: dict) -> Path:
    pdb = folder / export["pdb"]
    state = export["state"]
    if pdb.is_file() or state is None:
        return pdb
    kind = residue_type(evidence["input_residue"])
    if not valid_state(kind, state):
        raise ValueError(tr("error.invalid_state", residue=kind, state=state))
    default = SIDECHAIN_STATES[kind][0] if kind in SIDECHAIN_STATES else None
    base = next((folder / item["pdb"] for item in evidence["exports"]
                 if (kind == "HIS" or item["state"] == default) and (folder / item["pdb"]).is_file()), None)
    if base is None:
        raise FileNotFoundError(tr("review.no_source_pdb"))
    source = Structure(base)
    target = source.residues[ResidueID(evidence["chain"], evidence["residue_number"], evidence["insertion_code"])]
    neighbors = list(source.residues_by_name.values())
    with TemporaryDirectory(prefix=".hischeck-", dir=folder) as temp:
        staged = Path(temp) / pdb.name
        if kind == "HIS":
            sites, _ = ring_geometry(target)
            export_pdb(source, neighbors, target, staged, state, hydrogens_by_state(sites)[state])
        else:
            oxygen = evidence.get("protonation_oxygen", {}).get("input_atom")
            write_sidechain_pdb(source, neighbors, target, staged, state, oxygen)
        staged.replace(pdb)
    return pdb


def valid_state(kind: str, state: str | None) -> bool:
    return state in HIS_PROTONS if kind == "HIS" else state in SIDECHAIN_STATES.get(kind, ())


def read_reviews(directory: Path, rows: list[dict]) -> dict[str, str]:
    confirmed = {}
    for row in rows:
        key = row["directory"]
        path = directory / key / "review.json"
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            continue
        state = data.get("state") if isinstance(data, dict) else None
        kind = residue_type(row.get("input_residue") or "HIS")
        if not isinstance(state, str) or not valid_state(kind, state):
            raise ValueError(tr("review.invalid_record", path=path))
        confirmed[key] = state
    return confirmed


def confirm_review(pdb: Path, state: str, *, kind: str = "HIS") -> None:
    if not valid_state(kind, state) or not pdb.is_file():
        raise ValueError(tr("review.unavailable_variant"))
    with TemporaryDirectory(prefix=".hischeck-", dir=pdb.parent) as temp:
        staged = Path(temp) / "review.json"
        staged.write_text(json.dumps({"state": state, "pdb": pdb.name}, indent=2) + "\n", encoding="utf-8")
        staged.replace(pdb.parent / staged.name)


def cancel_review(folder: Path) -> None:
    (folder / "review.json").unlink(missing_ok=True)


def export_workbook(path: Path, rows: list[dict], confirmed: dict[str, str]) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Font

    by_chain = defaultdict(list)
    for row in rows:
        by_chain[row["chain"]].append(row)
    workbook = Workbook()
    workbook.remove(workbook.active)
    for chain, sites in by_chain.items():
        label = re.sub(r"[\\/*?:\[\]']", "_", chain)
        worksheet = workbook.create_sheet(tr("review.chain", chain=label)[:31] if label else tr("review.no_chain"))
        worksheet.append(result_columns())
        for row in sites:
            worksheet.append(result_values(row, confirmed.get(row["directory"])))
        worksheet.auto_filter.ref = worksheet.dimensions
        worksheet.freeze_panes = "A2"
        for column in ("A", "B", "C", "D"):
            worksheet.column_dimensions[column].width = 15
        for cells in worksheet.iter_rows():
            for cell in cells:
                cell.data_type = "s"
                cell.font = Font(name="Microsoft YaHei", size=11, bold=cell.row == 1)
    with TemporaryDirectory(prefix=".hischeck-", dir=path.parent) as temp:
        staged = Path(temp) / path.name
        workbook.save(staged)
        staged.replace(path)
