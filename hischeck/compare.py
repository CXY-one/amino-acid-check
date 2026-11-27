"""Compare a completed run with a manually annotated worksheet."""

import argparse
from collections import Counter
import csv
from pathlib import Path
import re

from .chemistry import HIS_ALIAS, HIS_PROTONS
from .i18n import languages, set_language, tr
from .review import result_columns
from .structure import HIS_NAMES, residue_type


def manual_labels(path: Path, sheet: str | None, *, chain: str | None = None) -> list[dict]:
    from openpyxl import load_workbook

    workbook = load_workbook(path, read_only=True, data_only=True)
    records = []
    seen = set()
    try:
        worksheet = workbook[sheet] if sheet else workbook.active
        rows = worksheet.values
        header = next(rows, ())
        columns = {str(label).strip(): index for index, label in enumerate(header) if label}
        review_language = next((language for language in languages() if set(result_columns(language)) <= columns.keys()), None)
        if review_language:
            site_column, name_column, status_column, judgment_column = result_columns(review_language)
        elif {"序号", "名称", "判断"} <= columns.keys():
            name_column, judgment_column = "名称", "判断"
        else:
            raise ValueError(tr("compare.bad_workbook"))
        for row_number, row in enumerate(rows, 2):
            def read(column: str):
                index = columns.get(column)
                return row[index] if index is not None and index < len(row) else None

            if str(read(name_column)).strip().upper() not in HIS_NAMES:
                continue
            if review_language:
                if str(read(status_column)).strip() != tr("status.checked", language=review_language):
                    continue
                site = str(read(site_column)).strip()
                site_chain, separator, number = site.rpartition(":")
                if not separator:
                    raise ValueError(tr("compare.bad_site", row=row_number, site=site))
                if chain is not None and site_chain != (chain or "_"):
                    continue
            else:
                number = read("序号")
                if isinstance(number, float) and number.is_integer():
                    number = int(number)
            match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", str(number).strip())
            state = str(read(judgment_column)).strip().upper()
            state = HIS_ALIAS.get(state, state)
            if not match:
                raise ValueError(tr("compare.bad_number", row=row_number, number=number))
            if state not in HIS_PROTONS:
                raise ValueError(tr("compare.bad_state", row=row_number, state=state))
            key = int(match[1]), match[2]
            if key in seen:
                raise ValueError(tr("compare.duplicate", row=row_number, number=number))
            seen.add(key)
            records.append({"residue_number": key[0], "insertion_code": key[1], "manual_state": state,
                            "manual_note": str(read("依据") or ""), "xlsx_row": row_number})
    finally:
        workbook.close()
    if not records:
        raise ValueError(tr("compare.no_labels"))
    return records


def compare(summary: Path, workbook: Path, chain: str, sheet: str | None) -> list[dict]:
    fields = {"chain", "residue_number", "insertion_code", "status", "recommendation", "top_candidate",
              "exported_states", "reasons", "directory"}
    with summary.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        if not fields.issubset(reader.fieldnames or []):
            raise ValueError(tr("compare.bad_summary"))
        rows = {(int(row["residue_number"]), row["insertion_code"]): row for row in reader
                if row["chain"] == chain and residue_type(row.get("input_residue") or "HIS") == "HIS"}
    result = []
    for record in manual_labels(workbook, sheet, chain=chain):
        key = record["residue_number"], record["insertion_code"]
        row = rows.get(key, {})
        recommendation = row.get("recommendation", "")
        status = row.get("status", "missing")
        exported = row.get("automatic_states", row.get("exported_states", ""))
        retained = record["manual_state"] in exported.split(";")
        if status == "recommended":
            outcome = "recommended_match" if recommendation == record["manual_state"] else "recommended_diff"
        elif status == "review":
            outcome = "review_includes_manual" if retained else "review_excludes_manual"
        else:
            outcome = status
        result.append({"site": f"{chain or '_'}:{key[0]}{key[1]}", **record, "status": status,
                       "recommendation": recommendation, "top_candidate": row.get("top_candidate", ""),
                       "comparison": outcome, "exported_states": exported, "reasons": row.get("reasons", ""),
                       "directory": row.get("directory", "")})
    return result


def main(argv: list[str] | None = None) -> int:
    set_language()
    parser = argparse.ArgumentParser(description=tr("compare.description"))
    parser.add_argument("summary", type=Path)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("--chain", required=True, help=tr("compare.chain"))
    parser.add_argument("--sheet", help=tr("compare.sheet"))
    parser.add_argument("-o", "--output", type=Path, required=True, help=tr("compare.output"))
    args = parser.parse_args(argv)
    try:
        records = compare(args.summary, args.workbook, "" if args.chain == "_" else args.chain, args.sheet)
        with args.output.open("x", encoding="utf-8-sig", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    except (ImportError, KeyError, OSError, ValueError) as error:
        parser.exit(2, f"hischeck.compare: {error}\n")
    counts = Counter(row["comparison"] for row in records)
    print(tr("compare.count", count=len(records)) + " " + ", ".join(f"{outcome}: {count}" for outcome, count in sorted(counts.items())))
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
