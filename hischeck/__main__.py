"""Command-line entry point for PDB analysis and local exports."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from math import isfinite
from tempfile import TemporaryDirectory

from . import __version__
from .analysis import HEAVY_MIN, METAL_SEARCH, Settings, analyze
from .chemistry import Role, load_overrides
from .i18n import set_language, tr
from .propka_estimates import TITRATABLE, estimate_pkas
from .sidechain_analysis import export_sidechain
from .structure import RESIDUE_TYPES, SiteResidue, Structure, export_pdb


DEFAULTS = Settings()
DEFAULT_OUTPUT_RADIUS = 6.0
PROGRESS_PREFIX = "hischeck-progress "


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def export_histidine(source: Structure, neighbors: list[SiteResidue], target: SiteResidue, folder: Path,
                     settings: Settings, overrides: dict[tuple[str, str], Role], all_states: bool) -> dict:
    evidence, hydrogens = analyze(source, target, settings, overrides)
    recommended = evidence["recommendation"]
    automatic = [candidate["state"] for candidate in evidence["candidates"]
                 if candidate["allowed"] and not any(item["severe"] for item in candidate["clashes"])
                 and (not recommended or candidate["state"] == recommended)]
    evidence["exports"] = []
    if all_states:
        evidence["automatic_states"] = automatic
    states = hydrogens if all_states else automatic
    for state in states:
        filename = f"site_{state}.pdb" if all_states or not recommended else "site.pdb"
        export_pdb(source, neighbors, target, folder / filename, state, hydrogens[state])
        evidence["exports"].append({"pdb": filename, "state": state})
    if not automatic:
        export_pdb(source, neighbors, target, folder / "site_input.pdb", None, {})
        evidence["exports"].append({"pdb": "site_input.pdb", "state": None})
    return evidence


def report_progress(args: argparse.Namespace, phase: str, completed: int = 0, total: int = 0) -> None:
    if getattr(args, "progress_json", False):
        print(PROGRESS_PREFIX + json.dumps({"phase": phase, "completed": completed, "total": total}), flush=True)


def analyze_file(args: argparse.Namespace) -> tuple[Path, Counter]:
    report_progress(args, "reading")
    source_path = args.pdb.resolve()
    output = args.output.resolve()
    if output.exists():
        raise ValueError(tr("error.output_exists", output=output))
    source = Structure(source_path, args.model)
    chains = {"" if chain == "_" else chain for chain in args.chain or []}
    missing_chains = chains - {key.chain for key in source.residues}
    if missing_chains:
        raise ValueError(tr("error.missing_chains", chains=', '.join(sorted(c or '_' for c in missing_chains))))
    overrides = load_overrides(args.roles)
    selected_types = set(getattr(args, "residue", None) or ("HIS",))
    sites = [residue for residue in source.residues.values()
             if residue.kind in selected_types and (not chains or residue.id.chain in chains)]
    settings = Settings(
        search_radius=args.search_radius,
        ph=args.ph,
        pka=DEFAULTS.pka if args.pka is None else args.pka,
        hbond_max=getattr(args, "hbond_max", DEFAULTS.hbond_max),
        s_bond_max=getattr(args, "s_bond_max", DEFAULTS.s_bond_max),
        carboxylate_max=getattr(args, "carboxylate_max", DEFAULTS.carboxylate_max),
    )
    render = getattr(args, "render", False)
    all_states = getattr(args, "all_states", False)
    with_propka = getattr(args, "propka", False)
    if render:
        try:
            from .preview import render_site
        except ImportError as error:
            raise RuntimeError(tr("error.no_pillow")) from error
    fields = ("site", "chain", "residue_number", "insertion_code", "input_residue", "status", "recommendation",
              "top_candidate", "decision_margin", "reasons", "exported_states", "neighbor_residues", "directory")
    if all_states:
        fields += ("automatic_states",)
    if with_propka:
        fields += ("propka_pka",)
        if any(site.kind in TITRATABLE for site in sites):
            report_progress(args, "propka")
            pkas, propka_info = estimate_pkas(source)
        else:
            pkas, propka_info = {}, {"status": "no_titratable_targets"}
        propka_info["estimated_targets"] = sum(site.id in pkas for site in sites)
    summary = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".hischeck-", dir=output.parent) as temp:
        staged_output = Path(temp) / "results"
        staged_output.mkdir()
        report_progress(args, "sites", total=len(sites))
        for completed, target in enumerate(sites, start=1):
            neighbors = source.neighborhood(target, args.radius)
            folder = staged_output / target.id.directory_name(target.kind)
            folder.mkdir()
            if target.kind == "HIS":
                evidence = export_histidine(source, neighbors, target, folder, settings, overrides, all_states)
            else:
                evidence = export_sidechain(source, neighbors, target, folder, settings, overrides)
            if with_propka:
                evidence["propka_pka"] = pkas.get(target.id)
            evidence["neighborhood_radius_A"] = args.radius
            evidence["neighborhood"] = [{"site": residue.id.label, "residue": residue.name, "altloc": residue.altloc}
                                        for residue in neighbors]
            if render:
                for export in evidence["exports"]:
                    image = Path(export["pdb"]).with_suffix(".png").name
                    render_site(folder / export["pdb"], evidence, export["state"], folder / image)
                    export["image"] = image
            write_json(folder / "evidence.json", evidence)
            row = {key: evidence.get(key) for key in fields if key in evidence}
            exported = [item["state"] for item in evidence["exports"] if item["state"] is not None]
            row.update(reasons=";".join(evidence["reasons"]), exported_states=";".join(exported),
                       neighbor_residues=len(neighbors) - 1, directory=folder.name)
            if all_states:
                row["automatic_states"] = ";".join(evidence.get("automatic_states", []))
            summary.append(row)
            report_progress(args, "sites", completed, len(sites))
        report_progress(args, "saving")
        counts = Counter(row["status"] for row in summary)
        with (staged_output / "summary.csv").open("w", encoding="utf-8-sig", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fields)
            writer.writeheader()
            writer.writerows(summary)
        write_json(staged_output / "run.json", {
            "hischeck_version": __version__, "input_pdb": str(source_path), "input_sha256": source.sha256,
            "model": source.model_number, "models_in_input": source.model_count, "target_chains": sorted(chains),
            "target_residues": [kind for kind in RESIDUE_TYPES if kind in selected_types],
            "search_radius_A": args.search_radius, "neighborhood_radius_A": args.radius,
            "hbond_max_A": settings.hbond_max, "s_bond_max_A": settings.s_bond_max,
            "carboxylate_max_A": settings.carboxylate_max,
            "render_previews": render, "all_states": all_states, "pH": args.ph,
            "His_reference_pKa": settings.pka if args.ph is not None else None,
            "role_overrides": {f"{residue}:{atom}": role.label for (residue, atom), role in overrides.items()},
            "site_count": len(summary), "status_counts": dict(counts),
            **({"propka": propka_info} if with_propka else {}),
        })
        if output.exists():
            raise ValueError(tr("error.output_now_exists", output=output))
        staged_output.rename(output)
    return output, counts


def main(argv: list[str] | None = None) -> int:
    set_language()
    parser = argparse.ArgumentParser(description=tr("cli.description"))
    parser.add_argument("pdb", type=Path, help=tr("cli.input"))
    parser.add_argument("-o", "--output", required=True, type=Path, help=tr("cli.output"))
    parser.add_argument("--residue", type=str.upper, choices=RESIDUE_TYPES, action="append",
                        help=tr("cli.residue"))
    parser.add_argument("--chain", action="append", help=tr("cli.chain"))
    parser.add_argument("--model", type=int, help=tr("cli.model"))
    parser.add_argument("--radius", type=float, default=DEFAULT_OUTPUT_RADIUS,
                        help=tr("cli.radius"))
    parser.add_argument("--search-radius", type=float, default=DEFAULTS.search_radius,
                        help=tr("cli.search_radius"))
    parser.add_argument("--hbond-max", type=float, default=DEFAULTS.hbond_max,
                        help=tr("cli.hbond_max"))
    parser.add_argument("--s-bond-max", type=float, default=DEFAULTS.s_bond_max,
                        help=tr("cli.s_bond_max"))
    parser.add_argument("--carboxylate-max", type=float, default=DEFAULTS.carboxylate_max,
                        help=tr("cli.carboxylate_max"))
    parser.add_argument("--render", action="store_true", help=tr("cli.render"))
    parser.add_argument("--progress-json", action="store_true", help=tr("cli.progress_json"))
    parser.add_argument("--all-states", action="store_true",
                        help=tr("cli.all_states"))
    parser.add_argument("--propka", action="store_true", help=tr("cli.propka"))
    parser.add_argument("--ph", type=float, help=tr("cli.ph"))
    parser.add_argument("--pka", type=float, help=tr("cli.pka"))
    parser.add_argument("--roles", type=Path, help=tr("cli.roles"))
    parser.add_argument("--version", action="version", version=f"hischeck {__version__}")
    args = parser.parse_args(argv)
    if not isfinite(args.radius) or args.radius <= 0:
        parser.error(tr("cli.bad_radius"))
    if any(not isfinite(value) or value < HEAVY_MIN for value in (args.hbond_max, args.s_bond_max)):
        parser.error(tr("cli.bad_hbond_limits", minimum=HEAVY_MIN))
    if not isfinite(args.carboxylate_max) or args.carboxylate_max < HEAVY_MIN:
        parser.error(tr("cli.bad_carboxylate_limit", minimum=HEAVY_MIN))
    minimum_search = max(args.hbond_max, args.s_bond_max, METAL_SEARCH)
    if not isfinite(args.search_radius) or args.search_radius < minimum_search:
        parser.error(tr("cli.bad_search_radius", minimum=minimum_search))
    if args.pka is not None and args.ph is None:
        parser.error(tr("cli.pka_requires_ph"))
    if any(value is not None and not 0 <= value <= 14 for value in (args.ph, args.pka)):
        parser.error(tr("cli.bad_ph"))
    try:
        output, counts = analyze_file(args)
    except (RuntimeError, ValueError, OSError) as error:
        parser.exit(2, f"hischeck: {error}\n")
    his_only = set(args.residue or ("HIS",)) == {"HIS"}
    kind = tr("cli.kind.his" if his_only else "cli.kind.sites")
    defaults = tr("cli.default_count", count=counts["default"]) if not his_only else ""
    print(tr("cli.result", total=sum(counts.values()), kind=kind, recommended=counts['recommended'], review=counts['review'], unresolved=counts['unresolved'], defaults=defaults))
    print(output / "summary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
