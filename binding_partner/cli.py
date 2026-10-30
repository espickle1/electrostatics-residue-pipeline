"""Command line: classify one residue across a list of structure files."""

from __future__ import annotations

import argparse
import dataclasses
import sys

from .analyze import analyze_files, make_reference
from .config import load_config
from .report import write_outputs
from .vote import vote


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="binding_partner",
        description="Distance-only binding-partner classification of one residue.",
    )
    p.add_argument("files", nargs="*", help="structure files (PDB / mmCIF); the reference file "
                   "is not added automatically, list it if it should vote")
    p.add_argument("--ref", required=True, help="reference structure that defines the residue")
    p.add_argument("--chain", required=True, help="chain of the residue in the reference file")
    p.add_argument("--resname", required=True, help="three-letter residue name, e.g. GLU")
    p.add_argument("--resnum", required=True, type=int, help="residue number in the reference file")
    p.add_argument("--icode", default="", help="insertion code, if any")
    p.add_argument("--config", help="YAML file with a partner: block of thresholds")
    p.add_argument("--biounit-glob", action="append", default=[], metavar="PATTERN",
                   help="filename pattern of files that are already a biological unit "
                   "(repeatable); adds to the config's biounit_globs")
    p.add_argument("--out", default="binding_partner_out", help="output directory")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cfg = load_config(args.config)
    if args.biounit_glob:
        cfg = dataclasses.replace(cfg, biounit_globs=cfg.biounit_globs + tuple(args.biounit_glob))
    ref = make_reference(args.ref, args.chain, args.resname, args.resnum, args.icode)
    files = args.files or [args.ref]
    results = analyze_files(files, ref, cfg)
    summary = vote(results, cfg)
    write_outputs(args.out, results, summary)

    for r in results:
        if r.status == "error":
            print(f"error: {r.path}: {r.detail}", file=sys.stderr)
    site = f"{args.chain}:{args.resname.upper()}{args.resnum}{args.icode}"
    rsa = "n/a" if summary.median_rsa is None else f"{summary.median_rsa:.2f}"
    print(f"Residue {site} (reference {args.ref})")
    print(f"Files: {summary.n_files} given, {summary.n_resolved_files} resolved, "
          f"{summary.n_entries} PDB entries")
    print(f"Relative SASA (median over entries): {rsa}")
    print(f"Result: {', '.join(summary.labels) or 'no resolved structure'}"
          + (f" [confidence {summary.confidence}]" if summary.confidence else ""))
    for row in summary.rows:
        name = f"{row.label}:{row.component}" if row.component else row.label
        print(f"  {name}: {row.n_entries}/{row.n_eligible} entries ({row.status})")
    print(f"Wrote files.csv, evidence.csv, consensus.csv, summary.json to {args.out}")
    return 0 if summary.n_resolved_files else 1
