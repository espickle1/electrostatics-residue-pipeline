"""Combine placed copies, files and PDB entries into one consensus."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from statistics import median

from .analyze import RESOLVED, FileResult, Key
from .config import Config

SUPPORTED, SINGLE_ENTRY, BELOW_THRESHOLD = "supported", "single_entry", "below_threshold"


@dataclass
class ConsensusRow:
    label: str
    component: str
    n_entries: int  # entries where the label is present
    n_eligible: int  # entries that could contain this partner
    fraction: float
    status: str
    confidence: str  # ok / low / empty when below threshold
    partners: str
    min_dist: float | None
    max_pairs: int | None
    entries: str


@dataclass
class Summary:
    n_files: int
    n_resolved_files: int
    n_entries: int
    median_rsa: float | None
    surface: bool | None  # median relative SASA at or above the cutoff
    labels: list[str]  # reported labels, ligand / ion as label:COMPONENT
    confidence: str
    rows: list[ConsensusRow] = field(default_factory=list)


@dataclass
class _Entry:
    present: set[Key] = field(default_factory=set)
    available: set[Key] = field(default_factory=set)
    rsa: list[float] = field(default_factory=list)
    partners: dict[Key, set[str]] = field(default_factory=dict)
    min_dist: dict[Key, float] = field(default_factory=dict)
    max_pairs: dict[Key, int] = field(default_factory=dict)


def _merge_entries(results: list[FileResult], cfg: Config) -> dict[str, _Entry]:
    """One entry per PDB id: a label is present if it is in any of the entry's files."""
    entries: dict[str, _Entry] = {}
    for result in results:
        entry = entries.setdefault(result.entry_id, _Entry())
        copies = len(result.sites)
        counts = Counter(key for site in result.sites for key in site.keys())
        entry.present |= {k for k, n in counts.items() if n / copies >= cfg.copy_fraction}
        entry.available |= result.available
        entry.rsa.extend(s.rsa for s in result.sites if s.rsa is not None)
        for site in result.sites:
            for p in site.partners:
                key = (p.label, p.component)
                entry.partners.setdefault(key, set()).add(p.partner_id)
                entry.min_dist[key] = min(entry.min_dist.get(key, p.min_dist), p.min_dist)
                entry.max_pairs[key] = max(entry.max_pairs.get(key, 0), p.n_pairs)
    return entries


def vote(results: list[FileResult], cfg: Config) -> Summary:
    resolved = [r for r in results if r.status == RESOLVED]
    entries = _merge_entries(resolved, cfg)
    n_entries = len(entries)
    single = n_entries < cfg.min_entries

    rows: list[ConsensusRow] = []
    for key in sorted(set().union(*(e.present for e in entries.values()))):
        present = [eid for eid, e in entries.items() if key in e.present]
        eligible = [eid for eid, e in entries.items() if key in e.available or key in e.present]
        fraction = len(present) / len(eligible)
        if len(present) >= cfg.min_entries and fraction >= cfg.min_entry_fraction:
            status, confidence = SUPPORTED, "ok"
        elif single:
            status, confidence = SINGLE_ENTRY, "low"
        else:
            status, confidence = BELOW_THRESHOLD, ""
        ids = sorted(set().union(*(entries[eid].partners.get(key, set()) for eid in present)))
        dists = [entries[eid].min_dist[key] for eid in present if key in entries[eid].min_dist]
        pairs = [entries[eid].max_pairs[key] for eid in present if key in entries[eid].max_pairs]
        rows.append(
            ConsensusRow(
                key[0], key[1], len(present), len(eligible), round(fraction, 3), status,
                confidence, ";".join(ids), min(dists) if dists else None,
                max(pairs) if pairs else None, ";".join(sorted(present)),
            )
        )
    rank = {SUPPORTED: 0, SINGLE_ENTRY: 1, BELOW_THRESHOLD: 2}
    rows.sort(key=lambda r: (rank[r.status], -r.n_entries, r.label))

    per_entry = [median(e.rsa) for e in entries.values() if e.rsa]
    median_rsa = median(per_entry) if per_entry else None
    surface = None if median_rsa is None else median_rsa >= cfg.surface_rsa_cutoff

    reported = [r for r in rows if r.status in (SUPPORTED, SINGLE_ENTRY)]
    if reported:
        labels = [f"{r.label}:{r.component}" if r.component else r.label for r in reported]
        confidence = "low" if any(r.confidence == "low" for r in reported) else "ok"
    elif n_entries:
        # Nothing in contact: accessibility alone decides none versus buried.
        baseline = "unknown" if surface is None else "none" if surface else "buried"
        status, confidence = (SINGLE_ENTRY, "low") if single else (SUPPORTED, "ok")
        rows.append(
            ConsensusRow(baseline, "", n_entries, n_entries, 1.0, status, confidence, "", None,
                         None, ";".join(sorted(entries)))
        )
        labels = [baseline]
    else:
        labels, confidence = [], ""
    return Summary(
        n_files=len(results),
        n_resolved_files=len(resolved),
        n_entries=n_entries,
        median_rsa=None if median_rsa is None else round(median_rsa, 3),
        surface=surface,
        labels=labels,
        confidence=confidence,
        rows=rows,
    )
