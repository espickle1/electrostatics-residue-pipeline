from binding_partner.analyze import FileResult, Site
from binding_partner.config import Config
from binding_partner.contacts import Partner
from binding_partner.vote import BELOW_THRESHOLD, SINGLE_ENTRY, SUPPORTED, vote

CFG = Config()


def partner(label, component="", pid="B", pairs=5, dist=3.5):
    return Partner(label, pid, component, pairs, dist, "")


def site(*partners, rsa=0.3):
    return Site("A", "GLU", 2, 1.0, rsa, list(partners), "" if partners else "none")


def entry(eid, *sites, available=()):
    r = FileResult(path=f"{eid}.cif", entry_id=eid, status="resolved", sites=list(sites))
    r.available = set(available)
    return r


def row(summary, label, component=""):
    (match,) = [r for r in summary.rows if r.label == label and r.component == component]
    return match


def test_label_in_two_entries_is_supported():
    results = [entry(e, site(partner("other_protein"))) for e in ("1AAA", "2BBB")]
    summary = vote(results, CFG)
    assert row(summary, "other_protein").status == SUPPORTED
    assert summary.labels == ["other_protein"]
    assert summary.confidence == "ok"


def test_single_entry_is_reported_with_low_confidence():
    summary = vote([entry("1AAA", site(partner("other_protein")))], CFG)
    assert row(summary, "other_protein").status == SINGLE_ENTRY
    assert summary.labels == ["other_protein"] and summary.confidence == "low"


def test_denominator_is_entries_that_could_contain_the_partner():
    fmn = ("ligand", "FMN")
    with_ligand = [
        entry(e, site(partner("ligand", "FMN")), available=[fmn]) for e in ("1AAA", "2BBB")
    ]
    apo = [entry(f"{i}CCC", site(), available=[]) for i in range(6)]
    summary = vote(with_ligand + apo, CFG)
    r = row(summary, "ligand", "FMN")
    assert (r.n_entries, r.n_eligible, r.status) == (2, 2, SUPPORTED)  # 2/8 would fail 30%


def test_both_conditions_must_hold():
    fmn = ("ligand", "FMN")
    # 2 present out of 10 eligible: two entries, but only 20% of those that could show it
    results = [entry(e, site(partner("ligand", "FMN")), available=[fmn]) for e in ("1AAA", "2BBB")]
    results += [entry(f"{i}XXX", site(), available=[fmn]) for i in range(8)]
    assert row(vote(results, CFG), "ligand", "FMN").status == BELOW_THRESHOLD


def test_ligand_identity_must_recur():
    results = [
        entry("1AAA", site(partner("ligand", "FMN")), available=[("ligand", "FMN")]),
        entry("2BBB", site(partner("ligand", "ATP")), available=[("ligand", "ATP")]),
    ]
    summary = vote(results, CFG)
    assert row(summary, "ligand", "FMN").status == BELOW_THRESHOLD
    assert row(summary, "ligand", "ATP").status == BELOW_THRESHOLD
    assert summary.labels == ["none"]  # nothing supported, accessible residue


def test_files_of_one_pdb_entry_count_once():
    ligand_file = entry("1AAA", site(partner("ligand", "FMN")), available=[("ligand", "FMN")])
    apo_file = entry("1AAA", site())  # e.g. a 3mer extract of the same entry
    summary = vote([ligand_file, apo_file], CFG)
    assert summary.n_entries == 1
    assert row(summary, "ligand", "FMN").n_entries == 1
    assert row(summary, "ligand", "FMN").status == SINGLE_ENTRY


def test_label_needs_half_of_the_placed_copies():
    four = [site(partner("same_protein"))] + [site()] * 3
    two = [site(partner("same_protein"))] * 2 + [site()] * 2
    assert vote([entry("1AAA", *four), entry("2BBB", *four)], CFG).labels == ["none"]
    summary = vote([entry("1AAA", *two), entry("2BBB", *two)], CFG)
    assert row(summary, "same_protein").status == SUPPORTED


def test_no_contact_is_none_or_buried_by_median_sasa():
    exposed = vote([entry(e, site(rsa=0.45)) for e in ("1AAA", "2BBB")], CFG)
    assert (exposed.labels, exposed.surface, exposed.median_rsa) == (["none"], True, 0.45)
    buried = vote([entry(e, site(rsa=0.05)) for e in ("1AAA", "2BBB")], CFG)
    assert (buried.labels, buried.surface) == (["buried"], False)


def test_contact_with_low_sasa_keeps_its_label():
    results = [entry(e, site(partner("other_protein"), rsa=0.02)) for e in ("1AAA", "2BBB")]
    summary = vote(results, CFG)
    assert summary.labels == ["other_protein"]
    assert summary.surface is False  # reported, not used as a gate


def test_unresolved_files_do_not_vote():
    dropped = FileResult(path="x.cif", entry_id="9ZZZ", status="position_absent")
    summary = vote([dropped, entry("1AAA", site(partner("other_protein")))], CFG)
    assert (summary.n_files, summary.n_resolved_files, summary.n_entries) == (2, 1, 1)
    assert vote([dropped], CFG).labels == []
