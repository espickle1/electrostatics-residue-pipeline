import dataclasses

import pytest
from conftest import partner_keys
from pdbbuild import A_SEQ, B_SEQ, LIGAND_ATOMS, PdbBuilder, biomt, cryst1

import binding_partner.analyze as analyze_module
from binding_partner.config import Config
from binding_partner.contacts import protein_label
from binding_partner.structure import ChainInfo


def chain_a():
    return PdbBuilder().peptide("A", A_SEQ)


def test_other_protein_for_a_different_sequence(analyze):
    b = chain_a().peptide("B", B_SEQ, y0=4.0)
    result = analyze(b)
    assert result.status == "resolved"
    assert [s.chain for s in result.sites] == ["A"]
    (partner,) = [p for p in result.sites[0].partners if p.label == "other_protein"]
    assert partner.partner_id == "B"
    assert partner.n_pairs >= 3
    assert partner.min_dist <= 4.5


def test_same_protein_for_an_identical_sequence(analyze):
    b = chain_a().peptide("B", A_SEQ, y0=4.0)
    result = analyze(b)
    assert sorted(s.chain for s in result.sites) == ["A", "B"]
    assert all(("same_protein", "") in s.keys() for s in result.sites)
    assert ("other_protein", "") not in partner_keys(result)


def test_sequence_neighbors_are_not_partners(analyze):
    # Residues 1 and 3 touch residue 2 in a single chain; they must not become same_chain.
    result = analyze(chain_a())
    assert result.sites[0].partners == []


@pytest.mark.parametrize(
    "start, expected",
    [(7, False), (8, True)],  # |i - j| <= 5 is local packing; 6 is the first that counts
)
def test_same_chain_needs_a_sequence_gap_above_five(analyze, start, expected):
    b = PdbBuilder().peptide("A", A_SEQ[:6]).peptide("A", ["LEU"], x0=4.31, y0=4.0, start=start)
    keys = partner_keys(analyze(b))
    assert (("same_chain", "") in keys) is expected


@pytest.mark.parametrize(
    "resname, atoms, expected",
    [
        ("FMN", LIGAND_ATOMS, {("ligand", "FMN")}),
        ("NAG", LIGAND_ATOMS, set()),  # glycans never count
        ("SO4", LIGAND_ATOMS, set()),  # crystallization agent
        ("HOH", LIGAND_ATOMS, set()),
        ("ZN", [("ZN", "ZN", (5.77, -5.5, 0.0))], {("ion", "ZN")}),  # one pair is enough
        ("XYZ", LIGAND_ATOMS[:2], set()),  # two pairs is below the ligand minimum
    ],
)
def test_hetero_groups(analyze, resname, atoms, expected):
    b = chain_a().group("A", resname, 101, atoms)
    assert partner_keys(analyze(b)) == expected


def test_ion_minimum_pairs_is_configurable(analyze):
    b = chain_a().group("A", "ZN", 101, [("ZN", "ZN", (5.77, -5.5, 0.0))])
    cfg = dataclasses.replace(Config(), min_pairs_ion=2)
    assert partner_keys(analyze(b, cfg)) == set()


def test_modified_residue_inside_the_chain_is_not_a_ligand(analyze):
    seq = ["ALA", "GLU", "LYS", "MSE", "LEU", "SER", "GLY", "VAL"]
    b = PdbBuilder().peptide("A", seq, het=("MSE",)).peptide("B", B_SEQ, y0=4.0)
    assert ("ligand", "MSE") not in partner_keys(analyze(b))


def test_highest_occupancy_altloc_is_used(analyze):
    def build(occ_a, occ_b):
        b = PdbBuilder().peptide("A", A_SEQ, no_cb={2})
        b.atom("CB", "GLU", "A", 2, (5.77, -1.5, 0.0), "C", occ=occ_a, alt="A")
        b.atom("CB", "GLU", "A", 2, (5.77, 9.0, 0.0), "C", occ=occ_b, alt="B")
        return b.group("A", "FMN", 101, LIGAND_ATOMS)

    # The ligand only reaches the CB of altloc A; altloc B puts it out of range.
    assert partner_keys(analyze(build(0.3, 0.7))) == set()
    assert partner_keys(analyze(build(0.7, 0.3))) == {("ligand", "FMN")}


def test_none_when_exposed_and_buried_when_not(analyze, monkeypatch):
    (site,) = analyze(chain_a()).sites
    assert site.partners == [] and site.rsa > 0.2 and site.baseline == "none"

    monkeypatch.setattr(analyze_module, "relative_sasa", lambda *a, **k: 0.05)
    (site,) = analyze(chain_a()).sites
    assert site.baseline == "buried"


def test_sasa_never_blocks_a_partner_call(analyze, monkeypatch):
    monkeypatch.setattr(analyze_module, "relative_sasa", lambda *a, **k: 0.01)
    b = chain_a().peptide("B", B_SEQ, y0=4.0)
    (site,) = analyze(b).sites
    assert site.rsa == 0.01
    assert [p.label for p in site.partners] == ["other_protein"]
    assert site.baseline == ""


def test_assembly_limits_partners_and_biounit_glob_skips_it(analyze, tmp_path):
    def build():
        return PdbBuilder(biomt("A")).peptide("A", A_SEQ).peptide("B", A_SEQ, y0=4.0)

    in_assembly = analyze(build(), name="pair.pdb")  # REMARK 350 lists chain A only
    assert in_assembly.model_kind == "assembly"
    assert partner_keys(in_assembly) == set()

    flagged = analyze(build(), Config(biounit_globs=("*pair*",)), name="pair.pdb")
    assert flagged.model_kind == "as_is"
    assert ("same_protein", "") in partner_keys(flagged)


def test_symmetry_mate_contact_is_crystal_only(analyze):
    b = PdbBuilder([cryst1("P 2")]).peptide("A", A_SEQ, x0=-3.0)
    result = analyze(b)
    assert partner_keys(result) == set()  # not in the (absent) assembly
    assert [(p.label, p.partner_id) for p in result.crystal_only] == [("same_protein", "A@sym1")]


def test_far_apart_chains_do_not_see_each_other_through_a_fake_cell(analyze):
    b = chain_a().peptide("B", B_SEQ, y0=40.0)
    assert partner_keys(analyze(b)) == set()


def chain_info(entity_id, seq):
    return ChainInfo("X", "protein", [], entity_id, seq)


def test_mmcif_entity_decides_even_when_sequences_match():
    cfg = Config()
    same, other = chain_info("1", "ACDEFG"), chain_info("2", "ACDEFG")
    assert protein_label(same, chain_info("1", "WWWWWW"), True, cfg) == "same_protein"
    assert protein_label(same, other, True, cfg) == "other_protein"


def test_pdb_falls_back_to_sequence_identity():
    cfg = Config()
    q = chain_info("", "ACDEFGHIKL")
    assert protein_label(q, chain_info("", "ACDEFGHIKL"), False, cfg) == "same_protein"
    assert protein_label(q, chain_info("", "ACDEFGHIKW"), False, cfg) == "same_protein"  # 90%
    assert protein_label(q, chain_info("", "ACDEFGHIWW"), False, cfg) == "other_protein"
