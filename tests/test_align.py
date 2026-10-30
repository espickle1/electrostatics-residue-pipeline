import gemmi
from pdbbuild import A_SEQ, PdbBuilder

from binding_partner.align import Placement, Placer, pair_identity
from binding_partner.analyze import analyze_file, make_reference
from binding_partner.config import Config
from binding_partner.structure import chain_infos, entry_id, load_structure

REF = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def test_place_follows_the_alignment_not_the_numbering():
    placer = Placer(REF, REF.index("K", 10), Config())  # the K of FVKSHF
    target = "GGGG" + REF  # four extra residues at the front shift every index
    out = placer.place(target)
    assert isinstance(out, Placement)
    assert out.index == REF.index("K", 10) + 4
    assert out.identity == 1.0


def test_place_reports_why_a_chain_is_rejected():
    ref_index = REF.index("K", 10)
    placer = Placer(REF, ref_index, Config())
    unrelated = placer.place("WWWWWWWWWWWWWWWWWW")
    assert unrelated in ("low_identity", "position_absent", "residue_mismatch")
    # Residue type differs at the aligned site (K -> E)
    mutated = REF[:ref_index] + "E" + REF[ref_index + 1 :]
    assert placer.place(mutated) == "residue_mismatch"
    # A fragment that stops before the site does not contain it
    assert placer.place(REF[: ref_index - 3]) == "position_absent"


def test_residue_type_check_can_be_switched_off():
    ref_index = REF.index("K", 10)
    mutated = REF[:ref_index] + "E" + REF[ref_index + 1 :]
    cfg = Config(require_same_residue=False)
    assert isinstance(Placer(REF, ref_index, cfg).place(mutated), Placement)


def test_pair_identity_counts_a_fragment_against_its_own_length():
    assert pair_identity("ACDEFGHIKL", "ACDEFGHIKL") == 1.0
    assert pair_identity("ACDEFG", "WWACDEFGWW") == 1.0
    assert pair_identity("", "ACD") == 0.0


def test_modified_residues_stay_in_the_protein_sequence(tmp_path):
    seq = ["ALA", "GLU", "LYS", "MSE", "LEU", "SER", "GLY", "VAL"]
    path = PdbBuilder().peptide("A", seq, het=("MSE",)).write(tmp_path / "m.pdb")
    (info,) = chain_infos(load_structure(path)[0]).values()
    assert info.kind == "protein"
    assert info.seq == "AEKMLSGV"


def test_entry_id_prefers_the_entry_id_over_the_block_name():
    st = gemmi.Structure()
    st.name = "H1_1RU7.cif"  # data block of an extract
    assert entry_id(st, "H1_1RU7_3mer.cif") == "H1_1RU7.CIF"
    st.info["_entry.id"] = "1ru7"
    assert entry_id(st, "H1_1RU7_3mer.cif") == "1RU7"


def test_file_without_the_site_is_dropped_with_a_reason(tmp_path):
    ref_path = PdbBuilder().peptide("A", A_SEQ).write(tmp_path / "ref.pdb")
    other = PdbBuilder().peptide("A", ["TRP", "PHE", "TYR", "HIS", "ARG", "ASP"]).write(
        tmp_path / "other.pdb"
    )
    result = analyze_file(other, make_reference(ref_path, "A", "GLU", 2), Config())
    assert result.status != "resolved"
    assert result.sites == []
