import csv
import json

import pytest
from pdbbuild import A_SEQ, B_SEQ, PdbBuilder

from binding_partner.analyze import make_reference
from binding_partner.cli import main
from binding_partner.config import Config, load_config


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def test_end_to_end_two_entries_and_one_unusable_file(tmp_path, capsys):
    paths = []
    for name in ("one.pdb", "two.pdb"):
        b = PdbBuilder().peptide("A", A_SEQ).peptide("B", B_SEQ, y0=4.0)
        paths.append(str(b.write(tmp_path / name)))
    unrelated = PdbBuilder().peptide("A", B_SEQ).write(tmp_path / "unrelated.pdb")
    out = tmp_path / "out"

    code = main(["--ref", paths[0], "--chain", "A", "--resname", "GLU", "--resnum", "2",
                 "--out", str(out), *paths, str(unrelated)])
    assert code == 0

    summary = json.loads((out / "summary.json").read_text())
    assert summary["labels"] == ["other_protein"]
    assert summary["confidence"] == "ok"
    assert (summary["n_files"], summary["n_resolved_files"], summary["n_entries"]) == (3, 2, 2)

    files = {r["entry_id"]: r for r in read_csv(out / "files.csv")}
    assert files["ONE"]["status"] == files["TWO"]["status"] == "resolved"
    assert files["UNRELATED"]["status"] != "resolved"

    consensus = read_csv(out / "consensus.csv")
    assert [(r["label"], r["status"], r["n_entries"], r["n_eligible"]) for r in consensus] == [
        ("other_protein", "supported", "2", "2")
    ]
    evidence = read_csv(out / "evidence.csv")
    assert {r["label"] for r in evidence} == {"other_protein"}
    assert all(float(r["min_dist"]) <= 4.5 and int(r["n_pairs"]) >= 3 for r in evidence)
    assert "Result: other_protein" in capsys.readouterr().out


def test_single_file_defaults_to_the_reference_and_reports_low_confidence(tmp_path):
    path = PdbBuilder().peptide("A", A_SEQ).write(tmp_path / "solo.pdb")
    out = tmp_path / "out"
    assert main(["--ref", str(path), "--chain", "A", "--resname", "GLU", "--resnum", "2",
                 "--out", str(out)]) == 0
    summary = json.loads((out / "summary.json").read_text())
    assert summary["labels"] == ["none"] and summary["confidence"] == "low"


def test_reference_residue_must_exist_and_match(tmp_path):
    path = PdbBuilder().peptide("A", A_SEQ).write(tmp_path / "r.pdb")
    with pytest.raises(ValueError, match="is GLU, not LYS"):
        make_reference(path, "A", "LYS", 2)
    with pytest.raises(ValueError, match="no residue 99"):
        make_reference(path, "A", "GLU", 99)
    with pytest.raises(ValueError, match="no protein chain 'Z'"):
        make_reference(path, "Z", "GLU", 2)


def test_config_yaml_overrides_defaults(tmp_path):
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(
        "partner:\n  contact_distance: 5.0\n  min_pairs: 4\n  ion_components: [zn, mg]\n"
        "  biounit_globs: ['*_3mer*']\n"
    )
    cfg = load_config(cfg_path)
    assert (cfg.contact_distance, cfg.min_pairs) == (5.0, 4)
    assert cfg.ion_components == frozenset({"ZN", "MG"})
    assert cfg.biounit_globs == ("*_3mer*",)
    assert cfg.same_chain_gap == Config().same_chain_gap


def test_unknown_config_key_is_an_error(tmp_path):
    cfg_path = tmp_path / "bad.yaml"
    cfg_path.write_text("partner:\n  contact_dist: 5.0\n")
    with pytest.raises(ValueError, match="contact_dist"):
        load_config(cfg_path)


def test_example_config_loads_and_matches_defaults():
    from pathlib import Path

    example = Path(__file__).parents[1] / "binding_partner" / "config.example.yaml"
    assert load_config(example) == Config()
