from __future__ import annotations

from pathlib import Path

import pytest

from apbs_run import assert_identical_grids, extract_grid_block, parse_elec_energy

FIXTURES = Path(__file__).parent / "fixtures"


def test_extract_grid_block():
    text = (FIXTURES / "apbs_complex.in").read_text(encoding="utf-8")
    grid = extract_grid_block(text)
    assert grid["dime"] == (97, 97, 97)
    assert grid["glen"] == (50.0, 50.0, 50.0)
    assert grid["gcent"] == (1.0, 2.0, 3.0)


def test_identical_grids_pass():
    inputs = {
        name: (FIXTURES / f"apbs_{name}.in").read_text(encoding="utf-8")
        for name in ("complex", "protein", "partner")
    }
    grid = assert_identical_grids(inputs)
    assert grid["dime"] == (97, 97, 97)


def test_identical_grids_fail():
    inputs = {
        "complex": (FIXTURES / "apbs_complex.in").read_text(encoding="utf-8"),
        "protein": (FIXTURES / "apbs_protein.in").read_text(encoding="utf-8"),
        "partner": (FIXTURES / "apbs_partner_bad.in").read_text(encoding="utf-8"),
    }
    with pytest.raises(AssertionError, match="grid mismatch"):
        assert_identical_grids(inputs)


def test_parse_elec_energy():
    text = (FIXTURES / "apbs_sample.log").read_text(encoding="utf-8")
    assert parse_elec_energy(text) == pytest.approx(-1234.5)
