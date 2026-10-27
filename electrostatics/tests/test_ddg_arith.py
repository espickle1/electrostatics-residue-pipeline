from __future__ import annotations

import pandas as pd
import pytest

from aggregate import classical_summary, two_row_summary
from apbs_run import ddg_elec


def test_ddg_elec_arithmetic():
    wt = (-100.0) - (-40.0) - (-30.0)  # -30
    var = (-90.0) - (-40.0) - (-30.0)  # -20
    assert ddg_elec(var, wt) == pytest.approx(10.0)


def test_two_row_summary():
    rows = [
        {"variant": "WT", "frame": "0", "quantity": "ddG_binding_elec", "value": -30.0, "unit": "kJ/mol"},
        {"variant": "A_E73Q", "frame": "0", "quantity": "ddG_binding_elec", "value": -20.0, "unit": "kJ/mol"},
        {"variant": "WT", "frame": "0", "quantity": "dE_int", "value": -5.0, "unit": "eV"},
        {"variant": "A_E73Q", "frame": "0", "quantity": "dE_int", "value": -4.0, "unit": "eV"},
    ]
    summary = two_row_summary(pd.DataFrame(rows))
    assert set(summary["variant"]) == {"WT", "A_E73Q"}
    wt_row = summary.set_index("variant").loc["WT"]
    var_row = summary.set_index("variant").loc["A_E73Q"]
    assert wt_row["dddG_elec_kJ_mol"] == 0.0
    assert var_row["dddG_elec_kJ_mol"] == 10.0
    assert var_row["ddE_int_eV"] == 1.0


def test_classical_summary_ranks_by_abs_dddg():
    rows = [
        {"variant": "WT", "frame": "0", "quantity": "ddG_binding_elec", "value": 0.0, "unit": "kJ/mol"},
        {"variant": "A_SMALL", "frame": "0", "quantity": "ddG_binding_elec", "value": 5.0, "unit": "kJ/mol"},
        {"variant": "A_LARGE", "frame": "0", "quantity": "ddG_binding_elec", "value": -40.0, "unit": "kJ/mol"},
        {"variant": "A_MID", "frame": "0", "quantity": "ddG_binding_elec", "value": 12.0, "unit": "kJ/mol"},
    ]
    summary = classical_summary(pd.DataFrame(rows))
    assert list(summary["variant"]) == ["A_LARGE", "A_MID", "A_SMALL", "WT"]
    assert list(summary["rank"]) == [1, 2, 3, 4]
    assert summary.set_index("variant").loc["A_LARGE", "dddG_elec_kJ_mol"] == -40.0
