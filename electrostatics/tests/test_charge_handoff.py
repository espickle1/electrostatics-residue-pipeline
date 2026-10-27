from __future__ import annotations

import pytest

from cluster import assert_charge_handoff


def test_charge_handoff_ok():
    assert_charge_handoff(-2.0, -2.0001, tol=1e-3)


def test_charge_handoff_fail():
    with pytest.raises(AssertionError, match="Charge handoff failed"):
        assert_charge_handoff(-2.0, -1.0, tol=1e-3)
