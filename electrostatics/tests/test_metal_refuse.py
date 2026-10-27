from __future__ import annotations

from cluster import refuse_if_unsupported_elements


def test_metal_refused():
    result = refuse_if_unsupported_elements({"C", "N", "O", "H", "Zn"})
    assert result["quantum_refused"] is True
    assert "Zn" in result["elements"]


def test_organic_ok():
    result = refuse_if_unsupported_elements({"C", "N", "O", "H", "S", "Cl"})
    assert result["quantum_refused"] is False
