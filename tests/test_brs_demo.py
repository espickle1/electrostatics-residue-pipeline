"""Checks on the real 1BRS file; skipped unless BRS_PDB points at a copy of it."""

import dataclasses
import os

import pytest

from binding_partner.analyze import analyze_file, make_reference
from binding_partner.config import Config

BRS = os.environ.get("BRS_PDB")
pytestmark = pytest.mark.skipif(not BRS, reason="set BRS_PDB to a 1BRS PDB file")


def classify(resname, resnum, cfg=None):
    ref = make_reference(BRS, "A", resname, resnum)
    result = analyze_file(BRS, ref, cfg or Config())
    (site,) = result.sites  # assembly 1 holds one barnase chain (A) and barstar (D)
    return site


@pytest.mark.parametrize("resname, resnum", [("ARG", 83), ("LYS", 27), ("HIS", 102), ("ARG", 87)])
def test_barnase_interface_residues_contact_barstar(resname, resnum):
    site = classify(resname, resnum)
    partners = {(p.label, p.partner_id) for p in site.partners}
    assert ("other_protein", "D") in partners


def test_buried_interface_residue_keeps_its_partner():
    # ARG87 is about 2% accessible on chain A alone but sits against barstar.
    site = classify("ARG", 87)
    assert site.rsa < 0.2
    assert any(p.label == "other_protein" for p in site.partners)


def test_glu73_needs_a_wider_shell_than_4_5_angstrom():
    # Only one non-water atom pair (OE1-ASP39 OD1, 4.49 A) lies within 4.5 A.
    assert not any(p.label == "other_protein" for p in classify("GLU", 73).partners)
    wide = classify("GLU", 73, dataclasses.replace(Config(), contact_distance=5.0))
    assert any(p.label == "other_protein" and p.partner_id == "D" for p in wide.partners)
