import pytest

from binding_partner.analyze import analyze_file, make_reference
from binding_partner.config import Config


@pytest.fixture
def analyze(tmp_path):
    """analyze(builder, cfg=None, name='x.pdb') -> FileResult for A:GLU2 of the built file."""

    def run(builder, cfg=None, name="x.pdb", ref_chain="A", ref_resname="GLU", ref_resnum=2):
        path = builder.write(tmp_path / name)
        ref = make_reference(path, ref_chain, ref_resname, ref_resnum)
        return analyze_file(path, ref, cfg or Config())

    return run


def partner_keys(result):
    """{(label, component)} over every placed copy."""
    return {(p.label, p.component) for s in result.sites for p in s.partners}
