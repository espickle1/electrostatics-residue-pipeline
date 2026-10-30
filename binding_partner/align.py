"""Place the reference residue in other chains by local sequence alignment."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from Bio import Align
from Bio.Align import substitution_matrices

from .config import Config

_ALIGNER = Align.PairwiseAligner(
    mode="local",
    substitution_matrix=substitution_matrices.load("BLOSUM62"),
    open_gap_score=-10.0,
    extend_gap_score=-0.5,
)


@lru_cache(maxsize=8192)
def _align(a: str, b: str) -> tuple[tuple, tuple, int, int]:
    """Aligned blocks of a and b, identical columns, aligned columns."""
    if not a or not b:
        return (), (), 0, 0
    best = next(iter(_ALIGNER.align(a, b)))
    blocks_a = tuple((int(s), int(e)) for s, e in best.aligned[0])
    blocks_b = tuple((int(s), int(e)) for s, e in best.aligned[1])
    matches = sum(
        a[i] == b[j]
        for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b)
        for i, j in zip(range(a0, a1), range(b0, b1))
    )
    columns = sum(e - s for s, e in blocks_a)
    return blocks_a, blocks_b, matches, columns


def pair_identity(a: str, b: str) -> float:
    """Identical aligned columns over the shorter sequence (fragments count in full)."""
    _, _, matches, _ = _align(a, b)
    return matches / min(len(a), len(b)) if a and b else 0.0


@dataclass(frozen=True)
class Placement:
    index: int  # index into the target chain's polymer residues
    identity: float  # identical columns / aligned columns against the reference chain


class Placer:
    def __init__(self, ref_seq: str, ref_index: int, cfg: Config):
        self.ref_seq = ref_seq
        self.ref_index = ref_index
        self.cfg = cfg

    def place(self, seq: str) -> Placement | str:
        """A Placement, or the reason the site cannot be placed on this chain."""
        blocks_r, blocks_t, matches, columns = _align(self.ref_seq, seq)
        for (r0, r1), (t0, _) in zip(blocks_r, blocks_t):
            if r0 <= self.ref_index < r1:
                index = t0 + self.ref_index - r0
                break
        else:
            return "position_absent"
        identity = matches / columns if columns else 0.0
        if identity < self.cfg.min_ref_identity:
            return "low_identity"
        if self.cfg.require_same_residue and seq[index] != self.ref_seq[self.ref_index]:
            return "residue_mismatch"
        return Placement(index, identity)
