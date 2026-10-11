"""Mismatch-only paired-read merger ported from MiLib 1.13.

Implements the MiXCR 3.0.13 default Opposite/Unweighted/MaxSubtraction mode.
Candidate order follows the original prefix/suffix motif union search; the
first accepted placement wins. This module does not align unmerged targets.
See the distribution LICENSE for the upstream source terms.
"""
from __future__ import annotations

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _overlap(first: str, qfirst: list[int], second: str, qsecond: list[int], offset: int):
    letters, qualities = [], []
    for position in range(min(0, offset), max(len(first), offset + len(second))):
        i, j = position, position - offset
        has_first, has_second = 0 <= i < len(first), 0 <= j < len(second)
        if not has_first:
            letter, quality = second[j], qsecond[j]
        elif not has_second:
            letter, quality = first[i], qfirst[i]
        elif first[i] == second[j]:
            letter, quality = first[i], max(qfirst[i], qsecond[j])
        else:
            letter = second[j] if qsecond[j] > qfirst[i] else first[i]
            quality = abs(qfirst[i] - qsecond[j])
        letters.append(letter)
        qualities.append(chr(quality + 33))
    return "".join(letters), "".join(qualities)


def merge_pair(seq1: str, qual1: str, seq2: str, qual2: str, *,
               minimal_overlap: int = 17, minimal_identity: float = 0.9) -> dict | None:
    """Return a merged sequence and Phred+33 quality, or None if unmergeable.

    R2 is reverse-complemented. The result retains the R1 strand, including
    when R2 is longer and temporarily becomes the primary search sequence.
    `offset` gives reverse-complemented R2's start relative to R1's start.
    Both sequences must use uppercase or lowercase A/C/G/T. MiLib's merger
    motif builder rejects ambiguous letters; this API reports them explicitly.
    """
    if minimal_overlap < 1 or not 0 <= minimal_identity <= 1:
        raise ValueError("Invalid overlap or identity threshold")
    seq1, seq2 = seq1.upper(), seq2.upper()
    for sequence, quality in ((seq1, qual1), (seq2, qual2)):
        if len(sequence) != len(quality):
            raise ValueError("Sequence and quality lengths differ")
        if set(sequence) - set("ACGT"):
            raise ValueError("Paired-read merger requires unambiguous A/C/G/T sequences")
        if any(not 33 <= ord(char) <= 126 for char in quality):
            raise ValueError("Expected Phred+33 quality characters")
    if min(len(seq1), len(seq2)) < minimal_overlap:
        return None
    first, second = seq1, reverse_complement(seq2)
    qfirst = [ord(char) - 33 for char in qual1]
    qsecond = [ord(char) - 33 for char in qual2[::-1]]
    swapped = len(second) > len(first)
    if swapped:
        first, second = second, first
        qfirst, qsecond = qsecond, qfirst
    motif_length = min(minimal_overlap, 62)
    # Java Math.round on a nonnegative value uses floor(x + 0.5).
    maximum_motif_errors = int(motif_length * (1.0 - minimal_identity) + 0.5)
    motif = [frozenset((second[i], second[len(second) - motif_length + i]))
             for i in range(motif_length)]
    for position in range(len(first) - motif_length + 1):
        if sum(first[position + i] not in motif[i] for i in range(motif_length)) > maximum_motif_errors:
            continue
        overlap = min(len(first) - position, len(second))
        candidates = [(position, position, 0, overlap)]
        end = position + motif_length
        overlap = min(end, len(second))
        candidates.append((min(end - len(second), 0), end - overlap,
                           max(0, len(second) - overlap), overlap))
        for offset, start1, start2, overlap in candidates:
            errors = sum(first[start1 + i] != second[start2 + i] for i in range(overlap))
            identity = (overlap - errors) / overlap
            if identity >= minimal_identity:
                sequence, quality = _overlap(first, qfirst, second, qsecond, offset)
                return {"sequence": sequence, "quality": quality, "overlap": overlap,
                        "identity": identity, "errors": errors,
                        "offset": -offset if swapped else offset,
                        "reversed": True, "score": 2 * overlap - 7 * errors}
    return None
