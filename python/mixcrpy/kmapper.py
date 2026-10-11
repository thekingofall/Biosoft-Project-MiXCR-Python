"""MiLib 1.13 KMapper seed voting and candidate filtering, implemented in Python.

Default MiXCR 3.0.13 k=5 parameters. Seed positions are original query
coordinates; hit ``from``/``to`` index the full seed_positions array, while
seed_offsets contains just that hit's slice. Original distribution terms apply.
"""
from functools import lru_cache
import struct

SEED_NOT_FOUND_OFFSET = -2147483647


def _f(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def _max_index(values):
    maximum = max(values)
    return max(i for i, value in enumerate(values) if value == maximum)


def _best_offset(offsets, max_delta=2):
    if len(offsets) == 1:
        return offsets[0][0]
    old, maximum, count = 0x1FFFFFFF, -2147483648, 0
    last_index = len(offsets) - 1
    counters = [0] * (max_delta + 1)
    for i in range(len(offsets)-1, -1, -1):
        current = offsets[i][0]
        if old-current > max_delta:
            if count < last_index-i:
                maximum, count = old-_max_index(counters), last_index-i
            old, last_index = current, i
            counters = [0] * (max_delta+1)
        counters[old-current] += 1
    if count < last_index+1:
        maximum = old-_max_index(counters)
    return maximum


@lru_cache(maxsize=32)
def _base(sequences):
    index = {}
    for gene_id, sequence in enumerate(sequences):
        if set(sequence) - set("ACGT"):
            raise ValueError("KMapper reference must use unambiguous A/C/G/T")
        if len(sequence) < 5:
            raise ValueError("KMapper reference length must be at least five")
        for position in range(len(sequence)-4):
            seed = sequence[position:position+5]
            # MiLib retains position zero but suppresses later homopolymer k-mers.
            if position and len(set(seed)) == 1:
                continue
            index.setdefault(seed, []).append((gene_id, position))
    return index


def map_candidates(sequence, genes, rng, start=0, end=None, kind="V", allowed_names=None,
                   floating_left=None, floating_right=None):
    """Map a query window to genes, consuming the supplied MiLib-compatible RNG.

    The caller supplies references in the same order as the original KAligner.
    D uses a separate local aligner and is outside this mapper's contract.
    """
    floating_left = kind != 'C' if floating_left is None else floating_left
    floating_right = kind == 'V' if floating_right is None else floating_right
    if kind not in {"V", "J", "C"}:
        raise ValueError("KMapper gene kind must be V, J, or C")
    end = len(sequence) if end is None else end
    if not 0 <= start <= end <= len(sequence):
        raise ValueError("Invalid query interval")
    if end-start <= 5:
        return []
    if set(sequence[start:end]) - set("ACGT"):
        raise ValueError("KMapper query requires unambiguous A/C/G/T")
    next_int = getattr(rng, "nextInt", None) or getattr(rng, "next_int")
    seeds, position = [start], start
    while True:
        position += next_int(7)+4
        if position >= end-5:
            break
        seeds.append(position)
    seeds.append(end-5)
    base = _base(tuple(g["sequence"] for g in genes))
    allowed = None if allowed_names is None else set(allowed_names)
    candidates = [[] for _ in genes]
    for seed_id, position in enumerate(seeds):
        for gene_id, ref_position in base.get(sequence[position:position+5], ()):
            if allowed is not None and genes[gene_id]["name"] not in allowed:
                continue
            candidates[gene_id].append((position-ref_position, seed_id))
    hits = []
    for gene_id in range(len(genes)-1, -1, -1):
        offsets = candidates[gene_id]
        if not offsets:
            continue
        offsets.sort()
        offset = _best_offset(offsets)
        score = sum(abs(value-offset) <= 2 for value, _ in offsets)
        left = max(offset, start)
        right = min(offset+len(genes[gene_id]["sequence"]), end)-5
        si_from, si_to = -1, -1
        for j in range(len(seeds)-1, -1, -1):
            if seeds[j] > right:
                continue
            if si_to == -1:
                si_to = j+1
            if seeds[j] >= left:
                continue
            si_from = j+1
            break
        si_from = max(0, si_from)
        if score < 1.5:
            continue
        hits.append({"gene":genes[gene_id], "id":gene_id, "offset":offset,
                     "score":float(score), "from":si_from, "to":si_to,
                     "seed_positions":seeds})
    mismatch, shift = _f(-.1), _f(-.3)
    max_score = 0.
    for hit in reversed(hits):
        offsets = [SEED_NOT_FOUND_OFFSET] * (hit["to"]-hit["from"])
        for offset, seed_id in reversed(candidates[hit["id"]]):
            if not hit["from"] <= seed_id < hit["to"] or abs(offset-hit["offset"]) > 2:
                continue
            i = seed_id-hit["from"]
            if offsets[i] != SEED_NOT_FOUND_OFFSET and abs(offsets[i]-hit["offset"]) <= abs(offset-hit["offset"]):
                continue
            offsets[i] = offset
        prev, c = -1, -1
        for i in range(len(offsets)-1, -1, -1):
            if offsets[i] == SEED_NOT_FOUND_OFFSET:
                continue
            if (c != -1 and offsets[c] != offsets[prev] and offsets[prev] != offsets[i]
                    and (offsets[c] < offsets[prev]) != (offsets[prev] < offsets[i])):
                offsets[prev], prev = SEED_NOT_FOUND_OFFSET, -1
            c, prev = prev, i
        score = 0.
        for offset in offsets:
            score = _f(score + (1. if offset != SEED_NOT_FOUND_OFFSET else mismatch))
        if floating_left:
            prev = -1
            for i in range(len(offsets)):
                if offsets[i] == SEED_NOT_FOUND_OFFSET:
                    continue
                if prev == -1:
                    prev = i
                    continue
                pair_score = _f(_f(1. + _f(abs(offsets[i]-offsets[prev])*shift)) + _f((i-prev-1)*mismatch))
                if pair_score <= 0:
                    offsets[prev], prev = SEED_NOT_FOUND_OFFSET, i
                    continue
                score = _f(score-_f(prev*mismatch))
                break
        if floating_right:
            prev = -1
            for i in range(len(offsets)-1, -1, -1):
                if offsets[i] == SEED_NOT_FOUND_OFFSET:
                    continue
                if prev == -1:
                    prev = i
                    continue
                pair_score = _f(_f(1. + _f(abs(offsets[i]-offsets[prev])*shift)) + _f((prev-1-i)*mismatch))
                if pair_score <= 0:
                    offsets[prev], prev = SEED_NOT_FOUND_OFFSET, i
                    continue
                # The original loop intentionally has no break at this point.
                score = _f(score-_f((len(offsets)-1-prev)*mismatch))
        c, prev = -1, -1
        for i in range(len(offsets)-1, -1, -1):
            if offsets[i] == SEED_NOT_FOUND_OFFSET:
                continue
            ref_position = seeds[i+hit["from"]]-offsets[i]
            if c == -1:
                c, prev = ref_position, i
                continue
            if c <= ref_position:
                offsets[i] = SEED_NOT_FOUND_OFFSET
                continue
            score = _f(score+_f(abs(offsets[i]-offsets[prev])*shift))
            c, prev = ref_position, i
        hit["score"], hit["seed_offsets"] = score, offsets
        max_score = max(max_score, score)
    threshold = max_score * _f(.7 if kind == "V" else .75)
    return [h for h in hits if h["score"] >= 1.5 and h["score"] >= threshold]
