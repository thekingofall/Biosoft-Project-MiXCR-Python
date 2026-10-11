"""CDR3 clone assembly port of MiXCR 3.0.13.

Derived from CloneAssembler, CloneAccumulator, AssemblerUtils,
CloneClusteringStrategy and RelativeConcentrationFilter. Upstream licensing
terms are retained in the distribution's LICENSE.mixcr file.
The runtime is Python; no Java subprocess or reference-output lookup is used.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from math import comb, floor
import re

DEFAULTS = {
    "minimalClonalSequenceLength": 12, "minimalQuality": 0,
    "badQualityThreshold": 20, "maxBadPointsPercent": .7,
    "mappingThreshold": "2of5", "qualityAggregationType": "Max",
    "separateByV": False, "separateByJ": False, "separateByC": False,
    "assemblingFeatures": ["CDR3"],
    "maximalPreClusteringRatio": 1.0,
    "preClusteringScoreFilteringRatio": 2.0,
    "preClusteringCountFilteringRatio": 2.0,
    "addReadsCountOnClustering": False,
    "relativeMinScores": {"v": .87, "j": .87, "c": .87},
    "cloneClusteringParameters": {
        "searchDepth": 2, "allowedMutationsInNRegions": 1,
        "searchParameters": "twoMismatchesOrIndels",
        "clusteringFilter": {"specificMutationProbability": .001},
    },
}
_CODES = {c: i for i, c in enumerate("AGCTNRYSWKMBDHV")}


def sequence_key(sequence):
    """MiLib Sequence.compareTo orders length, then nucleotide codes AGCT."""
    return len(sequence), tuple(_CODES[c] for c in sequence)


def mapping_variants_count(n, k):
    return 1 + sum(comb(n, i) * 3 ** i for i in range(1, k + 1))


def mapping_threshold(n, specification="2of5"):
    """Literal port of MappingThresholdCalculator, including its tail rule."""
    if isinstance(specification, int) or str(specification).isdigit():
        x = int(specification)
        variants = mapping_variants_count(x, x)
    else:
        match = re.fullmatch(r"(\d+)of(\d+)", str(specification))
        if not match:
            raise ValueError("mappingThreshold must be an integer or k-of-n written '2of5'")
        k, size = map(int, match.groups())
        if k > size:
            raise ValueError("mapping threshold k exceeds n")
        variants = mapping_variants_count(size, k)
    for size in range(min(n + 1, 1000)):
        k = 0
        while k < size and mapping_variants_count(size, k + 1) <= variants:
            k += 1
        if k == 1 and size > 1:
            break
        if size == n:
            return k
    return int(n <= variants // 3)


def mutation_path(major, minor, max_errors=2):
    """Banded edit distance with source-coordinate mutation positions.

    Equal-cost paths prefer substitution, deletion, insertion. Repetitive
    indel paths can differ from MiLib's sequence-tree enumeration.
    """
    if abs(len(major) - len(minor)) > max_errors:
        return None
    previous = {j: (j, tuple((0, "I", minor[k]) for k in range(j)))
                for j in range(min(len(minor), max_errors) + 1)}
    for i, a in enumerate(major, 1):
        current = {}
        for j in range(max(0, i-max_errors), min(len(minor), i+max_errors)+1):
            candidates = []
            if j and j-1 in previous:
                cost, path = previous[j-1]
                mismatch = a != minor[j-1]
                candidates.append((cost + mismatch, path + (((i-1, "S", minor[j-1]),) if mismatch else ())))
            if j in previous:
                cost, path = previous[j]
                candidates.append((cost+1, path+((i-1, "D", ""),)))
            if j and j-1 in current:
                cost, path = current[j-1]
                candidates.append((cost+1, path+((i, "I", minor[j-1]),)))
            if candidates:
                selected = min(candidates, key=lambda x: x[0])
                if selected[0] <= max_errors:
                    current[j] = selected
        previous = current
    return previous.get(len(minor), (None, None))[1]


def mutation_positions(major, minor, max_errors=2):
    path = mutation_path(major, minor, max_errors)
    return tuple(p for p, _, _ in path) if path is not None else None


def clonal_compatible(parts, other, mutations=None):
    """Port ClonalSequence.isCompatible, including homopolymers at gaps."""
    if len(parts) != len(other):
        return False
    if mutations is None:
        return all(len(a) == len(b) for a, b in zip(parts, other))
    pointers = []
    length = 0
    for part in parts:
        length += len(part)
        pointers.append(length)

    def count(kind, start, end):
        return sum(t == kind and start <= p < end for p, t, _ in mutations)

    def stretch(index):
        left, right = parts[index], parts[index+1]
        if not left or not right or left[-1] != right[0]:
            return 0, 0
        letter = left[-1]
        left_size = len(left)-len(left.rstrip(letter))
        right_size = len(right)-len(right.lstrip(letter))
        boundary = pointers[index]
        left_delta, right_delta = left_size, right_size
        for p, kind, to in mutations:
            if boundary-left_size <= p < boundary and kind != "D" and to != letter:
                left_delta = boundary-p
        for p, kind, to in reversed(mutations):
            if boundary <= p <= boundary+right_size and kind != "D" and to != letter:
                right_delta = p-boundary+1
        return left_delta, right_delta

    def compatible(index, additional):
        start = pointers[index-1] if index else 0
        end = pointers[index] + (index == len(parts)-1)
        delta = count("I", start, end)-count("D", start, end)
        requested = len(parts[index])+delta+additional-len(other[index])
        if not requested:
            return index == len(parts)-1 or compatible(index+1, 0)
        if index == len(parts)-1:
            return False
        boundary = pointers[index]
        if requested < 0:
            requested = -requested
            if count("I", boundary, boundary+1) >= requested:
                return compatible(index+1, -requested)
            left, right = stretch(index)
            if not left:
                return False
            insertions = count("I", boundary, boundary+right)
            if insertions >= requested or 2*count("D", boundary-left, boundary) >= requested-insertions:
                return compatible(index+1, -requested)
        else:
            left, right = stretch(index)
            if not right:
                return False
            deletions = count("D", boundary, boundary+right)
            if deletions >= requested or 2*count("I", boundary-left, boundary) >= requested-deletions:
                return compatible(index+1, requested)
        return False

    return compatible(0, 0)


def _extract_targets(alignment, features):
    targets = []
    for feature in features:
        target = alignment.get("features", {}).get(feature)
        if target is None and feature == "CDR3" and alignment.get("cdr3_nt"):
            target = {"sequence": alignment["cdr3_nt"], "quality": alignment.get("cdr3_quality"),
                      "n_regions": alignment.get("n_regions"), "contains_junction": True,
                      "contains_d": True, "cdr3_range": [0, len(alignment["cdr3_nt"])]}
        if target is None or target.get("sequence") is None:
            return None
        target = deepcopy(target)
        target["feature"] = feature
        target["sequence"] = target["sequence"].upper()
        if any(c not in _CODES for c in target["sequence"]):
            raise ValueError("Invalid nucleotide in feature " + feature)
        quality = target.get("quality")
        values = [ord(x)-33 for x in quality] if isinstance(quality, str) else quality
        if values is None or len(values) != len(target["sequence"]) or any(type(x) is not int or not 0 <= x <= 93 for x in values):
            raise ValueError("Invalid sequence/quality for feature " + feature)
        target["quality"] = "".join(chr(x+33) for x in values)
        if feature == "CDR3":
            target.setdefault("cdr3_range", [0, len(target["sequence"])])
            target.setdefault("contains_junction", True)
            target.setdefault("contains_d", True)
        targets.append(target)
    return targets


def _signature(alignment, opts):
    return tuple((alignment.get(gt+"_hits") or [{}])[0].get("gene")
                 if opts["separateBy"+gt.upper()] else "*" for gt in "vjc")


def _match_signature(signature, accumulator):
    for gt, gene in zip("vjc", signature):
        if gene == "*":
            continue
        hits = accumulator.get(gt+"_hits", [])
        if gene is None and hits:
            return False
        if gene is not None and hits and not any(h["gene"] == gene for h in hits):
            return False
    return True


def _accumulator_sort_key(acc):
    genes = tuple(acc[gt+"_hits"][0]["gene"] if acc[gt+"_hits"] else "" for gt in "vjc")
    return -_count(acc), genes, sequence_key(acc["sequence"])


def _qualities(alignment, sequence):
    quality = alignment.get("cdr3_quality")
    if quality is None:
        raise ValueError("Every extracted CDR3 requires cdr3_quality")
    values = [ord(x)-33 for x in quality] if isinstance(quality, str) else list(quality)
    if len(values) != len(sequence) or any(x < 0 or x > 93 for x in values):
        raise ValueError("Invalid CDR3 Phred quality or sequence/quality length mismatch")
    return values


def _accumulate(acc, alignment, quality, mapped=False):
    reads = alignment.get("number_of_reads", 1)
    acc["mapped_count" if mapped else "core_count"] += reads
    if mapped:
        return
    acc["quality_values"] = [max(a, b) for a, b in zip(acc["quality_values"], quality)]
    for gt in "vjc":
        for hit in alignment.get(gt+"_hits", []):
            gene = hit["gene"]
            acc["scores"][gt][gene] += hit["score"]
            acc["hit_templates"][gt].setdefault(gene, deepcopy(hit))


def _count(acc):
    return acc["core_count"] + acc["mapped_count"]


def clonal_banded_align(reference, query, right=False, width=5, global_alignment=False):
    """MiLib BandedLinearAligner semi-local, anchored at both left ends.

    Right mode reverses both strings and returns their original coordinates.
    Clipping and band dimensions follow the original wrapper and BandedMatrix.
    """
    original_m, original_n = len(reference), len(query)
    if right:
        reference, query = reference[::-1], query[::-1]
    if not global_alignment:
        limit = min(len(reference), len(query)) + width + 1
        reference, query = reference[:limit], query[:limit]
    m, n = len(reference), len(query)
    width = min(width, m, n)
    row_factor = 2*width + abs(n-m)
    delta = max(0, m-n) + width
    negative = -1073741824
    matrix = {(0, 0): 0}
    for j in range(1, row_factor-delta+1):
        matrix[0, j] = -12*j
    for i in range(1, delta+1):
        matrix[i, 0] = -12*i
    best, stop_i, stop_j = 0, 0, 0
    for i in range(m):
        for j in range(max(0, i-delta), min(i+row_factor-delta+1, n)):
            value = max(matrix.get((i,j), negative) + (5 if reference[i] == query[j] else -9),
                        matrix.get((i,j+1), negative) - 12,
                        matrix.get((i+1,j), negative) - 12)
            matrix[i+1, j+1] = value
            if value > best:
                best, stop_i, stop_j = value, i+1, j+1
    if global_alignment:
        stop_i, stop_j, best = m, n, matrix[m, n]
    i, j = stop_i, stop_j
    operations = []
    while i or j:
        value = matrix[i,j]
        if i and j and value == matrix.get((i-1,j-1), negative) + (5 if reference[i-1] == query[j-1] else -9):
            operations.append("M" if reference[i-1] == query[j-1] else "X")
            i -= 1; j -= 1
        elif i and value == matrix.get((i-1,j), negative)-12:
            operations.append("D"); i -= 1
        else:
            operations.append("I"); j -= 1
    if not right:
        operations.reverse()
    return {"score": best, "ref_start": original_m-stop_i if right else 0,
            "ref_end": original_m if right else stop_i,
            "query_start": original_n-stop_j if right else 0,
            "query_end": original_n if right else stop_j,
            "operations": "".join(operations)}


def annotate_clonal_hits(clone, reference):
    """Re-align V/J CDR3 intersections and then D within trimmed boundaries.

    V/J hit scores retain their assembled mean scores, as in CloneFactory;
    the independent alignment supplies CDR3-relative coordinates. This uses
    the original width-five banded linear algorithm.
    """
    from .alignment import align_clonal
    targets = clone.get("targets") or [{"feature": "CDR3", "sequence": clone["cdr3_nt"], "contains_d": True}]
    d_targets = [i for i, target in enumerate(targets)
                 if target.get("contains_d", target["feature"] == "CDR3")]
    if len(d_targets) > 1:
        raise ValueError("Several assembling features contain the D region")
    for gt in "vjc":
        for hit in clone[gt+"_hits"]:
            gene = reference.by_name[hit["gene"]]
            hit["alignments"] = []
            for target in targets:
                intersection = gene.get("feature_intersections", {}).get(target["feature"])
                if intersection is None and target["feature"] == "CDR3" and gt in "vj":
                    # Preserve custom references created before feature metadata.
                    lo = gene["points"].get("CDR3Begin", -1) if gt == "v" else 0
                    hi = len(gene["sequence"]) if gt == "v" else gene["points"].get("CDR3End", -1)
                    if lo >= 0 and hi >= lo:
                        intersection = {"range": [lo, hi], "mode": "left" if gt == "v" else "right"}
                alignment = None
                if intersection:
                    lo, hi = intersection["range"]
                    part = gene["sequence"][lo:hi]
                    mode = intersection["mode"]
                    alignment = clonal_banded_align(part, target["sequence"], right=mode == "right", global_alignment=mode == "global")
                    alignment["ref_start"] += lo
                    alignment["ref_end"] += lo
                    # A hit's score is the mean read score; the target's
                    # alignment score remains separately available here.
                    if len(targets) == 1:
                        for key in ("query_start", "query_end", "operations", "ref_start", "ref_end"):
                            hit[key] = alignment[key]
                hit["alignments"].append(alignment)
    if not d_targets:
        clone["d_hits"] = []
        return clone
    d_index = d_targets[0]
    sequence = targets[d_index]["sequence"]
    v, j = clone["v_hits"], clone["j_hits"]
    v_alignment = v[0]["alignments"][d_index] if v else None
    j_alignment = j[0]["alignments"][d_index] if j else None
    start = v_alignment["query_end"] if v_alignment else 0
    end = j_alignment["query_start"] if j_alignment else len(sequence)
    all_chains = {"TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"}
    vchains = set().union(*(reference.by_name[h["gene"]]["chain_set"] for h in v)) if v else all_chains
    jchains = set().union(*(reference.by_name[h["gene"]]["chain_set"] for h in j)) if j else all_chains
    chains = vchains & jchains or vchains | jchains
    clone["d_hits"] = (align_clonal(sequence, reference, ",".join(sorted(chains)), start, end)
                       if start < end else [])
    for hit in clone["d_hits"]:
        hit["alignments"] = [deepcopy(hit) if i == d_index else None for i in range(len(targets))]
    return clone


def assemble(alignments, options=None):
    """Assemble alignment dictionaries into clones and read assignments.

    Features are extracted from each record's features dictionary; CDR3
    retains its legacy field fallback. Optional
    clone_annotator(clone) may replace D and clonal alignment coordinates.
    """
    options = options or {}
    unknown = set(options) - (set(DEFAULTS) | {"assemblingFeatures", "clone_annotator"})
    if unknown:
        raise ValueError("Unknown assembly options: " + ", ".join(sorted(unknown)))
    opts = deepcopy(DEFAULTS)
    opts.update(options)
    if opts.get("clone_annotator") is not None and not callable(opts["clone_annotator"]):
        raise ValueError("clone_annotator must be callable")
    for key in ("minimalClonalSequenceLength", "minimalQuality", "badQualityThreshold"):
        if type(opts[key]) is not int or opts[key] < 0:
            raise ValueError(key + " must be a non-negative integer")
    if opts["minimalQuality"] > 93 or opts["badQualityThreshold"] > 93:
        raise ValueError("Quality thresholds must be Phred values from 0 to 93")
    if not isinstance(opts["maxBadPointsPercent"], (int, float)) or not 0 <= opts["maxBadPointsPercent"] <= 1:
        raise ValueError("maxBadPointsPercent must be between 0 and 1")
    for key in ("separateByV", "separateByJ", "separateByC", "addReadsCountOnClustering"):
        if type(opts[key]) is not bool:
            raise ValueError(key + " must be a boolean")
    for key in ("maximalPreClusteringRatio", "preClusteringScoreFilteringRatio", "preClusteringCountFilteringRatio"):
        value = opts[key]
        if not isinstance(value, (int, float)) or not value >= 0 or value == float("inf"):
            raise ValueError(key + " must be finite and non-negative")
        if key != "maximalPreClusteringRatio" and value == 0:
            raise ValueError(key + " must be positive")
    mapping_threshold(0, opts["mappingThreshold"])
    relative = opts["relativeMinScores"]
    if not isinstance(relative, dict) or set(relative) - set("vjc"):
        raise ValueError("relativeMinScores must contain only v, j and c")
    if any(not isinstance(value, (int, float)) or not 0 <= value <= 1 for value in relative.values()):
        raise ValueError("relativeMinScores values must be between 0 and 1")
    opts["relativeMinScores"] = {**DEFAULTS["relativeMinScores"], **relative}
    clustering = opts["cloneClusteringParameters"]
    if clustering is not None:
        allowed = {"searchDepth", "allowedMutationsInNRegions", "searchParameters", "clusteringFilter"}
        if not isinstance(clustering, dict) or set(clustering) - allowed:
            raise ValueError("Unknown or invalid cloneClusteringParameters")
        clustering = {**DEFAULTS["cloneClusteringParameters"], **clustering}
        if clustering["searchParameters"] != "twoMismatchesOrIndels":
            raise NotImplementedError("Only twoMismatchesOrIndels clustering search is implemented")
        for key in ("searchDepth", "allowedMutationsInNRegions"):
            if type(clustering[key]) is not int or clustering[key] < 0:
                raise ValueError(key + " must be a non-negative integer")
        filtering = clustering["clusteringFilter"]
        if not isinstance(filtering, dict) or set(filtering) - {"type", "specificMutationProbability"}:
            raise ValueError("Unknown or invalid clusteringFilter options")
        if filtering.get("type", "relativeConcentration") != "relativeConcentration":
            raise NotImplementedError("Only the relativeConcentration clustering filter is implemented")
        probability = filtering.get("specificMutationProbability", .001)
        if not isinstance(probability, (int, float)) or not 0 <= probability <= 1:
            raise ValueError("specificMutationProbability must be between 0 and 1")
        opts["cloneClusteringParameters"] = clustering
    features = opts["assemblingFeatures"]
    if not isinstance(features, (list, tuple)) or not features or any(not isinstance(f, str) or not f for f in features):
        raise ValueError("assemblingFeatures must be a non-empty list of feature names")
    if len(set(features)) != len(features):
        raise ValueError("assemblingFeatures contains a duplicate feature")
    if opts["qualityAggregationType"] != "Max":
        raise NotImplementedError("Only Max quality aggregation is implemented")
    alignments = list(alignments)
    accumulators = {}
    assignments = []
    deferred = []
    counters = Counter()
    limitations = set()
    seen_read_ids = set()
    for index, alignment in enumerate(alignments):
        read_ids = alignment.get("read_ids", [alignment.get("read_id", index)])
        reads = alignment.get("number_of_reads", len(read_ids))
        if not isinstance(read_ids, (list, tuple)) or not read_ids or any(type(r) is not int for r in read_ids):
            raise ValueError("read_ids must be a non-empty list of integer read identifiers")
        if type(reads) is not int or reads != len(read_ids) or len(set(read_ids)) != len(read_ids):
            raise ValueError("number_of_reads must match the distinct read_ids")
        if seen_read_ids.intersection(read_ids):
            raise ValueError("A read_id occurs in more than one input alignment")
        seen_read_ids.update(read_ids)
        alignment = dict(alignment, number_of_reads=reads)
        assignment = {"read_id": alignment.get("read_id", index),
                      "name": alignment.get("name", str(index)),
                      "read_ids": list(read_ids),
                      "clone_id": None, "status": "dropped"}
        assignments.append(assignment)
        targets = _extract_targets(alignment, features)
        sequence = "".join(t["sequence"] for t in targets) if targets is not None else ""
        if not targets or len(sequence) < opts["minimalClonalSequenceLength"]:
            assignment["reason"] = "failed_to_extract_cdr3" if list(features) == ["CDR3"] else "failed_to_extract_assembling_features"
            counters[assignment["reason"]] += reads
            continue
        parts = tuple(t["sequence"] for t in targets)
        signature = _signature(alignment, opts)
        quality = [ord(c)-33 for t in targets for c in t["quality"]]
        bad = [i for i, q in enumerate(quality) if q <= opts["badQualityThreshold"]]
        if len(bad) > len(sequence) * opts["maxBadPointsPercent"]:
            assignment["reason"] = "too_many_low_quality_points"
            counters[assignment["reason"]] += reads
        elif bad:
            deferred.append((alignment, sequence, parts, signature, quality, bad, assignment))
            counters["deferred"] += reads
        else:
            key = parts, signature
            if key not in accumulators:
                n_regions, offset, skipped = [], 0, 0
                for target in targets:
                    if not target.get("contains_junction", bool(target.get("n_regions"))):
                        skipped += len(target["sequence"])
                        continue
                    if target.get("n_regions") is None:
                        n_regions = None
                        break
                    if skipped:
                        limitations.add("N-region offsets preserve MiXCR 3.0.13's skip of preceding features without junctions.")
                    n_regions.extend([[lo+offset, hi+offset] for lo, hi in target["n_regions"]])
                    offset += len(target["sequence"])
                accumulators[key] = {
                    "sequence": sequence, "parts": parts, "signature": signature, "targets": targets,
                    "cdr3_nt": alignment.get("cdr3_nt", ""), "cdr3_aa": alignment.get("cdr3_aa", ""),
                    "core_count": 0, "mapped_count": 0, "quality_values": [0]*len(sequence),
                    "scores": {gt: Counter() for gt in "vjc"},
                    "hit_templates": {gt: {} for gt in "vjc"},
                    "d_hits": deepcopy(alignment.get("d_hits", [])),
                    "n_regions": n_regions,
                    "assignments": [],
                }
            acc = accumulators[key]
            _accumulate(acc, alignment, quality)
            assignment["status"] = "core"
            acc["assignments"].append(assignment)
    for acc in accumulators.values():
        acc["initial_core_count"] = acc["core_count"]
    for alignment, sequence, parts, signature, quality, bad, assignment in deferred:
        allowed = mapping_threshold(len(bad), opts["mappingThreshold"])
        candidates = []
        for (target_parts, target_signature), acc in accumulators.items():
            if signature != target_signature or not clonal_compatible(parts, target_parts):
                continue
            target = acc["sequence"]
            differences = [i for i, (a,b) in enumerate(zip(sequence,target)) if a != b]
            if len(differences) <= allowed and all(i in bad for i in differences):
                candidates.append((len(differences), sequence_key(target), acc))
        if not candidates:
            assignment["reason"] = "no_candidate_for_deferred_alignment"
            counters[assignment["reason"]] += alignment["number_of_reads"]
            continue
        best_distance = min(c[0] for c in candidates)
        candidates = sorted((c for c in candidates if c[0] == best_distance), key=lambda c: c[1])
        # CloneAssembler 3.0.13 iterates after weight crosses zero, selecting
        # the last candidate. Tree traversal order remains a parity boundary.
        if len(candidates) > 1:
            limitations.add("Ambiguous low-quality mapping uses AGCT sequence order; MiLib tree traversal tie order is unverified.")
        acc = candidates[-1][2]
        _accumulate(acc, alignment, quality, mapped=best_distance > 0)
        assignment["status"] = "mapped" if best_distance else "core"
        assignment["mismatches"] = best_distance
        acc["assignments"].append(assignment)
        counters["mapped_deferred"] += alignment["number_of_reads"]
    containers = {}
    for acc in accumulators.values():
        for gt in "vjc":
            sums = acc["scores"][gt]
            threshold = max(sums.values(), default=0) * opts["relativeMinScores"][gt]
            hits = []
            for gene, score in sums.items():
                if score >= threshold:
                    hit = deepcopy(acc["hit_templates"][gt][gene])
                    hit["score"] = floor(score * 10 / acc["core_count"] + .5) / 10
                    hits.append(hit)
            acc[gt+"_hits"] = sorted(hits, key=lambda h: (-h["score"], h["gene"]))
        containers.setdefault(acc["parts"], []).append(acc)
    survivors = []
    for group in containers.values():
        group.sort(key=_accumulator_sort_key)
        minimum_count = _count(group[0]) / opts["preClusteringCountFilteringRatio"]
        for i in range(len(group)-1):
            major = group[i]
            if major is None:
                continue
            signature = _signature(major, opts)
            threshold = int(_count(major) * opts["maximalPreClusteringRatio"])
            for j in range(i+1, len(group)):
                minor = group[j]
                if minor is not None and _count(minor) <= threshold and _match_signature(signature, minor):
                    major["core_count"] += minor["core_count"]
                    major["mapped_count"] += minor["mapped_count"]
                    for assignment in minor["assignments"]:
                        assignment["preclustered"] = True
                    major["assignments"].extend(minor["assignments"])
                    counters["clones_preclustered"] += 1
                    counters["reads_preclustered"] += _count(minor)
                    group[j] = None
        thresholds = {}
        for gt in "vj":
            eligible = [a[gt+"_hits"][0]["score"] for a in group
                        if a is not None and _count(a) >= minimum_count and a[gt+"_hits"]]
            thresholds[gt] = (max(eligible, default=0) / opts["preClusteringScoreFilteringRatio"]
                              if opts["separateBy"+gt.upper()] else 0)
        # Upstream intentionally iterates to length-1 here; retain its final
        # accumulator exclusion for source-version compatibility.
        for i, acc in enumerate(group):
            if acc is None:
                continue
            score_drop = i < len(group)-1 and any(acc[gt+"_hits"] and acc[gt+"_hits"][0]["score"] < thresholds[gt] for gt in "vj")
            quality_drop = min(acc["quality_values"]) < opts["minimalQuality"]
            if score_drop or quality_drop:
                for assignment in acc["assignments"]:
                    assignment.update(status="droppedwithclone", reason="clone_score_filter" if score_drop else "clone_quality_filter")
                counters["clones_dropped_score" if score_drop else "clones_dropped_quality"] += 1
                continue
            survivors.append(acc)
    survivors.sort(key=lambda a: (_count(a), sequence_key(a["sequence"])), reverse=True)
    clustering = opts["cloneClusteringParameters"]
    assigned = set()
    roots = []
    for acc in survivors:
        if id(acc) in assigned:
            continue
        roots.append(acc)
        assigned.add(id(acc))
        acc["children"] = []
        layer = [acc]
        if not clustering:
            continue
        for _ in range(clustering.get("searchDepth", 2)):
            next_layer = []
            for head in layer:
                for minor in survivors:
                    if id(minor) in assigned or _count(head) <= _count(minor):
                        continue
                    mutations = mutation_path(head["sequence"], minor["sequence"])
                    if mutations is None or not clonal_compatible(head["parts"], minor["parts"], mutations):
                        continue
                    probability = clustering.get("clusteringFilter", {}).get("specificMutationProbability", .001)
                    if _count(minor) > _count(head) * (probability * len(head["sequence"])) ** len(mutations):
                        continue
                    if not _match_signature(_signature(head, opts), minor):
                        continue
                    if head["n_regions"] is None:
                        limitations.add("Clustering candidates with missing N-region coordinates are retained as distinct clones.")
                        continue
                    # Source's 'continue out' counts affected N regions, once
                    # each. Range.containsBoundary includes both endpoints.
                    affected = sum(any(lo <= p <= hi for p, _, _ in mutations) for lo, hi in head["n_regions"])
                    if affected > clustering.get("allowedMutationsInNRegions", 1):
                        continue
                    if any(kind in ("I", "D") for _, kind, _ in mutations):
                        limitations.add("Indel clustering uses Python edit-path tie order; repetitive indel placement parity is unverified.")
                    assigned.add(id(minor))
                    acc["children"].append(minor)
                    next_layer.append(minor)
            layer = sorted(next_layer, key=lambda a: (_count(a), sequence_key(a["sequence"])), reverse=True)
    for root in roots:
        if opts["addReadsCountOnClustering"]:
            for minor in root["children"]:
                root["core_count"] += minor["core_count"]
                root["mapped_count"] += minor["mapped_count"]
    if not clustering or opts["addReadsCountOnClustering"]:
        roots.sort(key=_accumulator_sort_key)
    total_count = sum(_count(root) for root in roots)
    clones = []
    for clone_id, root in enumerate(roots):
        clone = {k: deepcopy(root[k]) for k in ("cdr3_nt", "cdr3_aa", "v_hits", "j_hits", "c_hits", "d_hits")}
        targets = deepcopy(root["targets"])
        pointer = 0
        cdr3_nt, cdr3_quality = "", ""
        for target in targets:
            size = len(target["sequence"])
            target["quality"] = "".join(chr(q+33) for q in root["quality_values"][pointer:pointer+size])
            pointer += size
            if target.get("cdr3_range") is not None:
                lo, hi = target["cdr3_range"]
                cdr3_nt, cdr3_quality = target["sequence"][lo:hi], target["quality"][lo:hi]
        from .alignment import translate
        clone.update(assembling_features=list(features), targets=targets,
                     clonal_sequence=root["sequence"], clonal_quality="".join(chr(q+33) for q in root["quality_values"]),
                     cdr3_nt=cdr3_nt, cdr3_aa=translate(cdr3_nt) if cdr3_nt else "")
        clone.update(clone_id=clone_id, count=_count(root), fraction=_count(root)/total_count,
                     quality=cdr3_quality,
                     core_count=root["core_count"], mapped_count=root["mapped_count"])
        if opts.get("clone_annotator"):
            clone = opts["clone_annotator"](clone)
        elif clone["d_hits"]:
            limitations.add("Clone D hits retain the founding alignment; clonal D re-alignment was not provided.")
        clones.append(clone)
        for assignment in root["assignments"]:
            assignment["clone_id"] = clone_id
        for child in root["children"]:
            for assignment in child["assignments"]:
                assignment.update(clone_id=clone_id, clustered=True,
                                  counted=bool(opts["addReadsCountOnClustering"]))
                counters["clustered_reads"] += len(assignment["read_ids"])
    expanded_assignments = []
    for assignment in assignments:
        read_ids = assignment.pop("read_ids")
        for read_id in read_ids:
            expanded_assignments.append(dict(assignment, read_id=read_id))
    assignments = expanded_assignments
    report = dict(counters)
    report.update(total_alignments=len(alignments), clones=len(clones),
                  total_reads=len(assignments),
                  assigned_reads=sum(a["clone_id"] is not None for a in assignments),
                  clone_count_sum=total_count,
                  dropped_reads=sum(a["clone_id"] is None for a in assignments),
                  limitations=sorted(limitations))
    return {"clones": clones, "read_assignments": assignments, "report": report}
