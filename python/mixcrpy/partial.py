"""MiXCR 3.0.13 partial CDR3 overlap assembly.

Ports PartialAlignmentsAssembler's seeded exact-overlap search, N-region
acceptance and quality sum. Multi-target records use the source-equivalent
ancillary TargetMerger and ordered-target aligner in partial_multi.py.
"""
from collections import Counter, defaultdict
from copy import deepcopy
from .alignment import _hits, _map_point, translate
from .kalignment import extend_hit, shift_homopolymer_indels
from .kmapper import map_candidates
from .randomutil import Well19937c, read_seed

DEFAULTS = dict(kValue=12, kOffset=-7, minimalAssembleOverlap=12,
                minimalNOverlap=7, maxLeftParts=256000, maxLeftMatches=6144,
                writePartial=True, overlappedOnly=False, allowChimeras=False,
                chains="ALL")


def _sequence(record):
    return record.get("oriented_sequence", record["sequence"])


def _quality(record):
    return record.get("oriented_quality", record["quality"])


def _point(record, reference, kind, point):
    hits = record.get(kind + "_hits", [])
    return _map_point(hits[0], reference, point) if hits else None


def _chains(record, reference, kind):
    return set().union(*(reference.by_name[h["gene"]]["chain_set"]
                        for h in record.get(kind + "_hits", [])))


def _sum_overlap(left, right, offset):
    """TargetMerger's collinear SumSubtraction overlap with max quality 45."""
    first, second = _sequence(left), _sequence(right)
    qfirst, qsecond = _quality(left), _quality(right)
    sequence, quality = [], []
    for position in range(min(0, offset), max(len(first), offset + len(second))):
        i, j = position, position-offset
        a, b = 0 <= i < len(first), 0 <= j < len(second)
        if not a:
            base, q = second[j], ord(qsecond[j])-33
        elif not b:
            base, q = first[i], ord(qfirst[i])-33
        elif first[i] == second[j]:
            base, q = first[i], min(45, ord(qfirst[i])+ord(qsecond[j])-66)
        else:
            qa, qb = ord(qfirst[i])-33, ord(qsecond[j])-33
            base, q = (second[j] if qb > qa else first[i]), abs(qa-qb)
        sequence.append(base); quality.append(chr(q+33))
    return "".join(sequence), "".join(quality)


def align_merged_target(sequence, quality, reference, chains="ALL"):
    """Single-target PartialAlignmentsAssemblerAligner, using raw batch scores.

    The merged read has empty description and id -1. J starts at the best V's
    end; this aligner has no VJ overlap-window subtraction or total-score gate.
    """
    from .alignment import chain_names
    wanted = chain_names(chains)
    by_type = {kind: sorted([g for g in reference.by_type[kind] if g["chain_set"] & wanted], key=lambda g:g["name"])
               for kind in "VDJC"}
    rng = Well19937c(read_seed(-1, "", sequence))
    current = wanted
    result = {"read_id": -1, "name": "", "sequence": sequence, "quality": quality,
              "oriented_sequence": sequence, "oriented_quality": quality,
              "reverse_complemented": False, "status": "aligned"}
    previous = None
    for kind in "VJC":
        start = result[previous + "_hits"][0]["query_end"] if previous and result[previous + "_hits"] else 0
        # getFilter0(EMPTY) returns null in the original, as does ALL.
        allowed = {g["name"] for g in by_type[kind] if g["chain_set"] & current} if current else None
        mapped = map_candidates(sequence, by_type[kind], rng, start, kind=kind, allowed_names=allowed)
        hits = [extend_hit(sequence, hit, start, kind=kind) for hit in mapped]
        hits = sorted((h for h in hits if h is not None), key=lambda h:-h["score"])
        if hits:
            threshold = max(40, hits[0]["score"] * .87)
            if len(hits)>7: threshold=max(threshold,hits[6]["score"])
            hits=[h for h in hits if h["score"]>=threshold]
        hits.sort(key=lambda h:(-h["score"],h["gene"]))
        current &= set().union(*(reference.by_name[h["gene"]]["chain_set"] for h in hits))
        result[kind.lower()+"_hits"] = hits
        previous = kind.lower()
    # Final V/J relative cut is maxHits=7; final C is maxHits=5.
    result["v_hits"] = result["v_hits"][:7]
    result["j_hits"] = result["j_hits"][:7]
    result["c_hits"] = result["c_hits"][:5]
    v, j = result["v_hits"], result["j_hits"]
    possible = (_chains(result, reference, "v") & _chains(result, reference, "j"))
    if not possible: possible = _chains(result, reference, "v") | _chains(result, reference, "j")
    result["d_hits"] = (_hits(sequence,[g for g in by_type["D"] if g["chain_set"] & possible],"D",v[0]["query_end"],j[0]["query_start"])
                         if v and j and v[0]["query_end"] <= j[0]["query_start"] else [])
    for kind in "vdjc":
        for hit in result[kind+"_hits"]:
            shift_homopolymer_indels(hit,reference.by_name[hit["gene"]]["sequence"],sequence)
    left, right = _point(result,reference,"v","CDR3Begin"), _point(result,reference,"j","CDR3End")
    cdr3 = sequence[left:right] if left is not None and right is not None and right>left else ""
    result.update(cdr3_start=left,cdr3_end=right,cdr3_nt=cdr3,cdr3_aa=translate(cdr3),
                  cdr3_quality=quality[left:right] if cdr3 else "", n_regions=[])
    if cdr3:
        ranges=[(v[0]["query_end"],j[0]["query_start"])]
        if result["d_hits"]:
            d=result["d_hits"][0]
            ranges=[(v[0]["query_end"],d["query_start"]),(d["query_end"],j[0]["query_start"])]
        result["n_regions"]=[[max(0,a-left),min(len(cdr3),b-left)] for a,b in ranges if b>a]
    return result


def assemble_partial(records, reference, options=None):
    """Overlap complementary partial alignments around CDR3.

    Input dictionaries retain existing alignment fields plus optional
    alignment_index/read_ids. Return alignments and the original report counters.
    """
    options=options or {}
    if set(options)-set(DEFAULTS): raise ValueError("Unknown partial assembly options")
    opts={**DEFAULTS,**options}
    for key in ("kValue","minimalAssembleOverlap","minimalNOverlap","maxLeftParts","maxLeftMatches"):
        if type(opts[key]) is not int or opts[key]<1: raise ValueError("Invalid partial option: "+key)
    if opts["kValue"]>=32: raise ValueError("kValue must be below 32")
    records=list(records)
    if any(len(record.get('targets',[record]))>1 for record in records):
        from .partial_multi import assemble_multiple
        return assemble_multiple(records,reference,opts)
    for record in records:
        if len(_sequence(record))!=len(_quality(record)): raise ValueError("Sequence/quality length mismatch")
    ids=[r.get("alignment_index",i) for i,r in enumerate(records)]
    if len(set(ids))!=len(ids): raise ValueError("Duplicate alignment_index")
    index=defaultdict(list); indexed=set(); not_indexed=set(); counters=Counter()
    k=opts["kValue"]
    for i,record in enumerate(records):
        if record.get("cdr3_nt"): continue
        v=record.get("v_hits",[])
        has_j_end=_point(record,reference,"j","CDR3End") is not None
        start=v[0]["query_end"]+opts["kOffset"] if v and not has_j_end else -1
        sequence=_sequence(record)
        if start<0 or start+k>=len(sequence):
            not_indexed.add(ids[i])
            if v and not has_j_end:counters["leftShortNRegion"]+=1
            continue
        for position in range(start,len(sequence)-k):
            seed=sequence[position:position+k]
            if set(seed)-set("ACGT"):
                counters["wildcardsInKMer"]+=1;continue
            index[seed].append((i,position))
        indexed.add(ids[i]);counters["leftParts"]+=1
        if counters["leftParts"]>=opts["maxLeftParts"]:break
    used=set(); output=[]
    for right_i,right in enumerate(records):
        right_id=ids[right_i]
        counters["totalProcessed"]+=1
        if right.get("cdr3_nt"):
            counters["containsCDR3"]+=1
            if not opts["overlappedOnly"]:output.append(deepcopy(right))
            continue
        if right_id in used:continue
        right_j=right.get("j_hits",[])
        is_right=bool(right_j) and _point(right,reference,"v","CDR3Begin") is None
        best=None
        if is_right:
            counters["rightParts"]+=1
            sequence=_sequence(right)
            stop=right_j[0]["query_start"]-opts["kOffset"]
            for position in range(min(stop,len(sequence)-k)):
                matches=index.get(sequence[position:position+k],[])
                if len(matches)>opts["maxLeftMatches"]:counters["maxRightMatchesLimitReached"]=1
                for match_i,(left_i,left_position) in enumerate(matches[:opts["maxLeftMatches"]]):
                    left=records[left_i];left_id=ids[left_i]
                    if left_id==right_id or left_id in used:continue
                    if not opts["allowChimeras"] and not (_chains(left,reference,"v") & _chains(right,reference,"j")):continue
                    if left.get("v_hits") and right.get("v_hits") and not ({h["gene"] for h in left["v_hits"]}&{h["gene"] for h in right["v_hits"]}):continue
                    left_sequence=_sequence(left);delta=left_position-position
                    begin=max(0,delta);end=min(len(left_sequence),len(sequence)+delta)
                    if left_sequence[begin:end]!=sequence[begin-delta:end-delta]:continue
                    overlap=end-begin
                    if best is None or overlap>best[0]:best=(overlap,left_i,delta,begin,end,matches,match_i)
        accepted=False
        if best and best[0]>=opts["minimalAssembleOverlap"]:
            overlap,left_i,delta,begin,end,matches,match_i=best
            removed=matches.pop(match_i)
            sequence,quality=_sum_overlap(records[left_i],right,delta)
            merged=align_merged_target(sequence,quality,reference,opts["chains"])
            v,j=merged["v_hits"],merged["j_hits"]
            if v and j:
                ndn_begin,ndn_end=v[0]["query_end"],j[0]["query_start"]
                regions=[(ndn_begin,ndn_end)] if ndn_end>ndn_begin else []
                if merged["d_hits"]:
                    d=merged["d_hits"][0]
                    regions=[part for a,b in regions for part in ((a,min(b,d["query_start"])),(max(a,d["query_end"]),b)) if part[1]>part[0]]
                n_length=sum(b-a for a,b in regions)
                overlap_n=sum(max(0,min(b,end)-max(a,begin)) for a,b in regions)
                minimum=min(opts["minimalNOverlap"],n_length)
                if overlap_n<minimum:counters["droppedSmallOverlapNRegion"]+=1
                elif minimum==0 and not (begin<=ndn_begin-1<end and begin<=ndn_end<end):counters["droppedNoNRegion"]+=1
                else:accepted=True
            if accepted:
                counters["overlapped"]+=1
                if delta<0 or len(_sequence(records[left_i]))>=len(_sequence(right))+delta:counters["overOverlapped"]+=1
                merged["read_ids"]=sorted(set(records[left_i].get("read_ids",[records[left_i].get("read_id",ids[left_i])])+right.get("read_ids",[right.get("read_id",right_id)])))
                merged["source_alignment_indices"]=[ids[left_i],right_id]
                merged["overlap_range"]=[begin,end]
                merged["read_id"]=merged["read_ids"][0]
                output.append(merged);used.update((ids[left_i],right_id))
            else:matches.append(removed)
        if not accepted and opts["writePartial"] and not opts["overlappedOnly"] and right_id in not_indexed:
            output.append(deepcopy(right));counters["partialAlignmentsAsIs"]+=1
    if opts["writePartial"] and not opts["overlappedOnly"]:
        # Leftover ordering is stable source order; original uses Trove bucket order.
        for i,record in enumerate(records):
            if ids[i] in indexed and ids[i] not in used:
                output.append(deepcopy(record));used.add(ids[i]);counters["partialAlignmentsAsIs"]+=1
    counters["outputAlignments"]=len(output);counters["kMerDiversity"]=len(index)
    return {"alignments":output,"report":dict(counters)}
