"""Independent Python nucleotide alignment and V(D)J annotation.

The scoring, gene features, V-then-J window, thresholds and chain filtering are
from MiXCR 3.0.13. V/J/C use the seeded MiLib KMapper and band-limited extension;
D uses exact local linear-gap dynamic programming on NumPy numeric arrays.
Supported runtime profile: human single-end reads, default linear scoring,
all receptor chains.
"""
from copy import deepcopy
import numpy as np
from .kmapper import map_candidates
from .kalignment import extend_hit, shift_homopolymer_indels
from .randomutil import Well19937c, read_seed, replace_wildcards

_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODONS = {a+b+c: _AMINO[i] for i, (a,b,c) in enumerate(
    (a,b,c) for a in _BASES for b in _BASES for c in _BASES)}
_CHAIN_NAMES={"TRA","TRB","TRG","TRD","IGH","IGK","IGL"}
_CHAIN_ALIASES={'ALL':_CHAIN_NAMES,'TCR':{'TRA','TRB','TRG','TRD'},
                'BCR':{'IGH','IGK','IGL'},'IG':{'IGH','IGK','IGL'},
                'TRAD':{'TRA','TRD'},'IGKL':{'IGK','IGL'}}


def chain_names(chains):
    names=set()
    for name in chains.upper().replace('+',',').split(','):
        name=name.strip()
        if name in _CHAIN_ALIASES:names.update(_CHAIN_ALIASES[name])
        elif name in _CHAIN_NAMES:names.add(name)
        else:raise ValueError('Unknown receptor chain: '+name)
    return names


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def translate(sequence):
    # MiLib TranslationParameters.FromCenter preserves both CDR3 anchors.
    length = len(sequence)//3
    left = (length+1)//2
    right = length-left
    first = "".join(_CODONS.get(sequence[i:i+3], "X") for i in range(0, left*3, 3))
    last = "".join(_CODONS.get(sequence[i:i+3], "X")
                   for i in range(len(sequence)-right*3, len(sequence), 3))
    return first + ("_" if len(sequence)%3 else "") + last


def linear_align(reference, query, *, floating_left=True, floating_right=True):
    """Return optimal score, half-open coordinates and reference-boundary map.

    Free ends on both sequences allow reads truncated inside reference genes.
    Fixed right boundary requires reaching either sequence end. A non-floating
    left boundary requires starting at the beginning of one sequence.
    """
    m, n = len(reference), len(query)
    if not m or not n:
        return None
    target = np.frombuffer(query.encode("ascii"), dtype=np.uint8)
    h = np.zeros((m+1, n+1), dtype=np.int32)
    offsets = np.arange(1, n+1, dtype=np.int32) * 12
    for i, base in enumerate(reference.encode("ascii"), 1):
        diag = h[i-1, :-1] + np.where(target == base, 5, -9)
        raw = np.maximum(diag, h[i-1, 1:] - 12)
        if floating_left:
            np.maximum(raw, 0, out=raw)
        # Prefix maximum solves horizontal linear-gap recurrence exactly.
        h[i, 1:] = np.maximum.accumulate(raw + offsets) - offsets
    if floating_right:
        i, j = np.unravel_index(int(h.argmax()), h.shape)
    else:
        lastrow = int(h[m, :].argmax())
        lastcol = int(h[:, n].argmax())
        i, j = (m, lastrow) if h[m, lastrow] >= h[lastcol, n] else (lastcol, n)
    score = int(h[i, j])
    if score <= 0:
        return None
    end_i, end_j = int(i), int(j)
    operations = []
    while i and j and (not floating_left or h[i,j] > 0):
        value = h[i,j]
        if value == h[i-1,j-1] + (5 if reference[i-1] == query[j-1] else -9):
            operations.append("M" if reference[i-1] == query[j-1] else "X")
            i -= 1; j -= 1
        elif value == h[i-1,j] - 12:
            operations.append("D"); i -= 1
        else:
            operations.append("I"); j -= 1
    start_i, start_j = int(i), int(j)
    operations.reverse()
    mapping = {int(i): int(j)}
    for op in operations:
        if op != "D": j += 1
        if op != "I": i += 1
        mapping[int(i)] = int(j)
    return {"score":score, "ref_start":start_i, "ref_end":end_i,
            "query_start":start_j, "query_end":end_j,
            "operations":"".join(operations), "boundary_map":mapping}


def _batch_scores(sequences, query, left, right):
    """Same DP recurrence for many references; traceback is deferred to hits."""
    count,n=len(sequences),len(query)
    if not count:return []
    lengths=np.array([len(s) for s in sequences],dtype=np.int32)
    maxlen=int(lengths.max())
    bases=np.zeros((count,maxlen),dtype=np.uint8)
    for k,seq in enumerate(sequences):bases[k,:len(seq)]=np.frombuffer(seq.encode('ascii'),dtype=np.uint8)
    target=np.frombuffer(query.encode('ascii'),dtype=np.uint8)
    previous=np.zeros((count,n+1),dtype=np.int32)
    best=np.zeros(count,dtype=np.int32)
    offsets=np.arange(1,n+1,dtype=np.int32)*12
    for i in range(1,maxlen+1):
        raw=np.maximum(previous[:,:-1]+np.where(bases[:,i-1,None]==target,5,-9),previous[:,1:]-12)
        if left:np.maximum(raw,0,out=raw)
        current=np.maximum.accumulate(raw+offsets,axis=1)-offsets
        valid=lengths>=i
        if right:
            best=np.maximum(best,np.where(valid,current.max(axis=1),0))
        else:
            best=np.maximum(best,np.where(valid,current[:,-1],0))
            best=np.maximum(best,np.where(lengths==i,current.max(axis=1),0))
        previous[:,1:]=current
    return best.tolist()


def _hits(query, genes, kind, start=0, end=None, preserve_tie_order=False):
    end = len(query) if end is None else end
    segment = query[start:end]
    if not segment:
        return []
    query_kmers = {segment[i:i+5] for i in range(len(segment)-4)}
    absolute = {"V":50, "J":40, "D":25, "C":40}[kind]
    hits = []
    candidates=[]
    for gene in genes:
        # All non-D candidates with at least one exact mapper-size seed are
        # evaluated; no top-candidate truncation based on k-mer frequency.
        if kind != "D" and not gene["kmers"].intersection(query_kmers):
            continue
        candidates.append(gene)
    scores=_batch_scores([g['sequence'] for g in candidates],segment,kind!='C',kind in ('V','D'))
    # Candidates below the absolute threshold cannot influence hit filtering.
    for gene, score in zip(candidates,scores):
        if score<absolute:continue
        result = linear_align(gene["sequence"], segment,
                              floating_left=kind != "C",
                              floating_right=kind in ("V", "D"))
        if result is None or result["score"] < absolute:
            continue
        if kind != "D" and result["query_end"] - result["query_start"] < 15:
            continue
        result["gene"] = gene["name"]
        result["query_start"] += start; result["query_end"] += start
        result["boundary_map"] = {i:j+start for i,j in result["boundary_map"].items()}
        hits.append(result)
    hits.sort(key=(lambda hit: -hit["score"]) if preserve_tie_order
              else (lambda hit: (-hit["score"], hit["gene"])))
    if hits:
        threshold = hits[0]["score"] * (0.85 if kind == "D" else 0.87)
        hits = [h for h in hits if h["score"] >= threshold][:3 if kind == "D" else 7]
    if preserve_tie_order:
        hits.sort(key=lambda hit:(-hit["score"],hit["gene"]))
    return hits


def _map_point(hit, reference, point):
    position = reference.by_name[hit["gene"]]["points"].get(point, -1)
    return hit["boundary_map"].get(position) if position >= 0 else None


def _seeded_hits(query,genes,kind,rng,start=0,end=None,allowed_names=None):
    genes=sorted(genes,key=lambda gene:gene['name'])
    mapped=map_candidates(query,genes,rng,start,end,kind,allowed_names=allowed_names)
    hits=[extend_hit(query,hit,start,end,kind) for hit in mapped]
    hits=[h for h in hits if h is not None]
    hits.sort(key=lambda h:-h['score'])
    if not hits:return []
    threshold=max(50 if kind=='V' else 40,hits[0]['score']*.87)
    if len(hits)>7:threshold=max(threshold,hits[6]['score'])
    return [h for h in hits if h['score']>=threshold]


def _validate_options(options, paired=False):
    options=dict(options or {})
    supported={'allowPartialAlignments','allowNoCDR3PartAlignments','minSumScore','maxHits'}
    if paired:supported.update(('no_merge','smartForceEdgeAlignments'))
    unknown=set(options)-supported
    if unknown:raise ValueError('Unsupported alignment options: '+', '.join(sorted(unknown)))
    for key,value in options.items():
        if key in ('minSumScore','maxHits'):
            if isinstance(value,bool) or not isinstance(value,(int,float)) or value<=0:
                raise ValueError(key+' must be positive')
            if key=='maxHits' and not isinstance(value,int):raise ValueError('maxHits must be an integer')
        elif not isinstance(value,bool):raise ValueError(key+' must be boolean')
    return options


def _align_sequence_partial(sequence, reference, chains, rng, options=None):
    """VDJCAlignerS.processPartial, including J-only rescue and DC RNG order."""
    options=options or {}
    wanted=chain_names(chains)
    genes={t:[g for g in reference.by_type[t] if g['chain_set']&wanted] for t in 'VDJC'}
    results=[]
    for reverse in (False,True):
        seq=reverse_complement(sequence) if reverse else sequence
        v=_seeded_hits(seq,genes['V'],'V',rng)
        vc=set().union(*(reference.by_name[h['gene']]['chain_set'] for h in v)) if v else wanted
        j=_seeded_hits(seq,genes['J'],'J',rng,max(0,v[0]['query_end']-3) if v else 0,
                       allowed_names={g['name'] for g in genes['J'] if g['chain_set']&vc})
        jc=set().union(*(reference.by_name[h['gene']]['chain_set'] for h in j)) if j else set()
        common=vc&jc
        if common and not options.get('allowChimeras',False):
            v=[h for h in v if reference.by_name[h['gene']]['chain_set']&common]
            j=[h for h in j if reference.by_name[h['gene']]['chain_set']&common]
        results.append(dict(reverse_complemented=reverse,oriented_sequence=seq,v_hits=v,j_hits=j))
    for result in results:
        seq=result['oriented_sequence'];v=result['v_hits'];j=result['j_hits'];c=[];d=[]
        if j:
            selected=set().union(*(reference.by_name[h['gene']]['chain_set'] for h in v+j))
            c=_seeded_hits(seq,genes['C'],'C',rng,j[0]['query_end'],
                           allowed_names={g['name'] for g in genes['C'] if g['chain_set']&selected})
        if v and j:
            dc=set().union(*(reference.by_name[h['gene']]['chain_set'] for h in c)) if c else wanted
            d=_hits(seq,[g for g in genes['D'] if g['chain_set']&dc],'D',v[0]['query_end'],j[0]['query_start'])
        total=(v[0]['score'] if v else 0)+(j[0]['score'] if j else 0)
        if options.get('includeDScore',False) and d:total+=d[0]['score']
        if options.get('includeCScore',False) and c:total+=c[0]['score']
        result.update(d_hits=d,c_hits=c,score=total)
    result=max(results,key=lambda r:r['score'])
    v,j=result['v_hits'],result['j_hits'];seq=result['oriented_sequence']
    if not v and not j:return {'status':'NoHits'}
    minimum=options.get('minSumScore',120)
    if result['score']<minimum:return {'status':'LowTotalScore'}
    delta=minimum-result['score'];maxhits=options.get('maxHits',5)
    for kind in 'vj':
        hits=result[kind+'_hits']
        if hits:result[kind+'_hits']=[h for h in hits if h['score']>=delta+hits[0]['score']][:maxhits]
    for kind in 'vdjc':result[kind+'_hits'].sort(key=lambda h:(-h['score'],h['gene']))
    v,j,d=result['v_hits'],result['j_hits'],result['d_hits']
    if not options.get('allowNoCDR3PartAlignments',False):
        available=False
        for hits,point,offset in ((v,'CDR3Begin',-3),(j,'CDR3End',3)):
            if hits:
                p=reference.by_name[hits[0]['gene']]['points'][point]+offset
                available|=hits[0]['ref_start']<=p<=hits[0]['ref_end']
        if not available:return {'status':'NoCDR3Parts'}
    for kind in 'vdjc':
        for hit in result[kind+'_hits']:shift_homopolymer_indels(hit,reference.by_name[hit['gene']]['sequence'],seq)
    left=_map_point(v[0],reference,'CDR3Begin') if v else None
    right=_map_point(j[0],reference,'CDR3End') if j else None
    nt=seq[left:right] if left is not None and right is not None and right>left else ''
    regions=[]
    if nt:
        bounds=[(v[0]['query_end'],d[0]['query_start']),(d[0]['query_end'],j[0]['query_start'])] if d else [(v[0]['query_end'],j[0]['query_start'])]
        regions=[[max(0,a-left),min(len(nt),b-left)] for a,b in bounds if b>a]
    result.update(status='aligned',cdr3_start=left,cdr3_end=right,cdr3_nt=nt,cdr3_aa=translate(nt),n_regions=regions)
    return result


def _align_sequence(sequence, reference, chains, rng, options=None):
    options=options or {}
    if options and options.get('allowPartialAlignments',False):
        return _align_sequence_partial(sequence,reference,chains,rng,options)
    wanted = chain_names(chains)
    by_type = {t:[g for g in reference.by_type[t] if wanted & g["chain_set"]] for t in "VDJC"}
    results = []
    any_v = False
    for reverse in (False, True):
        seq = reverse_complement(sequence) if reverse else sequence
        v = _seeded_hits(seq, by_type["V"], "V",rng)
        any_v |= bool(v)
        if not v: continue
        vchains = set().union(*(reference.by_name[h["gene"]]["chain_set"] for h in v))
        j = _seeded_hits(seq,by_type['J'],'J',rng,max(0,v[0]['query_end']-3),
                         allowed_names={g['name'] for g in by_type['J'] if g['chain_set'] & vchains})
        if not j: continue
        jchains = set().union(*(reference.by_name[h["gene"]]["chain_set"] for h in j))
        common = vchains & jchains
        v = [h for h in v if reference.by_name[h["gene"]]["chain_set"] & common]
        j = [h for h in j if reference.by_name[h["gene"]]["chain_set"] & common]
        if not v or not j: continue
        results.append((v[0]["score"]+j[0]["score"], reverse, seq, v, j, common))
    if not results:
        return {"status":"NoJHits" if any_v else "NoHits"}
    completed=[]
    # When both orientations have V/J, the original evaluates DC in target
    # order. Every C mapper call consumes RNG even when its result loses.
    for total,reverse,seq,v,j,common in results:
        c = _seeded_hits(seq,by_type['C'],'C',rng,j[0]['query_end'],
                         allowed_names={g['name'] for g in by_type['C'] if g['chain_set'] & common})
        # VDJCAlignerS.alignDC runs before calculateHits in 3.0.13. Its vHits and
        # jHits fields are still null, so getPossibleDLoci uses C hits or ALL.
        dchains = set().union(*(reference.by_name[h["gene"]]["chain_set"] for h in c)) if c else wanted
        d = _hits(seq, [g for g in by_type["D"] if g["chain_set"] & dchains], "D",
                  v[0]["query_end"], j[0]["query_start"])
        completed.append((total,reverse,seq,v,j,common,d,c))
    total, reverse, seq, v, j, common,d,c = max(completed, key=lambda item:item[0])
    minimum=options.get('minSumScore',120);maxhits=options.get('maxHits',5)
    if total < minimum:
        return {"status":"LowTotalScore"}
    v = [h for h in v if h["score"] + j[0]["score"] >= minimum][:maxhits]
    j = [h for h in j if h["score"] + v[0]["score"] >= minimum][:maxhits]
    for hits in (v,j,d,c):
        hits.sort(key=lambda h:(-h['score'],h['gene']))
        for hit in hits:shift_homopolymer_indels(hit,reference.by_name[hit['gene']]['sequence'],seq)
    left, right = _map_point(v[0],reference,"CDR3Begin"), _map_point(j[0],reference,"CDR3End")
    cdr3 = seq[left:right] if left is not None and right is not None and right > left else ""
    n_regions = []
    if cdr3:
        vend, jstart = v[0]["query_end"], j[0]["query_start"]
        regions = [(vend,d[0]["query_start"]),(d[0]["query_end"],jstart)] if d else [(vend,jstart)]
        n_regions = [[max(0,a-left),min(len(cdr3),b-left)] for a,b in regions if b > a]
    return {"status":"aligned", "reverse_complemented":reverse,
            "oriented_sequence":seq, "v_hits":v,"j_hits":j,"d_hits":d,"c_hits":c,
            "score":total,"cdr3_start":left,"cdr3_end":right,
            "cdr3_nt":cdr3,"cdr3_aa":translate(cdr3),"n_regions":n_regions}


def align_read(read_id, name, sequence, quality, reference, chains="ALL", options=None):
    options=_validate_options(options)
    sequence = sequence.upper()
    if len(sequence) != len(quality):
        raise ValueError("Sequence and quality lengths differ")
    if not set(sequence) <= set("ACGTRYKMSWBDHVN"):
        raise ValueError("Unsupported nucleotide character")
    raw_sequence,raw_quality=sequence,quality
    sequence,quality=replace_wildcards(read_id,sequence,quality)
    seed=read_seed(read_id,name,sequence)
    key = (sequence,chains,seed,tuple(sorted((options or {}).items())))
    if key not in reference._cache:
        if len(reference._cache) >= 512:
            reference._cache.pop(next(iter(reference._cache)))
        reference._cache[key] = _align_sequence(sequence,reference,chains,Well19937c(seed),options)
    result = deepcopy(reference._cache[key])
    oriented_quality = quality[::-1] if result.get("reverse_complemented") else quality
    result.update(read_id=read_id,name=name,sequence=raw_sequence,quality=raw_quality,
                  oriented_quality=oriented_quality)
    for kind in "vdjc": result.setdefault(kind+"_hits",[])
    result.setdefault("cdr3_nt", ""); result.setdefault("cdr3_aa", "")
    result.setdefault("n_regions",[])
    result["cdr3_quality"] = (oriented_quality[result["cdr3_start"]:result["cdr3_end"]]
                              if result["cdr3_nt"] else "")
    return result


def align_clonal(sequence, reference, chains="TRB", start=0, end=None):
    """D-gene realignment within a clonal CDR3 (same D local score thresholds)."""
    wanted = chain_names(chains)
    genes=[g for g in reference.by_type["D"] if wanted & g["chain_set"]]
    if reference.gene_hashes:
        # CloneAssembler stores the available genes in java.util.HashMap.
        # Capacity derives from the original aligner's selected references.
        # D ties are truncated in this map iteration order.
        capacity=16
        while reference.used_gene_count>capacity*.75:capacity*=2
        def bucket(gene):
            value=reference.gene_hashes[gene["name"]]&0xffffffff
            return (value^(value>>16))&(capacity-1)
        genes.sort(key=bucket)
    return _hits(sequence,genes,"D",start,end,preserve_tie_order=True)
