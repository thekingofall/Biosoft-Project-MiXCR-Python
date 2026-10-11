"""Multi-target PartialAlignmentsAssembler and dedicated V/J/C aligner."""
from collections import Counter,defaultdict
from copy import deepcopy
from .alignment import chain_names,_hits
from .paired import _seeded,_custom_align,_finalize,_map_boundary
from .randomutil import Well19937c,paired_read_seed
from .targetmerger import target_views,merge_ancillary


def _combine(per_target):
    groups={};n=len(per_target)
    for index,hits in enumerate(per_target):
        for hit in hits:
            group=groups.setdefault(hit['gene'],dict(gene=hit['gene'],score=0,alignments=[None]*n))
            group['alignments'][index]=hit;group['score']+=hit['score']
    return sorted(groups.values(),key=lambda h:(-h['score'],h['gene']))


def align_merged_targets(targets,reference,chains='ALL'):
    """PartialAlignmentsAssemblerAligner: chain/order constrained targets."""
    wanted=chain_names(chains)
    genes={k:[g for g in reference.by_type[k] if g['chain_set']&wanted] for k in 'VDJC'}
    rng=Well19937c(paired_read_seed(-1,[('',t['sequence'],t['quality']) for t in targets]))
    n=len(targets);last=0;first_j=-1;last_v=-1;current=wanted
    result=dict(status='aligned',read_id=-1,name='',targets=[dict(sequence=t['sequence'],quality=t['quality']) for t in targets])
    for ki,kind in enumerate('VJC'):
        per=[[] for _ in targets]
        for index in range(last,n):
            start=0
            if ki:
                previous=result['VJC'[ki-1].lower()+'_hits']
                if previous and previous[0]['alignments'][index]:start=previous[0]['alignments'][index]['query_end']
            allowed={g['name'] for g in genes[kind] if g['chain_set']&current} if current else None
            per[index]=_seeded(targets[index]['sequence'],genes[kind],rng,kind,start,allowed)
            if per[index]:
                last=index
                if kind=='V':last_v=index
                if kind=='J' and first_j==-1:first_j=index
        current&=set().union(*(reference.by_name[h['gene']]['chain_set'] for hits in per for h in hits))
        result[kind.lower()+'_hits']=_combine(per)
    v,j=result['v_hits'],result['j_hits']
    if v and j:
        for index in range(1,first_j+1):
            ja=j[0]['alignments'][index]
            if ja is None or ja['query_start']<8:continue
            for hit in v:
                if hit['alignments'][index] is not None or hit['alignments'][index-1] is None:continue
                gene=reference.by_name[hit['gene']];begin=gene['feature_points']['FR3Begin']
                if begin<0:continue
                a=_custom_align(gene['sequence'],targets[index]['sequence'][:ja['query_start']],begin)
                if a['score']>=40:
                    a['gene']=hit['gene'];hit['alignments'][index]=a;hit['score']+=a['score']
    v.sort(key=lambda h:(-h['score'],h['gene']))
    result['v_hits']=v=[h for h in v if h['score']>=v[0]['score']*.87][:7] if v else []
    if v and j:
        for index in range(last_v,n-1):
            va=v[0]['alignments'][index]
            if va is None or len(targets[index]['sequence'])-va['query_end']<8:continue
            for hit in j:
                if hit['alignments'][index] is not None or hit['alignments'][index+1] is None:continue
                a=_custom_align(reference.by_name[hit['gene']]['sequence'],targets[index]['sequence'],query_from=va['query_end'],reverse=True)
                if a['score']>=40:
                    a['gene']=hit['gene'];hit['alignments'][index]=a;hit['score']+=a['score']
    j.sort(key=lambda h:(-h['score'],h['gene']))
    result['j_hits']=j=[h for h in j if h['score']>=j[0]['score']*.87][:7] if j else []
    result['d_hits']=[]
    if v and j:
        vc=set().union(*(reference.by_name[h['gene']]['chain_set'] for h in v))
        jc=set().union(*(reference.by_name[h['gene']]['chain_set'] for h in j));dc=vc&jc or vc|jc
        for index,target in enumerate(targets):
            va,ja=v[0]['alignments'][index],j[0]['alignments'][index]
            if va is None or ja is None:continue
            per=[[] for _ in targets]
            per[index]=_hits(target['sequence'],[g for g in genes['D'] if g['chain_set']&dc],'D',va['query_end'],ja['query_start'])
            result['d_hits']=_combine(per);break
    c=result['c_hits'];result['c_hits']=[h for h in c if h['score']>=c[0]['score']*.87][:5] if c else []
    # The assembler shifts homopolymers after its overlap acceptance gate.
    _finalize(result,reference,normalize=False)
    result['sequence']=result['oriented_sequence'];result['quality']=result['oriented_quality']
    return result


def _point(view,reference,kind,point):
    hit=view['_best'].get(kind)
    if hit is None:return None
    points=reference.by_name[hit['gene']]['points']
    if point=='VEndTrimmed':return hit['query_end'] if points['CDR3Begin']>=0 and hit['ref_end']>=points['CDR3Begin']-3 else None
    if point=='JBeginTrimmed':return hit['query_start'] if points['CDR3End']>=0 and hit['ref_start']<=points['CDR3End']+3 else None
    return _map_boundary(hit,reference.by_name[hit['gene']]['points'][point])


def _left_index(views,reference):
    for index,view in enumerate(views):
        if _point(view,reference,'j','CDR3End') is not None:continue
        if _point(view,reference,'v','VEndTrimmed') is not None:return index
    return -1


def _right_index(views,reference):
    for index,view in enumerate(views):
        if _point(view,reference,'v','CDR3Begin') is not None:continue
        if _point(view,reference,'j','JBeginTrimmed') is not None:return index
    if len(views)!=2:return -1
    def length(view,kind):
        a=view['_best'].get(kind);return a['query_end']-a['query_start'] if a else 0
    def empty(view):return not any(view[k+'_hits'] for k in 'vdjc')
    if length(views[0],'v')==len(views[0]['sequence']) and empty(views[1]):return 1
    if length(views[1],'j')+length(views[1],'c')==len(views[1]['sequence']) and empty(views[0]):return 0
    return -1


def assemble_multiple(records,reference,opts):
    from .partial import _chains
    views=[target_views(r) for r in records]
    lefts=[_left_index(v,reference) for v in views];rights=[_right_index(v,reference) for v in views]
    ids=[r.get('alignment_index',i) for i,r in enumerate(records)]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate alignment_index')
    index=defaultdict(list);indexed=set();not_indexed=set();counters=Counter();k=opts['kValue']
    for i,record in enumerate(records):
        if record.get('cdr3_nt'):continue
        li=lefts[i]
        if li<0:not_indexed.add(ids[i]);continue
        target=views[i][li];start=_point(target,reference,'v','VEndTrimmed')+opts['kOffset'];seq=target['sequence']
        if start<0 or start+k>=len(seq):
            not_indexed.add(ids[i]);counters['leftShortNRegion']+=1;continue
        for position in range(start,len(seq)-k):
            seed=seq[position:position+k]
            if set(seed)-set('ACGT'):counters['wildcardsInKMer']+=1;continue
            index[seed].append((i,position))
        indexed.add(ids[i]);counters['leftParts']+=1
        if counters['leftParts']>=opts['maxLeftParts']:break
    used=set();output=[]
    for ri,right in enumerate(records):
        rid=ids[ri];counters['totalProcessed']+=1
        if right.get('cdr3_nt'):
            counters['containsCDR3']+=1
            if not opts['overlappedOnly']:output.append(deepcopy(right))
            continue
        if rid in used:continue
        rt=rights[ri];accepted=False;blacklist=set()
        if rt>=0:
            counters['rightParts']+=1;sequence=views[ri][rt]['sequence']
            stop=_point(views[ri][rt],reference,'j','JBeginTrimmed')
            stop=len(sequence) if stop is None else stop-opts['kOffset']
            while True:
                best=None
                for position in range(min(stop,len(sequence)-k)):
                    matches=index.get(sequence[position:position+k],[])
                    if len(matches)>opts['maxLeftMatches']:counters['maxRightMatchesLimitReached']=1
                    for mi,(li,lp) in enumerate(matches[:opts['maxLeftMatches']]):
                        left=records[li]
                        if ids[li]==rid or ids[li] in used or li in blacklist:continue
                        if not opts['allowChimeras'] and not (_chains(left,reference,'v')&_chains(right,reference,'j')):continue
                        if left.get('v_hits') and right.get('v_hits') and not ({h['gene'] for h in left['v_hits']}&{h['gene'] for h in right['v_hits']}):continue
                        ls=views[li][lefts[li]]['sequence'];delta=lp-position;begin=max(0,delta);end=min(len(ls),len(sequence)+delta)
                        if ls[begin:end]!=sequence[begin-delta:end-delta]:continue
                        overlap=end-begin
                        if best is None or overlap>best[0]:best=(overlap,li,delta,begin,end,matches,mi)
                if best is None or best[0]<opts['minimalAssembleOverlap']:break
                overlap,li,delta,begin,end,matches,mi=best;removed=matches.pop(mi)
                merged_targets=merge_ancillary(views[li],views[ri],lefts[li],rt,delta,begin,end,reference)
                if merged_targets is None:blacklist.add(li);continue
                if delta<0 or len(views[li][lefts[li]]['sequence'])>=len(sequence)+delta:counters['overOverlapped']+=1
                ot=next((i for i,t in enumerate(merged_targets) if 'overlap_range' in t),None)
                if ot is None:matches.append(removed);break
                merged=align_merged_targets(merged_targets,reference,opts['chains'])
                v,j=merged['v_hits'],merged['j_hits'];a=v[0]['alignments'][ot] if v else None;b=j[0]['alignments'][ot] if j else None
                if a is None or b is None:matches.append(removed);break
                ndn_begin,ndn_end=a['query_end'],b['query_start'];regions=[(ndn_begin,ndn_end)] if ndn_end>ndn_begin else []
                d=merged['d_hits'][0]['alignments'][ot] if merged['d_hits'] else None
                if d:regions=[part for x,y in regions for part in ((x,min(y,d['query_start'])),(max(x,d['query_end']),y)) if part[1]>part[0]]
                begin,end=merged_targets[ot]['overlap_range'];n=sum(y-x for x,y in regions)
                overlap_n=sum(max(0,min(y,end)-max(x,begin)) for x,y in regions);minimum=min(opts['minimalNOverlap'],n)
                if overlap_n<minimum:counters['droppedSmallOverlapNRegion']+=1
                elif minimum==0 and not (begin<=ndn_begin-1<end and begin<=ndn_end<end):counters['droppedNoNRegion']+=1
                else:accepted=True
                if not accepted:matches.append(removed);break
                counters['overlapped']+=1;used.update((ids[li],rid))
                _finalize(merged,reference)
                merged['read_ids']=sorted(set(records[li].get('read_ids',[records[li]['read_id']])+right.get('read_ids',[right['read_id']])))
                merged['read_id']=merged['read_ids'][0];merged['source_alignment_indices']=[ids[li],rid];merged['overlap_range']=[begin,end];merged['overlap_target']=ot
                output.append(merged);break
        if not accepted and opts['writePartial'] and not opts['overlappedOnly'] and rid in not_indexed:
            output.append(deepcopy(right));counters['partialAlignmentsAsIs']+=1
    if opts['writePartial'] and not opts['overlappedOnly']:
        for i,record in enumerate(records):
            if ids[i] in indexed and ids[i] not in used:
                output.append(deepcopy(record));used.add(ids[i]);counters['partialAlignmentsAsIs']+=1
    counters['outputAlignments']=len(output);counters['kMerDiversity']=len(index)
    return dict(alignments=output,report=dict(counters))
