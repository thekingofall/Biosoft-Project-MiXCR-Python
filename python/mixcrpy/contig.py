"""MiXCR 3.0.13 FullSeqAssembler point variants and contig reconstruction.

Source: assembler/fullseq/FullSeqAssembler.java and MiLib QualityTrimmer.
Upstream licensing terms are retained in LICENSE.mixcr. All runtime code is
Python. Point evidence represents variable-length sequence alleles on the
reference grid, with the assembling feature collapsed into one position.
"""
from collections import defaultdict
from math import ceil, log10

DEFAULTS = dict(branchingMinimalQualityShare=.1, branchingMinimalSumQuality=80,
    decisiveBranchingSumQualityThreshold=120, alignedSequenceEdgeDelta=3,
    alignmentEdgeRegionSize=7, minimalNonEdgePointsFraction=.25,
    outputMinimalQualityShare=.5, outputMinimalSumQuality=50,
    minimalContigLength=20, alignedRegionsOnly=False,
    trimmingParameters=dict(averageQualityThreshold=10, windowSize=8))


def call_point(packed, reads, is_feature=False, options=None):
    """Literal callVariantsForPoint; packed = variant_id<<8 | edge<<7 | Q."""
    p = {**DEFAULTS, **(options or {})}
    groups = defaultdict(list)
    unassigned = set()
    for i in reads:
        if packed[i] == -1: unassigned.add(i)
        else: groups[packed[i] >> 8].append(i)
    if not groups: return [dict(state=-1, reads=set(reads), significant=0)]
    total = sum(packed[i] & 127 for ids in groups.values() for i in ids)
    variants = []; best = -1; best_quality = -1
    for allele, ids in sorted(groups.items()):
        quality = sum(packed[i] & 127 for i in ids)
        nonedge = sum(not packed[i] & 128 for i in ids)
        if nonedge / len(ids) >= p['minimalNonEdgePointsFraction'] and (
            quality >= p['branchingMinimalSumQuality'] and quality >= p['branchingMinimalQualityShare']*total
            or quality >= p['decisiveBranchingSumQualityThreshold']):
            variants.append(dict(state=allele << 8 | min(58, quality), reads=set(ids), significant=len(ids)))
        else:
            unassigned.update(ids)
            if quality > best_quality: best, best_quality = allele, quality
    if variants:
        for v in variants: v['reads'].update(unassigned)
        return variants
    if is_feature or best_quality >= p['outputMinimalQualityShare']*total and best_quality >= p['outputMinimalSumQuality']:
        error = 1-best_quality/total if total else float('nan')
        quality = best_quality if error == 0 else min(0 if total == 0 else int(-10*log10(error)), best_quality)
        return [dict(state=best << 8 | min(58, quality), reads=set(reads), significant=1)]
    return [dict(state=-1, reads=set(reads), significant=0)]


def _cluster(branches, npoints):
    branches.sort(key=lambda b:b['count'])
    observed = [set() for _ in range(npoints)]
    for i in range(len(branches)-1, -1, -1):
        branch = branches[i]; novel = False
        for j,state in enumerate(branch['states']):
            allele = (state & 0xffffffff) >> 8
            if allele not in observed[j]: observed[j].add(allele); novel=True
        if novel: continue
        weights=[]
        for other in branches[i+1:]:
            quality=sum(min(a & 255,b & 255) for a,b in zip(branch['states'],other['states'])
                        if (a & 0xffffffff)>>8 != (b & 0xffffffff)>>8)
            weights.append(10**(-quality/10))
        total=sum(weights)
        if total:
            for other,weight in zip(branches[i+1:],weights): other['count']+=branch['count']*weight/total
        del branches[i]
    branches.sort(key=lambda b:-b['count'])


def _trim(q, left, right, step, rise, threshold, window):
    if not q:return -1 if step==1 else 0
    window=min(right-left,window); limit=ceil(threshold*window)
    end=position=left if step==1 else right-1
    total=0
    for _ in range(window):total+=q[position];position+=step
    if rise == (total >= limit):
        i=0
        while i<window:
            position-=step;end-=step
            if end<0 or end>=len(q):
                while left<=position<right and (rise != (q[position]<threshold)):position-=step
                return position
            total+=q[end]-q[position]
            if rise != (total>=limit):
                while (rise != (q[position]<threshold)) and i<window:position-=step;i+=1
                return position
            i+=1
        return position-step
    while (rise != (total>=limit)) and left<=position<right:
        total+=q[position]-q[end];end+=step;position+=step
    unsuccessful=rise != (total>=limit)
    position-=step
    while left<=position<right and (rise != (q[position]<threshold)):position-=step
    return -2-position if unsuccessful else position


def quality_islands(quality, initial=None, threshold=10, window=8):
    """MiLib calculateAllIslands / calculateIslandsFromInitialRange."""
    q=[ord(c)-33 for c in quality] if isinstance(quality,str) else quality
    absolute=lambda p:p if p>=-1 else -2-p
    ranges=[]
    def islands(start,step):
        while 0<=start<len(q):
            a=_trim(q,start if step==1 else 0,len(q) if step==1 else start,step,True,threshold,window)
            if a < -1:break
            b=absolute(_trim(q,a+step if step==1 else 0,len(q) if step==1 else a+step,step,False,threshold,window))
            if abs(b-a)>=window:
                if step==1:
                    if ranges and ranges[-1][1]==a+1:ranges[-1][1]=b+1
                    else:ranges.append([a+1,b+1])
                elif ranges and ranges[-1][0]==a:ranges[-1][0]=b
                else:ranges.insert(0,[b,a])
            start=b+step
    if initial is None:islands(0,1)
    else:
        a=absolute(_trim(q,0,initial[0],-1,False,threshold,window))
        b=absolute(_trim(q,initial[1],len(q),1,False,threshold,window))+1
        ranges.append([a,b])
        if a>window-1:islands(a-1,-1)
        if b<=len(q)-window:islands(b,1)
    return ranges


def assemble_points(point_reads, cdr3, feature_position, count=None, options=None):
    """Assemble aligned point evidence into quality-trimmed contig branches.

    point_reads is one dict per read: global position -> (sequence, Q, edge).
    Empty sequence encodes a deletion; a multi-base allele encodes an insertion
    or the complete CDR3. Returns branches with target sequence/quality/ranges,
    position maps and the CDR3 target/offset. Clone gene realignment is separate.
    """
    unknown=set(options or {})-set(DEFAULTS)
    if unknown:raise ValueError('Unknown contig options: '+', '.join(sorted(unknown)))
    p={**DEFAULTS,**(options or {})}
    if p['outputMinimalSumQuality']>p['branchingMinimalSumQuality']:raise ValueError('Output quality exceeds branching quality')
    if not point_reads:raise ValueError('Contig assembly requires read evidence')
    variants=list('AGCT')+['',cdr3]; ids={s:i for i,s in enumerate(variants)}
    coverage=defaultdict(int)
    for row in point_reads:
        if feature_position not in row:raise ValueError('Each read must cover the assembling feature')
        for position,(seq,q,edge) in row.items():
            if not 0<=q<=58:raise ValueError('Point quality outside MiLib range 0..58')
            if seq not in ids:ids[seq]=len(variants);variants.append(seq)
            coverage[position]+=1
    points=sorted(coverage,key=lambda i:(-coverage[i],-i))
    branches=[dict(count=float(len(point_reads) if count is None else count),reads=set(range(len(point_reads))),states=[])]
    for point in points:
        packed=[(ids[row[point][0]]<<8 | row[point][1] | (128 if row[point][2] else 0)) if point in row else -1 for row in point_reads]
        next_branches=[]
        for branch in branches:
            choices=call_point(packed,branch['reads'],point==feature_position,p)
            total=sum(v['significant'] for v in choices)
            for v in choices:
                next_branches.append(dict(count=branch['count'] if v['state']==-1 else branch['count']*v['significant']/total,
                    reads=v['reads'],states=branch['states']+[v['state']]))
        branches=next_branches
    _cluster(branches,len(points))
    outputs=[]
    for branch in branches:
        states=sorted((point,state) for point,state in zip(points,branch['states']) if state!=-1)
        if any(point==feature_position and state>>8!=5 for point,state in states):continue
        targets=[];feature_target=-1;feature_offset=-1;block=[]
        for i,(point,state) in enumerate(states):
            seq=variants[state>>8];q=state & 255
            if point==feature_position:feature_target=len(targets);feature_offset=sum(len(x[1]) for x in block)
            block.append((point,seq,q))
            if i+1==len(states) or states[i+1][0]!=point+1:
                mapping=[pos for pos,s,_ in block for _ in s]
                targets.append(dict(sequence=''.join(s for _,s,_ in block),quality=''.join(chr(q+33)*len(s) for _,s,q in block),
                    position_map=mapping,range=[mapping[0],mapping[-1]+1] if mapping else [block[0][0],point+1]))
                block=[]
        if feature_target<0:continue
        clean=[];new_feature=-1;new_offset=-1
        for i,target in enumerate(targets):
            initial=[feature_offset,feature_offset+len(cdr3)] if i==feature_target else None
            trimming=p['trimmingParameters']
            ranges=quality_islands(target['quality'],initial,trimming['averageQualityThreshold'],trimming['windowSize']) if trimming else [[0,len(target['sequence'])]]
            for start,end in ranges:
                is_feature=i==feature_target and start<=feature_offset and end>=feature_offset+len(cdr3)
                if end-start<p['minimalContigLength'] and not is_feature:continue
                mapping=target['position_map'][start:end]
                if is_feature:new_feature=len(clean);new_offset=feature_offset-start
                clean.append(dict(sequence=target['sequence'][start:end],quality=target['quality'][start:end],
                    range=[mapping[0],mapping[-1]+1] if mapping else target['range'],position_map=mapping))
        outputs.append(dict(count=branch['count'],targets=clean,cdr3_target=new_feature,cdr3_offset=new_offset,cdr3_length=len(cdr3)))
    return outputs


def _mutations(hit):
    r=hit['ref_start']; mutations=[]
    operations=hit.get('operations')
    if operations is None:
        raise ValueError('Contig evidence requires alignment operations')
    for op in operations:
        if op!='M':mutations.append((r,op))
        if op!='I':r+=1
    if r!=hit['ref_end']:raise ValueError('Invalid contig alignment operations')
    return mutations


def _query_to_ref(hit,position):
    """Original Mutations.convertToSeq1Position, retaining insertion encoding."""
    position+=hit['ref_start']-hit['query_start']
    seq2p=previous1=previous2=0
    for seq1p,op in _mutations(hit):
        insertion=op=='I'
        if op=='D':seq2p-=1
        if insertion:seq1p-=1;seq2p+=1
        seq2p+=seq1p-previous1
        if seq2p==position and insertion:return -2-(position-previous2+previous1)
        if seq2p>=position:return position-previous2+previous1
        previous1,previous2=seq1p,seq2p
    return position-previous2+previous1


def _ref_to_query(hit,position):
    result=position
    for p,op in _mutations(hit):
        if p>position:break
        if op=='D':
            if p==position:break
            result-=1
        if op=='I':result+=1
    return result+hit['query_start']-hit['ref_start']


def point_evidence(records, clone, reference, options=None):
    """FullSeqAssembler.toPointSequences for CDR3-cluster read alignments.

    Accepts single or multiple oriented targets and per-target hit alignments.
    The clone's leading V/J genes define the shared reference coordinate grid.
    Default subcloning region is absent, as in original default parameters.
    """
    p={**DEFAULTS,**(options or {})}
    if not clone.get('v_hits') or not clone.get('j_hits'):raise ValueError('Contig requires V and J clone hits')
    vgene=clone['v_hits'][0]['gene'];jgene=clone['j_hits'][0]['gene']
    length_v=reference.by_name[vgene]['points']['CDR3Begin']
    joffset=reference.by_name[jgene]['points']['CDR3End']
    if min(length_v,joffset)<0:raise ValueError('Clone V/J genes lack CDR3 boundary')
    feature=1024+length_v;rows=[]
    for record in records:
        cdr3=record.get('cdr3_nt','');cdr3q=record.get('cdr3_quality','')
        if not cdr3 or len(cdr3q)!=len(cdr3):raise ValueError('Contig reads must retain CDR3 sequence and quality')
        row={feature:(cdr3,min(ord(q)-33 for q in cdr3q),True)}
        targets=record.get('targets') or [dict(sequence=record.get('oriented_sequence',record['sequence']),
                                             quality=record.get('oriented_quality',record['quality']))]
        for it,target in enumerate(targets):
            sequence,quality=target['sequence'],target['quality']
            if len(sequence)!=len(quality):raise ValueError('Contig sequence/quality length mismatch')
            for kind,gene,refpoint,offset in [('v',vgene,length_v,1024),('j',jgene,joffset,feature+1-joffset)]:
                parent=next((h for h in record.get(kind+'_hits',[]) if h['gene']==gene),None)
                if parent is None:continue
                hit=parent.get('alignments',[parent])[it]
                if hit is None:continue
                if kind=='v' and hit['ref_start']>length_v:continue
                if kind=='j' and hit['ref_end']<joffset:continue
                _mutations(hit)
                point=_ref_to_query(hit,refpoint) if hit['ref_start']<=refpoint<=hit['ref_end'] else None
                if kind=='v':start=hit['query_start'] if p['alignedRegionsOnly'] else 0;end=point if point is not None else hit['query_end']
                else:start=point if point is not None else hit['query_start'];end=hit['query_end'] if p['alignedRegionsOnly'] else len(sequence)
                aligned_start=max(start,hit['query_start']);aligned_end=min(end,hit['query_end'])
                first=_query_to_ref(hit,aligned_start)
                if aligned_start==aligned_end:last=first
                else:
                    last=_query_to_ref(hit,aligned_end-1)
                    if first<0:first=-2-first
                    if last<0:last=-3-last
                    last+=1
                spans=[]
                shift=offset+hit['ref_start']-hit['query_start']
                spans.extend((i+shift,i,i+1) for i in range(start,aligned_start))
                spans.extend((i+offset,_ref_to_query(hit,i),_ref_to_query(hit,i+1)) for i in range(first,last))
                shift=offset+hit['ref_end']-hit['query_end']
                spans.extend((i+shift,i,i+1) for i in range(aligned_end,end))
                for position,a,b in spans:
                    if position==feature:raise ValueError('Contig flank intrudes into CDR3 position')
                    if a==b:
                        neighbors=[ord(quality[x])-33 for x in (a-1,a+1) if 0<=x<len(quality)]
                        q=sum(neighbors)//len(neighbors) if neighbors else 0
                    else:q=min(ord(x)-33 for x in quality[a:b])
                    # Original default has no subCloningRegion: all points are edge-marked.
                    candidate=(sequence[a:b],q,True)
                    if position not in row or q>row[position][1]:row[position]=candidate
        rows.append(row)
    return rows,feature


def assemble_contig(records, clone, reference, options=None):
    """Reconstruct and annotate one assigned CDR3 clone using original algorithms.

    Includes coverage-based V/J selection, point variants, branch clustering,
    quality trimming and CloneFactory-compatible V/D/J/C target alignments.
    Every supplied alignment must belong to this clone's assignment group.
    """
    records=list(records)
    for kind,point in [('v','CDR3Begin'),('j','CDR3End')]:
        if not any(reference.by_name[h['gene']]['points'][point]>=0 for h in clone.get(kind+'_hits',[])):
            return [_original_clone_branch(clone,'assembly_canceled')]
    base=select_coverage_genes(records,clone,reference)
    rows,feature=point_evidence(records,base,reference,options)
    result=assemble_points(rows,clone['cdr3_nt'],feature,clone.get('count'),options)
    if not result:return [_original_clone_branch(clone,'empty_output_fallback')]
    for branch in result:
        annotate_contig(branch,clone,reference,base,options)
    return result


def _original_clone_branch(clone,status):
    from copy import deepcopy
    branch=dict(count=clone['count'],targets=[dict(sequence=clone['cdr3_nt'],quality=clone.get('quality',clone.get('clonal_quality','')))],
        cdr3_target=0,cdr3_offset=0,cdr3_length=len(clone['cdr3_nt']),reconstruction_status=status)
    for kind in 'vdjc':branch[kind+'_hits']=deepcopy(clone.get(kind+'_hits',[]))
    return branch


def select_coverage_genes(records,clone,reference):
    """CommandAssembleContigs selects V/J by union of covered reference bases."""
    from copy import deepcopy
    base=deepcopy(clone)
    for kind,point in [('v','CDR3Begin'),('j','CDR3End')]:
        candidates=[h for h in clone[kind+'_hits'] if reference.by_name[h['gene']]['points'][point]>=0]
        if not candidates:raise ValueError('Clone lacks a compatible '+kind.upper()+' gene')
        covered={h['gene']:set() for h in candidates}
        for record in records:
            for hit in record.get(kind+'_hits',[]):
                if hit['gene'] not in covered:continue
                for a in hit.get('alignments',[hit]):
                    if a:covered[hit['gene']].update(range(a['ref_start'],a['ref_end']))
        capacity=16
        while len(candidates)>capacity*.75:capacity*=2
        def order(item):
            i,h=item;value=reference.gene_hashes.get(h['gene'],0)&0xffffffff
            return ((value^(value>>16))&(capacity-1),i)
        ordered=[h for _,h in sorted(enumerate(candidates),key=order)]
        best=max(ordered,key=lambda h:len(covered[h['gene']]))
        base[kind+'_hits']=[deepcopy(best)]+[h for h in base[kind+'_hits'] if h['gene']!=best['gene']]
    return base


def _flank_align(reference,query,rstart,rlen,qstart,qlen,side='global',global_reference=False):
    from .kalignment import band_align
    width=2*rlen
    if side=='global':mode='global';added=0;reverse=False
    else:
        mode='added' if global_reference else 'semi';reverse=side=='right';added=width
        if global_reference and qlen>rlen:
            updated=min(qlen,rlen+width)
            if reverse:qstart+=qlen-updated
            qlen=updated;added=min(qlen,width+qlen-rlen)
    ri,qi,ops=band_align(reference[rstart:rstart+rlen],query[qstart:qstart+qlen],width=width,
        mode=mode,added2=min(added,qlen),reverse=reverse)
    return dict(ref_start=rstart+rlen-ri if reverse else rstart,ref_end=rstart+rlen if reverse else rstart+ri,
        query_start=qstart+qlen-qi if reverse else qstart,query_end=qstart+qlen if reverse else qstart+qi,
        operations=ops,score=5*ops.count('M')-9*ops.count('X')-12*(ops.count('D')+ops.count('I')))


def _move_alignment(hit,offset):
    from copy import deepcopy
    hit=deepcopy(hit)
    if hit is None:return None
    hit['query_start']+=offset;hit['query_end']+=offset
    hit.pop('boundary_map',None)
    return hit


def _merge_alignments(left,right):
    if left['ref_end']!=right['ref_start'] or left['query_end']!=right['query_start']:
        raise ValueError('Contig flank and CDR3 alignments are not adjacent')
    return dict(ref_start=left['ref_start'],ref_end=right['ref_end'],query_start=left['query_start'],
        query_end=right['query_end'],score=left['score']+right['score'],operations=left['operations']+right['operations'])


def annotate_contig(branch,clone,reference,base=None,options=None):
    """FullSeqAssembler.buildClone for default linear 5/-9/-12 scoring."""
    from copy import deepcopy
    from .assembly import annotate_clonal_hits
    p={**DEFAULTS,**(options or {})};base=base or clone
    vh=base['v_hits'][0];jh=base['j_hits'][0]
    vr=reference.by_name[vh['gene']];jr=reference.by_name[jh['gene']]
    length_v=vr['points']['CDR3Begin'];jo=jr['points']['CDR3End']
    bound=1024+length_v;rightbound=bound+1;jlen=len(jr['sequence'])-jo
    feature=branch['cdr3_target'];offset=branch['cdr3_offset'];length=branch['cdr3_length']
    val=[];jal=[]
    for i,target in enumerate(branch['targets']):
        lo,hi=target['range'];seq=target['sequence'];v=j=None
        vfloat=not p['alignedRegionsOnly'] and i==0
        if i==feature:
            vfrom=max(0,lo-1024)
            v=_flank_align(vr['sequence'],seq,vfrom,length_v-vfrom,0,offset,
                'right' if lo<1024 or vfloat else 'global',not vfloat)
            qstart=offset+length
            j=_flank_align(jr['sequence'],seq,jo,jlen if hi>=rightbound+jlen else hi-rightbound,
                qstart,len(seq)-qstart,'left' if hi>=rightbound+jlen else 'global',True)
        elif hi<bound:
            if hi>1024:
                vfrom=max(0,lo-1024)
                v=_flank_align(vr['sequence'],seq,vfrom,hi-1024-vfrom,0,len(seq),
                    'right' if lo<1024 or vfloat else 'global',not vfloat)
        elif lo>rightbound:
            if lo<rightbound+jlen:
                jfrom=jo+lo-rightbound
                j=_flank_align(jr['sequence'],seq,jfrom,len(jr['sequence'])-jfrom if hi>=rightbound+jlen else hi-lo,
                    0,len(seq),'left' if hi>=rightbound+jlen else 'global',True)
        else:raise ValueError('Contig target crosses assembling-feature boundary')
        val.append(v);jal.append(j)
    recreated=deepcopy(clone)
    recreated['targets']=[dict(feature='CDR3',sequence=clone['cdr3_nt'],quality=branch['targets'][feature]['quality'][offset:offset+length],contains_d=True)]
    annotate_clonal_hits(recreated,reference)
    # CDR3 alignment for selected base genes is retained from the input clone.
    va=(vh.get('alignments') or [vh])[0];ja=(jh.get('alignments') or [jh])[0]
    val[feature]=_merge_alignments(val[feature],_move_alignment(va,offset))
    jal[feature]=_merge_alignments(_move_alignment(ja,offset),jal[feature])
    for kind in 'vdjc':
        hits=[]
        for hit in recreated[kind+'_hits']:
            moved=deepcopy(hit);old=(hit.get('alignments') or [hit])[0]
            moved['alignments']=[_move_alignment(old,offset) if i==feature else None for i in range(len(branch['targets']))]
            if kind in 'vj' and hit['gene']==(vh if kind=='v' else jh)['gene']:
                moved['alignments']=val if kind=='v' else jal
                if hit['gene']!=recreated[kind+'_hits'][0]['gene']:moved['score']=recreated[kind+'_hits'][0]['score']+1
            for key in ('ref_start','ref_end','query_start','query_end','operations','boundary_map'):moved.pop(key,None)
            hits.append(moved)
        branch[kind+'_hits']=sorted(hits,key=lambda h:(-h['score'],h['gene']))
    return branch
