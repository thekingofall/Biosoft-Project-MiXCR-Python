"""MiXCR 3.0.13 VDJCAlignerPVFirst paired-target alignment in Python.

Implements Opposite targets, joint hit scores, edge-forced remapping, FR3/J
rescue, chimera resolution, C/D alignment, and alignment-guided overlap.
"""
from copy import deepcopy
import numpy as np
from .alignment import (align_read, chain_names, reverse_complement, translate,
                        _hits, _align_sequence, _validate_options)
from .kmapper import map_candidates
from .kalignment import extend_hit, shift_homopolymer_indels
from .randomutil import Well19937c, paired_read_seed, replace_wildcards
from .merge import merge_pair, _overlap


def _point(reference,gene,name):
    record=reference.by_name[gene]
    if name in record['points']:return record['points'][name]
    extra=record.get('feature_points',{})
    if name in extra:return extra[name]
    raise ValueError('Reference annotation missing '+gene+':'+name)


def _map_boundary(hit,point):
    return hit['boundary_map'].get(point,hit['boundary_map'].get(str(point)))


def _with_mapping(hit):
    r,q=hit['ref_start'],hit['query_start'];mapping={r:q}
    for op in hit['operations']:
        if op!='I':r+=1
        if op!='D':q+=1
        mapping[r]=q
    hit['boundary_map']=mapping
    return hit


def _custom_align(reference,query,ref_from=0,query_from=0,reverse=False):
    """AlignerCustom semi-local with boundSeq1=false, boundSeq2=true."""
    source=reference[ref_from:];target=query[query_from:]
    if reverse:source=source[::-1];target=target[::-1]
    m,n=len(source),len(target)
    h=np.zeros((m+1,n+1),dtype=np.int32)
    h[0,:]=-12*np.arange(n+1,dtype=np.int32)
    nt=np.frombuffer(target.encode('ascii'),dtype=np.uint8)
    offsets=12*np.arange(1,n+1,dtype=np.int32)
    for i,base in enumerate(source.encode('ascii'),1):
        raw=np.maximum(h[i-1,:-1]+np.where(nt==base,5,-9),h[i-1,1:]-12)
        h[i,1:]=np.maximum.accumulate(raw+offsets)-offsets
    i,j=np.unravel_index(int(h.argmax()),h.shape);i=int(i);j=int(j)
    ei,ej=i,j;score=int(h[i,j]);ops=[]
    while j:
        current=h[i,j]
        if i and current==h[i-1,j-1]+(5 if source[i-1]==target[j-1] else -9):
            ops.append('M' if source[i-1]==target[j-1] else 'X');i-=1;j-=1
        elif i and current==h[i-1,j]-12:ops.append('D');i-=1
        else:ops.append('I');j-=1
    if reverse:
        rs,re=ref_from+m-ei,ref_from+m-i
        qs,qe=query_from+n-ej,query_from+n-j
    else:
        rs,re=ref_from+i,ref_from+ei;qs,qe=query_from+j,query_from+ej;ops.reverse()
    return _with_mapping(dict(score=score,ref_start=rs,ref_end=re,query_start=qs,
                              query_end=qe,operations=''.join(ops)))


def _seeded(query,genes,rng,kind,start=0,allowed=None,floating_left=None,floating_right=None):
    genes=sorted(genes,key=lambda g:g['name'])
    mapped=map_candidates(query,genes,rng,start,len(query),kind,allowed,
                          floating_left,floating_right)
    hits=[extend_hit(query,hit,start,len(query),kind,floating_left,floating_right) for hit in mapped]
    hits=[h for h in hits if h is not None]
    hits.sort(key=lambda h:-h['score'])
    if not hits:return []
    threshold=max(40,hits[0]['score']*.87)
    if len(hits)>7:threshold=max(threshold,hits[6]['score'])
    return [h for h in hits if h['score']>=threshold]


def _chains(hits,reference,empty_all=True):
    if not hits:return chain_names('ALL') if empty_all else set()
    return set().union(*(reference.by_name[h['gene']]['chain_set'] for h in hits))


def _allowed(genes,hits,reference):
    selected=_chains(hits,reference)
    return {g['name'] for g in genes if g['chain_set']&selected}


def _paired_hits(first,second,reference):
    groups={}
    for index,hits in enumerate((first,second)):
        for hit in hits:
            group=groups.setdefault(hit['gene'],dict(gene=hit['gene'],alignments=[None,None]))
            group['alignments'][index]=hit
    # java.util.HashMap iteration precedes a stable score-only sort.
    capacity=16
    while len(groups)>.75*capacity:capacity*=2
    def bucket(hit):
        value=reference.gene_hashes[hit['gene']]&0xffffffff
        return (value^(value>>16))&(capacity-1)
    result=sorted(groups.values(),key=bucket)
    for hit in result:hit['score']=sum(a['score'] for a in hit['alignments'] if a)
    return result


def _sum_filter(hits,kind,maxhits=5):
    hits.sort(key=lambda h:-h['score'])
    if not hits:return []
    minimum=max(50 if kind=='V' else 40,int(hits[0]['score']*.87))
    return [h for h in hits if h['score']>=minimum][:maxhits]


def _chimera(first,second,kind):
    if not first or not second:return False
    if not any(h['score']>120 for h in first) or not any(h['score']>120 for h in second):return False
    if {h['gene'] for h in first}&{h['gene'] for h in second}:return False
    if kind=='V':
        (second if first[0]['ref_end']>second[0]['ref_end'] else first).clear()
    else:
        (first if first[0]['ref_start']>second[0]['ref_start'] else second).clear()
    return True


def _initial(targets,genes,rng,reference,options):
    v=[_seeded(t['sequence'],genes['V'],rng,'V') for t in targets]
    if options.get('smartForceEdgeAlignments',True):
        force_right=force_left=False
        for paired in _sum_filter(_paired_hits(*v,reference),'V',options.get('maxHits',5)):
            left,right=paired['alignments']
            if left is None or right is None:continue
            if left['query_end']!=len(targets[0]['sequence']) and right['query_start']==0:force_right=True
            elif left['query_end']==len(targets[0]['sequence']) and right['query_start']!=0:force_left=True
        if force_right:v[0]=_seeded(targets[0]['sequence'],genes['V'],rng,'V',floating_right=False)
        if force_left:v[1]=_seeded(targets[1]['sequence'],genes['V'],rng,'V',floating_left=False)
    vchimera=_chimera(*v,'V')
    vhits=_paired_hits(*v,reference)
    for hit in vhits:
        left,right=hit['alignments']
        if right is not None or left['query_end']<len(targets[0]['sequence'])-5:continue
        begin=_point(reference,hit['gene'],'FR3Begin')
        if begin<0:continue
        rescue=_custom_align(reference.by_name[hit['gene']]['sequence'],targets[1]['sequence'],begin)
        if rescue['score']<40:continue
        rescue['gene']=hit['gene'];hit['alignments'][1]=rescue;hit['score']+=rescue['score']
    vhits=_sum_filter(vhits,'V',options.get('maxHits',5))
    j=[]
    for index,target in enumerate(targets):
        vhit=vhits[0]['alignments'][index] if vhits else None
        start=0
        if vhit:
            if vhit['ref_end']<=_point(reference,vhit['gene'],'FR3Begin') or vhit['query_end']==len(target['sequence']):
                j.append([]);continue
            start=max(0,vhit['query_end']-3)
        j.append(_seeded(target['sequence'],genes['J'],rng,'J',start,_allowed(genes['J'],vhits,reference)))
    jchimera=_chimera(*j,'J')
    jhits=_paired_hits(*j,reference)
    for hit in jhits:
        left,right=hit['alignments']
        if left is not None or right['query_start']>5:continue
        start=max(0,vhits[0]['alignments'][0]['query_end']-3) if vhits and vhits[0]['alignments'][0] else 0
        rescue=_custom_align(reference.by_name[hit['gene']]['sequence'],targets[0]['sequence'],query_from=start,reverse=True)
        if rescue['score']<40:continue
        rescue['gene']=hit['gene'];hit['alignments'][0]=rescue;hit['score']+=rescue['score']
    jhits=_sum_filter(jhits,'J',options.get('maxHits',5))
    if not options.get('allowChimeras',False):
        common=_chains(vhits,reference)&_chains(jhits,reference)
        if common:
            vhits=[h for h in vhits if reference.by_name[h['gene']]['chain_set']&common]
            jhits=[h for h in jhits if reference.by_name[h['gene']]['chain_set']&common]
    return dict(targets=targets,v_hits=vhits,j_hits=jhits,
                v_chimera=vchimera,j_chimera=jchimera,
                score=(vhits[0]['score'] if vhits else 0)+(jhits[0]['score'] if jhits else 0))


def _combine(first,second,reference):
    hits=_paired_hits(first,second,reference)
    hits.sort(key=lambda h:(-h['score'],h['gene']))
    return hits


def _dc(helper,genes,rng,reference):
    v,j=helper['v_hits'],helper['j_hits'];targets=helper['targets']
    c=[[],[]]
    if j:
        allowed=_allowed(genes['C'],v+j,reference)
        for index,target in enumerate(targets):
            ja=j[0]['alignments'][index]
            if index==0 and ja is None:continue
            c[index]=_seeded(target['sequence'],genes['C'],rng,'C',ja['query_end'] if ja else 0,allowed)
    elif v and v[0]['alignments'][0] and v[0]['alignments'][1] is None:
        c[1]=_seeded(targets[1]['sequence'],genes['C'],rng,'C',0,_allowed(genes['C'],v,reference))
    helper['c_hits']=_combine(*c,reference)
    d=[[],[]]
    if v and j:
        groups=(v,j,helper['c_hits'])
        chains=set.intersection(*(_chains(group,reference) for group in groups))
        if not chains:chains=set.union(*(_chains(group,reference,False) for group in groups))
        dgenes=[g for g in genes['D'] if g['chain_set']&chains]
        for index,target in enumerate(targets):
            va,ja=v[0]['alignments'][index],j[0]['alignments'][index]
            if va is None or ja is None or va['query_end']>=ja['query_start']:continue
            d[index]=_hits(target['sequence'],dgenes,'D',va['query_end'],ja['query_start'])
    helper['d_hits']=_combine(*d,reference)


def _raw_position(hit,position):
    if position<hit['ref_start'] or position>hit['ref_end']:return -1
    r,q=hit['ref_start'],hit['query_start']
    for op in hit['operations']:
        if r==position and op!='I':return -2-q if op=='D' else q
        if op!='I':r+=1
        if op!='D':q+=1
    return q if r==position else -1


def _alignment_merge(result,reference,rng,chains,options=None):
    targets=result['targets']
    for kind in 'vjc':
        hits=result[kind+'_hits']
        if not hits:continue
        if not all(any(h['alignments'][i] is not None for h in hits) for i in (0,1)):continue
        left,right=hits[0]['alignments']
        if left is None or right is None:continue
        start=max(left['ref_start'],right['ref_start']);end=min(left['ref_end'],right['ref_end'])
        if end<=start:continue
        delta=_raw_position(left,start)-_raw_position(right,start)
        if delta!=_raw_position(left,end)-_raw_position(right,end):continue
        offset0=max(delta,0);offset1=max(-delta,0)
        overlap=min(len(targets[0]['sequence'])-offset0,len(targets[1]['sequence'])-offset1)
        if overlap<=0:continue
        mismatches=sum(targets[0]['sequence'][offset0+i]!=targets[1]['sequence'][offset1+i] for i in range(overlap))
        if (overlap-mismatches)/overlap<float(np.float32(.9)):
            remove=1 if left['score']>right['score'] else 0
            hits[0]['alignments'][remove]=None
            hits[0]['score']=sum(a['score'] for a in hits[0]['alignments'] if a)
            # VDJCObject construction resorts the modified hit list.
            hits.sort(key=lambda h:(-h['score'],h['gene']))
            result['alignment_merge_conflict']=kind
            return result
        first,second=targets
        if delta<0:first,second=second,first;delta=-delta
        sequence,quality=_overlap(first['sequence'],[ord(x)-33 for x in first['quality']],
                                  second['sequence'],[ord(x)-33 for x in second['quality']],delta)
        merged=_align_sequence(sequence,reference,chains,rng,options)
        if merged['status']!='aligned':return result
        return _single_result(merged,sequence,quality,'alignment',dict(overlap=overlap,errors=mismatches,offset=delta))
    return result


def _single_result(result,sequence,quality,method,merge_info=None):
    result=deepcopy(result)
    oriented_quality=quality[::-1] if result.get('reverse_complemented') else quality
    result['oriented_quality']=oriented_quality
    result['targets']=[dict(sequence=result.get('oriented_sequence',sequence),quality=oriented_quality)]
    result['number_of_reads']=1
    result['paired_ends']=2
    result['merge_method']=method
    if merge_info:result['merge_info']=merge_info
    for kind in 'vdjc':
        for hit in result.get(kind+'_hits',[]):
            hit['alignments']=[deepcopy(hit)]
    result['cdr3_quality']=oriented_quality[result['cdr3_start']:result['cdr3_end']] if result.get('cdr3_nt') else ''
    return result


def _finalize(result,reference,normalize=True):
    targets=result['targets']
    for kind in 'vdjc':
        hits=result.setdefault(kind+'_hits',[])
        hits.sort(key=lambda h:(-h['score'],h['gene']))
        for hit in hits:
            for index,alignment in enumerate(hit['alignments']):
                if alignment and normalize:
                    shift_homopolymer_indels(alignment,reference.by_name[hit['gene']]['sequence'],targets[index]['sequence'])
    candidates=[]
    for index,target in enumerate(targets):
        if not result['v_hits'] or not result['j_hits']:continue
        vh,jh=result['v_hits'][0],result['j_hits'][0]
        va,ja=vh['alignments'][index],jh['alignments'][index]
        if va is None or ja is None:continue
        left=_map_boundary(va,_point(reference,vh['gene'],'CDR3Begin'))
        right=_map_boundary(ja,_point(reference,jh['gene'],'CDR3End'))
        if left is None or right is None or right<=left:continue
        quality=target['quality'][left:right]
        candidates.append((min(map(ord,quality)),index,left,right))
    result.update(cdr3_nt='',cdr3_aa='',cdr3_quality='',cdr3_start=None,cdr3_end=None,n_regions=[])
    if candidates:
        _,index,left,right=max(candidates,key=lambda c:c[0])
        target=targets[index];nt=target['sequence'][left:right]
        result.update(cdr3_nt=nt,cdr3_aa=translate(nt),cdr3_quality=target['quality'][left:right],
                      cdr3_start=left,cdr3_end=right,cdr3_target=index)
        va=result['v_hits'][0]['alignments'][index];ja=result['j_hits'][0]['alignments'][index]
        dh=result['d_hits'][0]['alignments'][index] if result['d_hits'] else None
        regions=[(va['query_end'],dh['query_start']),(dh['query_end'],ja['query_start'])] if dh else [(va['query_end'],ja['query_start'])]
        result['n_regions']=[[max(0,a-left),min(len(nt),b-left)] for a,b in regions if b>a]
    else:index=0
    result['oriented_sequence']=targets[index]['sequence'];result['oriented_quality']=targets[index]['quality']
    for kind in 'vdjc':
        for hit in result[kind+'_hits']:
            alignment=hit['alignments'][index]
            if alignment:
                for key,value in alignment.items():
                    if key not in ('gene','score'):hit[key]=deepcopy(value)
    return result


def align_pair(read_id,name1,seq1,qual1,name2,seq2,qual2,reference,chains='ALL',options=None):
    options=_validate_options(options,paired=True)
    for sequence,quality in ((seq1,qual1),(seq2,qual2)):
        if len(sequence)!=len(quality):raise ValueError('Sequence and quality lengths differ')
        if not set(sequence.upper())<=set('ACGTRYKMSWBDHVN'):raise ValueError('Unsupported nucleotide character')
    single_options={k:v for k,v in options.items() if k not in ('no_merge','smartForceEdgeAlignments')}
    seq1,qual1=replace_wildcards(read_id,seq1.upper(),qual1)
    seq2,qual2=replace_wildcards(read_id,seq2.upper(),qual2)
    if not options.get('no_merge',False):
        merged=merge_pair(seq1,qual1,seq2,qual2)
        if merged is not None:
            single=align_read(read_id,'',merged['sequence'],merged['quality'],reference,chains,single_options)
            result=_single_result(single,merged['sequence'],merged['quality'],'sequence',merged)
            result.update(read_id=read_id,name=name1,sequence=seq1,quality=qual1)
            return result
    rng=Well19937c(paired_read_seed(read_id,[(name1,seq1,qual1),(name2,seq2,qual2)]))
    wanted=chain_names(chains)
    genes={kind:[g for g in reference.by_type[kind] if g['chain_set']&wanted] for kind in 'VDJC'}
    target_sets=[
        [dict(sequence=seq1,quality=qual1),dict(sequence=reverse_complement(seq2),quality=qual2[::-1])],
        [dict(sequence=seq2,quality=qual2),dict(sequence=reverse_complement(seq1),quality=qual1[::-1])]]
    helpers=[_initial(targets,genes,rng,reference,options) for targets in target_sets]
    partial=options.get('allowPartialAlignments',False)
    if partial:
        for helper in helpers:
            _dc(helper,genes,rng,reference)
            for kind in 'dc':
                if options.get('include'+kind.upper()+'Score',False) and helper[kind+'_hits']:
                    helper['score']+=helper[kind+'_hits'][0]['score']
    result=max(helpers,key=lambda h:h['score'])
    v,j=result['v_hits'],result['j_hits']
    minimum=options.get('minSumScore',120)
    if not v and not j or not partial and (not v or not j):result['status']='NoHits' if not v and not j else 'NoVHits' if not v else 'NoJHits'
    else:
        if not partial:_dc(result,genes,rng,reference)
        if result['score']<minimum:result['status']='LowTotalScore'
        else:
            difference=minimum-result['score'];maxhits=options.get('maxHits',5)
            if v:
                vmin=int(max(v[0]['score']*.87,difference+v[0]['score']))
                result['v_hits']=[h for h in v if h['score']>=vmin][:maxhits]
            if j:
                jmin=int(difference+j[0]['score'])
                result['j_hits']=[h for h in j if h['score']>=jmin][:maxhits]
            result['status']='aligned' if partial or any(v[0]['alignments'][i] and j[0]['alignments'][i] for i in (0,1)) else 'VAndJOnDifferentTargets'
            if partial and not options.get('allowNoCDR3PartAlignments',False):
                # Paired mode tests coverage towards the junction; single mode
                # uses availability of the outer CDR3-adjacent points.
                available=False
                for kind,point,offset in (('v','CDR3Begin',3),('j','CDR3End',-3)):
                    hits=sorted(result[kind+'_hits'],key=lambda h:(-h['score'],h['gene']))
                    if not hits:continue
                    bound=_point(reference,hits[0]['gene'],point)+offset
                    for hit in hits[0]['alignments']:
                        if hit:available|=(hit['ref_end']>=bound if kind=='v' else hit['ref_start']<=bound)
                if not available:result['status']='NoCDR3Parts'
    if result['status']=='aligned':
        for kind in 'vdjc':result[kind+'_hits'].sort(key=lambda h:(-h['score'],h['gene']))
        if not options.get('no_merge',False):result=_alignment_merge(result,reference,rng,chains,options)
        if len(result['targets'])==2:result=_finalize(result,reference)
    else:
        for kind in 'vdjc':result[kind+'_hits']=[]
        result.update(cdr3_nt='',cdr3_aa='',cdr3_quality='',n_regions=[])
    result.update(read_id=read_id,name=name1,sequence=seq1,quality=qual1,number_of_reads=1,paired_ends=2)
    return result
