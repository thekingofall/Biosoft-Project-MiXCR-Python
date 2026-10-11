"""AlignedTarget and TargetMerger operations for partial-read assembly."""
from copy import deepcopy
from .paired import _raw_position,_with_mapping,_paired_hits
from .kalignment import band_align
from .merge import merge_pair
from .alignment import reverse_complement


def _trove_order(keys,hashes):
    """TObjectLongHashMap insertion/probing and descending bucket traversal."""
    # PrimeFinder capacity values from the original bundled Trove library.
    primes=(23,31,37,43,47,67,79,89,97,137,163,179,197,277,311,331,359,379,397,433,557,599,631,673,719,761,797,877,953,1039,1117,1201,1277,1361,1439,1523,1597,1759,1907,2081,2237,2411,2557,2729,2879,3049,3203,3527,3821,4177,4481,4831,5119,5471,5779,6101,6421,7057,7643,8363,8963,9677,10243)
    table=[None]*23;size=0
    def insert(table,key):
        value=hashes[key]&0x7fffffff;pos=value%len(table);step=1+value%(len(table)-2)
        while table[pos] is not None:
            if table[pos]==key:return False
            pos=(pos-step)%len(table)
        table[pos]=key;return True
    for key in keys:
        if not insert(table,key):continue
        size+=1
        if size>len(table)//2:
            capacity=next((p for p in primes if p>=2*len(table)),None)
            if capacity is None:raise ValueError('Reference exceeds supported Trove capacity table')
            old=table;table=[None]*capacity
            for entry in reversed(old):
                if entry is not None:insert(table,entry)
    return [key for key in reversed(table) if key is not None]


def target_views(record):
    targets=record.get('targets',[dict(sequence=record.get('oriented_sequence',record['sequence']),quality=record.get('oriented_quality',record['quality']))])
    result=[]
    for index,target in enumerate(targets):
        view=dict(sequence=target['sequence'],quality=target['quality'],_origin=id(record),_index=index,_best={},_parent=record)
        for kind in 'vdjc':
            hits=record.get(kind+'_hits',[]);view[kind+'_hits']=[]
            for pos,hit in enumerate(hits):
                alignment=hit.get('alignments',[hit])[index]
                if pos==0:view['_best'][kind]=alignment
                if alignment:
                    a=deepcopy(alignment);a['gene']=hit['gene'];view[kind+'_hits'].append(a)
        result.append(view)
    return result


def forced_merge(left,right,offset,reference):
    from .partial import _sum_overlap
    if offset<0:return forced_merge(right,left,-offset,reference)
    sequence,quality=_sum_overlap(left,right,offset)
    result=dict(sequence=sequence,quality=quality,_origin=None)
    for kind in 'vjc':
        result[kind+'_hits']=[]
        names={h['gene'] for h in left[kind+'_hits']+right[kind+'_hits']}
        for name in names:
            a=next((h for h in left[kind+'_hits'] if h['gene']==name),None)
            b=next((h for h in right[kind+'_hits'] if h['gene']==name),None)
            if a and b and _raw_position(a,b['ref_start'])!=b['query_start']+offset:
                if a['score']>b['score']:b=None
                else:a=None
            if a and b:
                if a['ref_start']<b['ref_start']:rs,qs=a['ref_start'],a['query_start']
                else:rs,qs=b['ref_start'],b['query_start']+offset
                if a['ref_end']>b['ref_end']:re,qe=a['ref_end'],a['query_end']
                else:re,qe=b['ref_end'],b['query_end']+offset
            elif a:rs,re,qs,qe=a['ref_start'],a['ref_end'],a['query_start'],a['query_end']
            else:rs,re,qs,qe=b['ref_start'],b['ref_end'],b['query_start']+offset,b['query_end']+offset
            _,_,ops=band_align(reference.by_name[name]['sequence'][rs:re],sequence[qs:qe])
            score=5*ops.count('M')-9*ops.count('X')-12*(ops.count('I')+ops.count('D'))
            result[kind+'_hits'].append(_with_mapping(dict(gene=name,score=score,ref_start=rs,ref_end=re,query_start=qs,query_end=qe,operations=ops)))
        result[kind+'_hits'].sort(key=lambda h:(-h['score'],h['gene']))
    result['d_hits']=[]
    if 'overlap_range' in left:result['overlap_range']=left['overlap_range'][:]
    if 'overlap_range' in right:result['overlap_range']=[p+offset for p in right['overlap_range']]
    return result


def merge_targets(left,right,reference):
    """Return (merged, score, conflict); sequence fallback uses 20 nt / 95%."""
    for kind in 'vjc':
        if not left[kind+'_hits'] or not right[kind+'_hits']:continue
        if left.get('_origin') is not None and left.get('_origin')==right.get('_origin'):
            hits=left['_parent'][kind+'_hits']
            top=hits[0];a=top['alignments'][left['_index']];b=top['alignments'][right['_index']]
        else:
            pairs=_paired_hits(left[kind+'_hits'],right[kind+'_hits'],reference)
            pairs.sort(key=lambda h:-h['score']);a,b=pairs[0]['alignments']
        if not a or not b:continue
        start,end=max(a['ref_start'],b['ref_start']),min(a['ref_end'],b['ref_end'])
        if end<=start:continue
        delta=_raw_position(a,start)-_raw_position(b,start)
        if delta!=_raw_position(a,end)-_raw_position(b,end):continue
        x,y=max(delta,0),max(-delta,0)
        overlap=min(len(left['sequence'])-x,len(right['sequence'])-y)
        errors=sum(left['sequence'][x+i]!=right['sequence'][y+i] for i in range(overlap))
        if (overlap-errors)/overlap<0.8500000238418579:return None,0,True
        return forced_merge(left,right,delta,reference),2*overlap-7*errors,False
    merged=merge_pair(left['sequence'],left['quality'],reverse_complement(right['sequence']),right['quality'][::-1],minimal_overlap=20,minimal_identity=.95)
    if merged is None:return None,0,False
    return forced_merge(left,right,merged['offset'],reference),merged['score'],False


def order_targets(targets,reference):
    best={}
    for kind in 'vdjc':
        scores={}
        for target in targets:
            for hit in target[kind+'_hits']:scores[hit['gene']]=scores.get(hit['gene'],0)+hit['score']
        if scores:best[kind]=max(_trove_order(scores,reference.gene_hashes),key=scores.get)
    ordered=[]
    for target in targets:
        for index,kind in enumerate('vdjc'):
            hit=next((h for h in target[kind+'_hits'] if h['gene']==best.get(kind)),None)
            if hit:ordered.append(((index,hit['ref_start']),target));break
    return [t for _,t in sorted(ordered,key=lambda item:item[0])]


def merge_ancillary(left,right,left_id,right_id,delta,begin,end,reference):
    central=forced_merge(left[left_id],right[right_id],delta,reference)
    central['overlap_range']=[begin,end]
    groups=[left[:left_id]+right[:right_id],left[left_id+1:]+right[right_id+1:]]
    # The source intentionally keeps the best score across successive merges.
    best_score=0
    for group in groups:
        while True:
            best=None
            for i,target in enumerate(group):
                merged,score,conflict=merge_targets(target,central,reference)
                if conflict:return None
                if score>best_score:best=(i,merged);best_score=score
            if best is None:break
            i,central=best;group.pop(i)
    for group in groups:
        while True:
            change=False
            for i in range(len(group)):
                for j in range(i+1,len(group)):
                    merged,score,conflict=merge_targets(group[i],group[j],reference)
                    if conflict:return None
                    if merged:
                        group[i]=merged;group.pop(j);change=True;break
                if change:break
            if not change:break
    return order_targets(groups[0]+[central]+groups[1],reference)
