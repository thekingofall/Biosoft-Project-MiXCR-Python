"""Python port of MiLib 1.13 KAlignmentHit and linear band extensions."""
NEGATIVE = -1073741824
MISSING = -2147483647


def band_align(reference, query, width=2, mode="global", added1=0, added2=0, reverse=False):
    """Origin-anchored DP, with MiLib's unequal-length band and tie order."""
    if reverse:
        reference=reference[::-1];query=query[::-1]
    m,n=len(reference),len(query)
    width=min(width,m,n)
    rowfactor=2*width+abs(n-m)
    delta=max(0,m-n)+width
    h={(0,0):0}
    for j in range(1,rowfactor-delta+1):h[0,j]=-12*j
    for i in range(1,delta+1):h[i,0]=-12*i
    maximum=0;stop_i=stop_j=0
    for i in range(m):
        rowmax=NEGATIVE
        for j in range(max(0,i-delta),min(i+rowfactor-delta+1,n)):
            score=max(h.get((i,j),NEGATIVE)+(5 if reference[i]==query[j] else -9),
                      h.get((i,j+1),NEGATIVE)-12,h.get((i+1,j),NEGATIVE)-12)
            h[i+1,j+1]=score
            if score>maximum:maximum=score;stop_i=i+1;stop_j=j+1
            rowmax=max(rowmax,score)
        if mode=='semi' and rowmax-i*5 < -1000:break
    if mode=='global':
        stop_i,stop_j=m,n
    elif mode=='added':
        maximum=NEGATIVE
        for i in range(m-added1,m+1):
            score=h.get((i,n),NEGATIVE)
            if score>maximum:maximum=score;stop_i=i;stop_j=n
        for j in range(n-added2,n+1):
            score=h.get((m,j),NEGATIVE)
            if score>maximum:maximum=score;stop_i=m;stop_j=j
    i,j=stop_i,stop_j
    operations=[]
    while i or j:
        score=h.get((i,j),NEGATIVE)
        if i and j and score==h.get((i-1,j-1),NEGATIVE)+(5 if reference[i-1]==query[j-1] else -9):
            operations.append('M' if reference[i-1]==query[j-1] else 'X');i-=1;j-=1
        elif i and score==h.get((i-1,j),NEGATIVE)-12:
            operations.append('D');i-=1
        elif j and score==h.get((i,j-1),NEGATIVE)-12:
            operations.append('I');j-=1
        else:raise RuntimeError('Unreachable banded traceback cell')
    if not reverse:operations.reverse()
    return stop_i,stop_j,''.join(operations)


def _trim(operations,reverse=False):
    work=operations[::-1] if reverse else operations
    score=0;minimum=1;stop=0
    for index,op in enumerate(work,1):
        score+=5 if op=='M' else -9 if op=='X' else -12
        if score<=0 and score<=minimum:minimum=score;stop=index
    if not stop:return operations,''
    return (operations[:-stop],operations[-stop:]) if reverse else (operations[stop:],operations[:stop])


def extend_hit(query,hit,start=0,end=None,kind='V',floating_left=None,floating_right=None):
    """Convert one KMapper hit to a scored, traced gene alignment."""
    end=len(query) if end is None else end
    floating_left=kind!='C' if floating_left is None else floating_left
    floating_right=kind=='V' if floating_right is None else floating_right
    reference=hit['gene']['sequence']
    positions=hit['seed_positions'][hit['from']:hit['to']]
    anchors=[(p,o) for p,o in zip(positions,hit['seed_offsets']) if o!=MISSING]
    if not anchors:return None
    p,offset=anchors[0]
    qfrom=offset-2;qlen=p-qfrom;qadd=5
    rfrom=start-offset-2;rlen=p-offset-rfrom;radd=5
    if qfrom<start:
        change=start-qfrom;qlen-=change;qadd=max(0,min(qlen,qadd-change));qfrom=start
    if rfrom<0:
        rlen+=rfrom;radd=max(0,min(rlen,radd+rfrom));rfrom=0
    ri,qi,ops=band_align(reference[rfrom:rfrom+rlen],query[qfrom:qfrom+qlen],
                         mode='semi' if floating_left else 'added',added1=radd,added2=qadd,reverse=True)
    refstart=rfrom+rlen-ri;querystart=qfrom+qlen-qi
    operations=ops+'M'*5
    previous_p,previous_offset=p,offset
    refend=p-offset+5;queryend=p+5
    for p,offset in anchors[1:]:
        qfrom=previous_p+5;qlen=p-qfrom
        rfrom=qfrom-previous_offset;rlen=p-offset-rfrom
        if rlen<0 or qlen<0:
            qfrom-=5;rfrom-=5;qlen+=5;rlen+=5
            if not operations.endswith('MMMMM'):raise RuntimeError('KMapper seed overlap violated')
            operations=operations[:-5]
        if rlen>0 or qlen>0:
            _,_,ops=band_align(reference[rfrom:rfrom+rlen],query[qfrom:qfrom+qlen])
            operations+=ops
        operations+='M'*5
        refend=p-offset+5;queryend=p+5
        previous_p,previous_offset=p,offset
    qfrom=queryend;qlen=len(reference)-refend+2;qadd=5
    rfrom=refend;rlen=end-queryend+2;radd=5
    if qfrom+qlen>end:
        qadd-=qfrom+qlen-end;qlen=end-qfrom;qadd=max(0,min(qadd,qlen))
    if rfrom+rlen>len(reference):
        radd-=rfrom+rlen-len(reference);rlen=len(reference)-rfrom;radd=max(0,min(radd,rlen))
    ri,qi,ops=band_align(reference[rfrom:rfrom+rlen],query[qfrom:qfrom+qlen],
                         mode='semi' if floating_right else 'added',added1=radd,added2=qadd)
    refend=rfrom+ri;queryend=qfrom+qi;operations+=ops
    if floating_right:
        operations,trim=_trim(operations,reverse=True)
        refend-=sum(op!='I' for op in trim);queryend-=sum(op!='D' for op in trim)
    if floating_left:
        operations,trim=_trim(operations)
        refstart+=sum(op!='I' for op in trim);querystart+=sum(op!='D' for op in trim)
    r,q=refstart,querystart;mapping={r:q}
    for op in operations:
        if op!='I':r+=1
        if op!='D':q+=1
        mapping[r]=q
    if (r,q)!=(refend,queryend):raise RuntimeError('Banded alignment coordinate mismatch')
    score=5*operations.count('M')-9*operations.count('X')-12*(operations.count('I')+operations.count('D'))
    return {'gene':hit['gene']['name'],'score':score,'ref_start':refstart,'ref_end':refend,
            'query_start':querystart,'query_end':queryend,'operations':operations,'boundary_map':mapping}


def shift_homopolymer_indels(hit, reference, query):
    """MiXCR CommandAlign postprocessing: shift equivalent gaps to the left."""
    mutations=[]
    r,q=hit['ref_start'],hit['query_start']
    for op in hit['operations']:
        if op!='M':mutations.append([r,op,reference[r] if op=='D' else query[q]])
        if op!='I':r+=1
        if op!='D':q+=1
    previous=hit['ref_start']
    for mutation in mutations:
        position,op,base=mutation
        if op=='X':previous=position+1;continue
        while position>previous and reference[position-1]==base:position-=1
        mutation[0]=position
        previous=position+(op=='D')
    r=hit['ref_start'];operations=[]
    for position,op,base in mutations:
        operations.extend('M'*(position-r));operations.append(op)
        r=position+(op!='I')
    operations.extend('M'*(hit['ref_end']-r))
    hit['operations']=''.join(operations)
    r,q=hit['ref_start'],hit['query_start'];mapping={r:q}
    for op in operations:
        if op!='I':r+=1
        if op!='D':q+=1
        mapping[r]=q
    hit['boundary_map']=mapping
    return hit
