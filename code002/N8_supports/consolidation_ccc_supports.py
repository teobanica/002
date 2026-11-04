import itertools,collections,math
from functools import lru_cache
def quotient(r,s): return tuple((a-b)%30 for a,b in zip(r,s))

def is_rot_cycle(vals,p):
    step=30//p; vals=sorted(v%30 for v in vals)
    return any(sorted((r+k*step)%30 for k in range(p))==vals for r in range(step))

@lru_cache(None)
def A_supports(q):
    q=tuple(q); out=[]
    def rec(rem,pairs):
        if not rem:
            out.append(tuple(sorted(tuple(sorted(x)) for x in pairs))); return
        a=min(rem)
        for b in sorted(rem):
            if b>a and (q[a]-q[b])%30==15:
                rec(rem-{a,b},pairs+[(a,b)])
    rec(set(range(8)),[])
    return tuple(sorted(set(out)))

@lru_cache(None)
def B_supports(q):
    q=tuple(q); out=set(); I=range(8)
    for pair in itertools.combinations(I,2):
        if (q[pair[0]]-q[pair[1]])%30 != 15: continue
        rem=[i for i in I if i not in pair]
        for tri in itertools.combinations(rem,3):
            tri2=tuple(i for i in rem if i not in tri)
            if is_rot_cycle([q[i] for i in tri],3) and is_rot_cycle([q[i] for i in tri2],3):
                out.add((tuple(sorted((tuple(sorted(tri)),tuple(sorted(tri2))))),tuple(sorted(pair))))
    return tuple(sorted(out))

@lru_cache(None)
def C_supports(q):
    q=tuple(q); out=set(); I=range(8)
    for five in itertools.combinations(I,5):
        tri=tuple(i for i in I if i not in five)
        if is_rot_cycle([q[i] for i in five],5) and is_rot_cycle([q[i] for i in tri],3):
            out.add((tuple(sorted(five)),tuple(sorted(tri))))
    return tuple(sorted(out))

def bits(q):
    z=0
    if A_supports(tuple(q)): z|=1
    if B_supports(tuple(q)): z|=2
    if C_supports(tuple(q)): z|=4
    return z

def masks_for_blocks(r,s,blocks):
    out=[]
    for block in blocks:
        a=block[0]
        for c in block[1:]:
            m=0
            for rr,cc in ((r,c),(s,c),(r,a),(s,a)):
                if rr and cc: m ^= 1<<((rr-1)*7+(cc-1))
            if m: out.append(m)
    return tuple(out)

def rank2(masks):
    piv=[0]*49; rk=0
    for x in masks:
        while x:
            p=x.bit_length()-1
            if piv[p]: x^=piv[p]
            else:
                piv[p]=x; rk+=1
                if rk==49:return 49
                break
    return rk

def tangent_check(rows):
    base=[]; opt=[]
    for r,s in itertools.combinations(range(8),2):
        q=quotient(rows[r],rows[s]); b=bits(q)
        if b==1:
            aa=A_supports(q); assert len(aa)==1
            base += masks_for_blocks(r,s,aa[0])
        elif b==4:
            cc=C_supports(q); assert len(cc)==2
            opt.append((masks_for_blocks(r,s,cc[0]),masks_for_blocks(r,s,cc[1])))
        else: raise AssertionError(('unexpected profile',b))
    assert len(opt)==12
    for choice in range(1<<12):
        m=list(base)
        for k,o in enumerate(opt):m.extend(o[(choice>>k)&1])
        assert rank2(m)==49
    return 1<<12

