"""Exact local regular locus at a root-valued dephased square.
At a fixed point each feasible prime-cycle support gives homogeneous linear
phase equations. Their finite union is the exact local regular locus. The
algorithm first intersects the row spaces common to all supports per pair,
then enumerates the remaining alternatives inside this smaller space.
"""
from phase_algebra import *
from functools import lru_cache

def cycle_partitions(q, parts):
    def rec(indices,ps):
        if not ps:
            yield ();return
        a=indices[0]
        for p in sorted(set(ps)):
            for others in combinations(indices[1:],p-1):
                block=(a,)+others
                if isinstance(q[a],tuple):
                    if any(q[j][1:]!=q[a][1:] for j in block):continue
                    if sorted((q[j][0]-q[a][0])%1 for j in block)!=[F(k,p) for k in range(p)]:continue
                elif sorted((q[j]-q[a])%1 for j in block)!=[F(k,p) for k in range(p)]:continue
                rest=tuple(j for j in indices if j not in block);pp=list(ps);pp.remove(p)
                for tail in rec(rest,tuple(pp)):yield tuple(sorted((block,)+tail))
    return sorted(set(rec(tuple(range(8)),tuple(parts))))

def basis(M):
    if not M.rows:return S.zeros(0,M.cols)
    B,piv=M.T.rref();return M[list(piv),:] if piv else S.zeros(0,M.cols)

def canonical(M):
    if not M.rows:return S.zeros(0,M.cols)
    B,piv=M.rref();return B[:len(piv),:]

def key(M):return tuple(tuple(F(x) for x in M.row(i)) for i in range(M.rows))

def null(M):
    ns=M.nullspace();return S.Matrix.hstack(*ns) if ns else S.zeros(M.cols,0)

def local_locus(H):
    formal=isinstance(H[0][0],(list,tuple))
    H=[[tuple([F(x[0])%1]+list(map(F,x[1:]))) if formal else F(x)%1 for x in row] for row in H]
    records=[];common=[]
    for i,j in combinations(range(8),2):
        q=quotient(H[i],H[j]) if formal else [(a-b)%1 for a,b in zip(H[i],H[j])]
        parts=[]
        for tag,ps in [('A',[2]*4),('B',[2,3,3]),('C',[3,5])]:
            parts += [(tag,p) for p in cycle_partitions(q,ps)]
        assert parts
        mats=[]
        for tag,part in parts:
            rows=[]
            for b in part:
                for c in b[1:]:
                    v=[0]*8;v[c]=1;v[b[0]]=-1;rows.append(v)
            mats.append(S.Matrix(rows))
        # A row equation common to all support spaces annihilates every nullspace.
        C=null(S.Matrix.vstack(*(null(A).T for A in mats))).T
        D=S.zeros(8,49)
        for c in range(1,8):
            if i:D[c,(i-1)*7+c-1]=1
            if j:D[c,(j-1)*7+c-1]=-1
        common.append(C*D)
        records.append((i,j,parts,mats,D))
    G=canonical(S.Matrix.vstack(*common));N=null(G)
    print('Common equations rank',G.rows,'remaining',N.cols,flush=True)
    alternatives=[]; blockers=[]
    for i,j,parts,mats,D in records:
        reduced=[A*D*N for A in mats]
        if N.cols==1 and all(any(x for x in A) for A in reduced):
            cert=[]
            for (tag,part),A,E in zip(parts,mats,reduced):
                k=next(k for k in range(E.rows) if E[k,0])
                row=(A*D)[k,:]
                det=S.Matrix.vstack(G,row).det()
                assert det
                cert.append({'profile':tag,'cycles':part,'row':enc([F(x) for x in row]),'determinant':enc(F(det))})
            blockers.append({'pair':[i,j],'supports':cert})
        unique={key(canonical(A)):canonical(A) for A in reduced}
        alternatives.append((len(unique),i,j,list(unique.values()),parts))
    alternatives.sort(key=lambda x:x[0])
    spaces={():S.zeros(0,N.cols)}
    counts=[]
    for ct,i,j,alts,parts in alternatives:
        new={}
        for L in spaces.values():
            for A in alts:
                B=canonical(S.Matrix.vstack(L,A));new[key(B)]=B
        # Keep maximal solution spaces (minimal constraint spaces).
        keep=[]
        for B in sorted(new.values(),key=lambda x:x.rows):
            if any(canonical(S.Matrix.vstack(B,A)).rows==B.rows for A in keep):continue
            keep.append(B)
        spaces={key(B):B for B in keep}
        counts.append([i,j,len(parts),ct,len(spaces)])
    out=[]
    for B in spaces.values():
        V=N*null(B);dim=V.cols
        # Clear each column denominator: same local subspace, global torus.
        for c in range(dim):
            den=S.ilcm(*[x.q for x in V[:,c]])
            V[:,c]*=den
        sq=[]
        for i in range(8):
            row=[]
            for j in range(8):
                v=[F(0)]*dim if i==0 or j==0 else [F(x) for x in V[(i-1)*7+j-1,:]]
                row.append(((H[i][j][0] if formal else H[i][j]),)+tuple(v))
            sq.append(row)
        typ=Counter(profiles(quotient(sq[i],sq[j])) for i,j in combinations(range(8),2))
        assert '' not in typ
        out.append({'dimension':dim,'matrix':enc(sq),'profiles':dict(typ)})
    return {'common_rank':G.rows,'common_equations':enc([[F(x) for x in row] for row in G.tolist()]),'common_nullspace':enc([[F(x) for x in row] for row in N.tolist()]),'pair_counts':counts,'one_parameter_blockers':blockers,'families':out}

if __name__=='__main__':
    x=read('mixed_finite_squares.json')[0]
    result=local_locus([[F(v,12) for v in row] for row in x['matrix_mod12']])
    (R/'outputs/mixed_local_locus.json').write_text(json.dumps(result,indent=2))
    print('Exact local families:',[(f['dimension'],f['profiles']) for f in result['families']])
