"""Rebuild the seven representative families and exact maximality certificates."""
from phase_algebra import *
from math import factorial

def dec(H):return [[tuple(F(*x) for x in e) for e in row] for row in H]
def matching(q):
    c=Counter(q);ans=1
    while c:
        e=min(c);n=c.pop(e);o=((e[0]+F(1,2))%1,)+e[1:]
        assert c.pop(o,0)==n
        ans*=factorial(n)
    return ans

def invariant(H):
    return [dict(Counter(matching(tuple(quotient(a,b))) for a,b in combinations(K,2))) for K in (H,list(zip(*H)))]

targets=read('maximal_square_templates.json')
targets.append({'class_id':17,'row_family':'transpose_additional_4','base_dimension':4,'square':[list(x) for x in zip(*targets[4]['square'])]})
targets.append({'class_id':16,'row_family':'additional_2','base_dimension':2,'square':read('generic_squares_27.json')['templates'][0]['square']})
(R/'outputs/generic_cover_targets.json').write_text(json.dumps(targets,indent=2))
targets.append({'class_id':None,'row_family':'transpose_additional_2','base_dimension':2,'square':[list(x) for x in zip(*targets[6]['square'])]})
(R/'outputs/row_column_closed_targets.json').write_text(json.dumps(targets,indent=2))
certs=[];invariants=[]
for ti,t in enumerate(targets):
    H=dec(t['square']);d=len(H[0][0])-1
    assert all(profiles(quotient(a,b))=='A' for K in (H,list(zip(*H))) for a,b in combinations(K,2))
    eq=[]
    for i,j in combinations(range(8),2):
        q=quotient(H[i],H[j]);groups={}
        for c,e in enumerate(q):groups.setdefault(e[1:],[]).append((c,1 if e[0]==0 else -1));assert e[0] in [0,F(1,2)]
        for group in groups.values():
            assert sum(s for c,s in group)==0
            row=[0]*49
            for c,s in group:
                if c and i:row[(i-1)*7+c-1]+=s
                if c and j:row[(j-1)*7+c-1]-=s
            if any(row):eq.append(row)
    A=S.Matrix(eq);V=S.Matrix([H[i][j][1:] for i in range(1,8) for j in range(1,8)])
    assert A*V==S.zeros(A.rows,d) and V.rank()==d
    ir=list(A.T.rref()[1]);ic=list(A[ir,:].rref()[1]);det=A[ir,ic].det()
    assert len(ir)==len(ic)==49-d and det
    certs.append({'target':ti,'dimension':d,'rank':len(ir),'equations':eq,'minor_rows':ir,'minor_columns':ic,'minor_determinant':int(det)})
    invariants.append({'target':ti,'dimension':d,'row_column_matching_distributions':invariant(H)})
    print('Family',ti,'dimension',d,'rank',len(ir),'minor',det,'invariants',invariant(H),flush=True)
assert invariants[1]['row_column_matching_distributions']==invariants[2]['row_column_matching_distributions']
seven=[invariants[i] for i in [0,1,3,4,5,6,7]]
assert len({str(t['row_column_matching_distributions']) for t in seven})==7
# Same exact certificate as the preliminary computation, up to equation sign/order.
(R/'outputs/family_rank_certificates.json').write_text(json.dumps(certs,indent=2))
(R/'outputs/all_matching_invariants.json').write_text(json.dumps(invariants,indent=2))
(R/'outputs/seven_families.json').write_text(json.dumps([{'name':n,'target_index':i,**targets[i]} for n,i in zip(['F5a','F5b','F5c','F4','F4T','F2','F2T'],[0,1,3,4,5,6,7])],indent=2))
