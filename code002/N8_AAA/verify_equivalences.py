"""Check every saved equivalence entry by entry, independently of the search."""
from phase_algebra import *
targets=read('generic_cover_targets.json')
sources=json.loads((R/'inputs/equivalence_expanded_sources.json').read_text())
certs=read('generic_equivalences_expanded.json');assert len(certs)==len(sources)==1519
def matrix(t):return [[tuple(F(*v) for v in e) for e in r] for r in t['square']]
for k,(s,c) in enumerate(zip(sources,certs)):
    assert c['source']==k and c['target']>=0
    H=matrix(s);T=matrix(targets[c['target']]);p=c['rows'];q=c['columns'];ds=len(H[0][0])-1
    assert sorted(p)==sorted(q)==list(range(8))
    values=[(F(e[0],2),)+tuple(map(F,e[1:])) for e in c['parameter_map']]
    for i in range(8):
        for j in range(8):
            lhs=tuple(H[p[i]][q[j]][a]-H[p[0]][q[j]][a]-H[p[i]][q[0]][a]+H[p[0]][q[0]][a] for a in range(ds+1))
            rhs=[T[i][j][0]]+[F(0)]*ds
            for a,v in enumerate(values):
                for b,x in enumerate(v):rhs[b]+=T[i][j][a+1]*x
            assert (lhs[0]-rhs[0])%1==0 and lhs[1:]==tuple(rhs[1:]),(k,i,j)
print('PASS: all 1519 affine-family embeddings checked entry by entry with rational phases.')
cert=certs[2]
assert cert['target']==1
A=S.Matrix([x[1:] for x in cert['parameter_map']]);assert abs(A.det())==1
print('PASS: the two equivalent five-parameter templates have an invertible integer parameter change.')
(R/'outputs/equivalence_check.json').write_text(json.dumps({'checked_embeddings':len(certs),'uncovered':0,'duplicate_maximal_template_map_determinant':int(A.det())},indent=2)+'\n')
