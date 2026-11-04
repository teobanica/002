"""Compact parameter conditions, regression, certificates and readable catalogue."""
from fractions import Fraction as Q
from itertools import combinations,product
from pathlib import Path
from collections import Counter,defaultdict
import json,hashlib,sys
from families import family
from census import deserialize,canonical,contains,solve_cached
from independent import canon
from analyze import four_rows,f4_match,valid_patterns
ROOT=Path(__file__).resolve().parent

def maximal(cs):
    out=[]
    for dim in sorted({len(c[2]) for c in cs},reverse=True):
        layer=[c for c in cs if len(c[2])==dim]
        out.extend(c for c in layer if not any(contains(b,c) for b in out))
    return set(out)

def point(q):return canonical(q,[() for _ in q])
def curve(c,v):return canon(c,[(x,) for x in v])

def expected_base(name):
    out=[]
    if name=='F1':
        for k in range(6):
            for j in range(6):
                if j in ((2*k)%6,(2*k+3)%6):continue
                if k%2==0:out.append(curve([Q(k,6),Q(j,6),0],[0,0,1]))
                else:out.append(curve([Q(k,6),0,Q(j,6)],[0,1,0]))
    elif name=='F2':
        for j in [1,2,4,5]:
            out.extend([curve([Q(j,6),0],[0,1]),curve([0,Q(j,6)],[1,0])])
    elif name=='F3':
        for j in [1,5]:out.extend([curve([Q(j,6),0],[0,1]),curve([0,Q(j,6)],[1,0])])
        for j in [2,4]:out.append(curve([0,Q(j,6)],[1,-1]))
        out.extend(point([Q(j,6),Q(k,6)]) for j,k in product([2,3,4],repeat=2) if j!=k)
    elif name=='F4':out=[curve([0],[1])]
    else:out=[point([Q(j,6)]) for j in range(6)]
    return set(out)

def check_witness(cs,Bs,w,t=None):
    if t is not None:
        cs=[[ (c+sum(b)*t)%1 for c,b in zip(row,bb)] for row,bb in zip(cs,Bs)]
        Bs=[[() for _ in range(7)] for _ in range(4)]
        ac=Q(w['f4_parameter']);ab=()
    else:ac=Q(w['parameter_constant']);ab=tuple(w['parameter_basis'])
    rr=[j-1 for j in w['rows']];cc=[j-1 for j in w['columns']];base=rr[0];anchor=cc[0]
    target=family('F4')
    for h,r in enumerate(rr[1:],1):
        for k,j in enumerate(cc):
            c=(cs[r][j]-cs[base][j]-cs[r][anchor]+cs[base][anchor])%1
            b=tuple(Bs[r][j][l]-Bs[base][j][l]-Bs[r][anchor][l]+Bs[base][anchor][l] for l in range(len(ab)))
            e=target['constants'][h][k];v=target['basis'][h][k][0]
            assert (c-Q(e,6)-v*ac)%1==0 and b==tuple(v*x for x in ab)

def fibre(comp,theta,d):
    c,B,p=comp;A=tuple(B[:d]);rhs=tuple((Q(v)-c[i])%1 for i,v in enumerate(theta))
    for pc,PB in solve_cached(A,rhs,len(p)):
        k=len(PB[0]) if PB else 0
        cc=[(x+sum(a*b for a,b in zip(row,pc)))%1 for x,row in zip(c,B)]
        BB=[tuple(sum(row[h]*PB[h][j] for h in range(len(p))) for j in range(k)) for row in B]
        yield canon(cc,BB)

def evaluate_fibre(name,theta):
    data=json.loads((ROOT/'results'/(name+'.json')).read_text());s=family(name);d=s['parameters']
    fs=maximal({x for r in data['maximal_components'] for x in fibre(deserialize(r),theta,d)})
    out=[]
    for comp in fs:
        cs,Bs=four_rows(s,comp);w=f4_match(cs,Bs)
        out.append(dict(dimension=len(comp[2]),constant=[str(v) for v in comp[0]],basis=comp[1],f4=bool(w)))
    return dict(theta=[str(Q(t)) for t in theta],components=out,finite_rows=sum(r['dimension']==0 for r in out),
                circles=sum(r['dimension']==1 for r in out),not_f4_finite_rows=sum(r['dimension']==0 and not r['f4'] for r in out))

def monomial(c,n=0):
    e=int(6*Q(c))%6; assert Q(e,6)==Q(c)%1
    base='1' if e==0 else ('omega' if e==1 else 'omega^'+str(e))
    if n==0:return base
    t='t' if n==1 else ('t^-1' if n==-1 else 't^'+str(n))
    return t if base=='1' else base+'*'+t

def main():
    lines=['COMPLETE FOURTH-ROW CENSUS FOR THE SIX AFFINE SEEDS',
           'omega = exp(2*pi*i/6); every t below ranges independently over the whole unit circle.',
           'Literal seed formulas appear below. Each record gives its seed parameters and new row x.',
           'The first new entry is 1. Labels denote fixed-seed joint solution components, not matrix equivalence classes.',
           'F4=all means every matrix in this component contains an equivalent F4 triple.',
           'F4=finite lists exactly the t values with that property; F4=none means none.', '']
    reference=json.loads((ROOT/'input_reference.json').read_text())['normalized_families']
    for name in ['F1','F2','F3','F4','F5','F6']:
        s=family(name);c=[];B=[]
        for row,bb in zip(s['constants'][1:],s['basis'][1:]):
            c.extend(Q(v-row[0],6)%1 for v in row)
            B.extend(tuple(x-y for x,y in zip(v,bb[0])) for v in bb)
        old=reference[name]
        assert canonical(c,B)==(tuple(Q(v,24) for v in old['constant_24ths']),tuple(tuple(v) for v in old['basis']),tuple(old['pivot']))
        for i,j in combinations(range(3),2):
            c=[Q(x-y,6)%1 for x,y in zip(s['constants'][i],s['constants'][j])]
            B=[tuple(x-y for x,y in zip(a,b)) for a,b in zip(s['basis'][i],s['basis'][j])]
            assert valid_patterns(c,B)
        lines.append(name+' seed (first row is all ones):')
        for row,bb in zip(s['constants'][1:],s['basis'][1:]):
            entries=[]
            for e,b in zip(row,bb):
                factors=[monomial(Q(e,6))] if e%6 else []
                factors += [('abc'[k] if v==1 else 'abc'[k]+'^'+str(v)) for k,v in enumerate(b) if v]
                entries.append('*'.join(factors) or '1')
            lines.append('  ('+', '.join(entries)+')')
        lines.append('')
    verification=[]
    for name in ['F1','F2','F3','F4','F5','F6']:
        data=json.loads((ROOT/'results'/(name+'.json')).read_text());an=json.loads((ROOT/'results'/(name+'_analysis.json')).read_text());f4=json.loads((ROOT/'results'/(name+'_f4_loci.json')).read_text());s=family(name);d=s['parameters']
        actual=maximal({canon(deserialize(r)[0][:d],deserialize(r)[1][:d]) for r in data['maximal_components']})
        assert actual==expected_base(name),(name,len(actual),len(expected_base(name)))
        inds=json.loads((ROOT/'results'/(name+'_independent.json')).read_text());assert inds['match']
        for i,(r,a,f) in enumerate(zip(data['maximal_components'],an['branches'],f4['branches']),1):
            comp=deserialize(r);cs,Bs=four_rows(s,comp)
            # Independently sum Laurent coefficients in Q(omega).
            powers=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
            for u,v in combinations(range(4),2):
                sums=defaultdict(lambda:[0,0])
                for j in range(7):
                    phase=(cs[u][j]-cs[v][j])%1;assert (phase*6).denominator==1
                    deg=tuple(x-y for x,y in zip(Bs[u][j],Bs[v][j]));z=powers[int(6*phase)]
                    for k in range(2):sums[deg][k]+=z[k]
                assert all(z==[0,0] for z in sums.values())
            if a['f4_witness'] is not None:check_witness(cs,Bs,a['f4_witness'])
            for w in f['special']:check_witness(cs,Bs,w,Q(w['t']))
            c,B,p=comp
            entry=[monomial(c[j],B[j][0] if p else 0) for j in range(d+7)]
            lines += [f'{name}-{i:03d}: '+('one circle' if p else 'one point'),
                      '  seed: '+', '.join('abc'[j]+'='+entry[j] for j in range(d)),
                      '  x = ('+', '.join(entry[d:])+')',
                      '  F4='+f['status']+(' at t in {'+', '.join(monomial(Q(z['t'])) for z in f['special'])+'}' if f['special'] else ''),'']
        verification.append(dict(name=name,components=len(data['maximal_components']),compact_parameter_conditions_match=True,direct_laurent_checks=True,f4_certificates_checked=True))
    # F4 pilot, independently recorded as its seven explicit formulas.
    P=[[0,4,3,4,2,0,1],[0,4,1,0,4,2,3],[0,5,2,2,4,2,5],[0,3,2,4,0,4,1],
       [0,0,0,3,2,4,3],[0,3,1,5,5,3,2],[0,5,3,3,3,1,0]]
    expected=set()
    for i,e in enumerate(P):
        b=[0]*7 if i==4 else [0,0,0,1,0,0,1]
        expected.add(canonical([Q(0)]+[Q(v,6) for v in e],[(1,)]+[(v,) for v in b]))
    assert expected=={deserialize(r) for r in json.loads((ROOT/'results/F4.json').read_text())['maximal_components']}
    # F1: 24 disjoint base circles, five rows per seed, two wholly F4.
    an=json.loads((ROOT/'results/F1_analysis.json').read_text());by={r['id']:r for r in an['branches']}
    assert len(an['base_loci'])==24
    for g in an['base_loci']:
        assert g['projection_dimension']==1 and len(g['branches'])==5
        assert sum(by[i]['f4_witness'] is not None for i in g['branches'])==2
    # F5/F6: count exact finite fibres and rows not already F4.
    fibres56={n:[evaluate_fibre(n,[Q(j,6)]) for j in range(6)] for n in ['F5','F6']}
    assert all(x['not_f4_finite_rows']<4 and x['circles']==0 for row in fibres56.values() for x in row)
    # F3: inspect all intersections of the six base curves and all extra point loci.
    an=json.loads((ROOT/'results/F3_analysis.json').read_text());curves=[g for g in an['base_loci'] if g['projection_dimension']==1]
    lookup={r['id']:r for r in an['branches']}
    assert len(curves)==6
    assert all(len(g['branches'])==4 and sum(lookup[i]['f4_witness'] is not None for i in g['branches'])==2 for g in curves)
    intersections=set()
    for g,h in combinations(curves,2):
        A=tuple(tuple(row) for row in g['equations']+h['equations']);b=tuple(Q(x) for x in g['rhs']+h['rhs'])
        for c,B in solve_cached(A,b,2):
            assert len(B[0])==0;intersections.add(tuple(c))
    extra={tuple(Q(x) for x in g['rhs']) for g in an['base_loci'] if g['projection_dimension']==0}
    special=[evaluate_fibre('F3',p) for p in sorted(intersections|extra)]
    unresolved=[r for r in special if r['circles'] or r['not_f4_finite_rows']>=4]
    D1={(Q(x,6),Q(y,6)) for x,y in [(0,j) for j in [1,2,4,5]]+[(j,0) for j in [1,2,4,5]]+[(1,2),(2,1),(4,5),(5,4)]}
    D2={(Q(x,6),Q(y,6)) for pool in ([2,3,4],[1,3,5]) for x,y in product(pool,repeat=2) if x!=y}
    D3={(Q(1,6),Q(1,6)),(Q(5,6),Q(5,6))}
    assert intersections|extra==D1|D2|D3 and len(D1)==len(D2)==12
    for f in special:
        pt=tuple(Q(v) for v in f['theta'])
        if pt in D1:assert (f['finite_rows'],f['circles'],f['not_f4_finite_rows'])==(5,0,2)
        elif pt in D2:assert (f['finite_rows'],f['circles'],f['not_f4_finite_rows'])==(7,0,0)
        else:assert pt in D3 and f['circles']==4 and f['finite_rows']==0
    reduction=dict(F1='all square completions contain an F4 triple by 5-row-list / 3-unknown-row pigeonhole',
                   F5_F6_fibres=fibres56,F3_base_curve_intersections=[[str(x) for x in t] for t in sorted(intersections)],
                   F3_special_fibres=special,F3_unresolved_seed_points=unresolved)
    (ROOT/'catalogue.txt').write_text('\n'.join(lines)+'\n')
    (ROOT/'verification.json').write_text(json.dumps(dict(checks=verification,F4_pilot_match=True,reduction=reduction),indent=2)+'\n')
    print('All six compact parameter conditions verified; all 333 components directly checked.')
    print('F4 pilot and every recorded F4 equivalence certificate verified.')
    print('F1,F4,F5,F6 square sectors covered by F4 or by the row-count argument.')
    print('F3 base intersections:',len(intersections),'unresolved seed points:',[(r['theta'],r['finite_rows'],r['circles'],r['not_f4_finite_rows']) for r in unresolved])
if __name__=='__main__':main()
