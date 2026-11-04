"""All F2 square completions, including every exceptional seed fibre."""
from common import *
from collections import Counter

def main():
    groups=defaultdict(list);points=set()
    for idx,rec in enumerate(data('F2')['maximal_components'],1):
        comp=deserialize(rec);c,B,p=comp
        if len(p)==0:points.add(c[:2]);continue
        assert len(p)==1 and any(row[0] for row in B[:2])
        base=canonical(c[:2],B[:2]);assert base[0]==c[:2] and base[1]==B[:2]
        groups[base].append((idx,comp))
    assert len(groups)==8 and all(len(g)==4 for g in groups.values())
    circles=[]
    for base,records in sorted(groups.items()):
        cs,Bs=seed_rows('F2',base[0],base[1])
        cs.extend(comp[0][2:] for idx,comp in records);Bs.extend(comp[1][2:] for idx,comp in records)
        cycles=verify_matrix(cs,Bs);w=f4_match(cs,Bs)
        circles.append(dict(base=serialize(base),fourth_row_ids=[idx for idx,c in records],
                            matrix=serialize_matrix(cs,Bs),cycles=cycles,f4_witness=w))
        print('F2 circle',base[0],base[1],'F4 witness',w,flush=True)
        assert w is not None
    # Derive every intersection; use exact affine congruences, not a grid.
    bases=list(groups)
    for i,j in combinations(range(len(bases)),2):
        c,B,p=bases[j];fixed=[h for h in range(2) if B[h]==(0,)]
        for comp in restrict(bases[i],fixed,[c[h] for h in fixed]):
            assert len(comp[2])==0
            points.add(comp[0])
    S={Q(j,6) for j in [1,2,4,5]}
    expected={tuple(Q(j,6) for j in q) for q in product(range(6),repeat=2) if Q(q[0],6) in S or Q(q[1],6) in S}
    assert points==expected and len(points)==32
    finite=[]
    for theta in sorted(points):
        fs=fibre('F2',theta);assert all(len(c[2])==0 for c in fs)
        cs0,Bs0=seed_rows('F2',theta,[(),()])
        cs=[c[0][2:] for c in fs];Bs=[c[1][2:] for c in fs]
        edges=[]
        for i,j in combinations(range(len(fs)),2):
            orth=not gram(cs[i],Bs[i],cs[j],Bs[j])
            reg=bool(valid_patterns(tuple((a-b)%1 for a,b in zip(cs[i],cs[j])),tuple(() for _ in range(7))))
            assert orth==reg
            if reg:edges.append((i,j))
        edge_set=set(edges)
        cliques=[s for s in combinations(range(len(fs)),4) if all(p in edge_set for p in combinations(s,2))]
        specializations={}
        for label,circle in enumerate(circles,1):
            bc,bB,bp=deserialize(circle['base']);free=bp[0]
            t=(theta[free]-bc[free])%1
            if any((bc[h]+bB[h][0]*t-theta[h])%1 for h in range(2)):continue
            matrix=circle['matrix']
            specialized=frozenset(tuple((Q(c)+row[0]*t)%1 for c,row in zip(rc,rb))
                                  for rc,rb in zip(matrix['constants'][3:],matrix['basis'][3:]))
            assert len(specialized)==4
            specializations.setdefault(specialized,[]).append(label)
        assert set(specializations)=={frozenset(cs[i] for i in s) for s in cliques}
        # Every edge belongs to one of these square cliques. This fixes the
        # entire graph, including the isolated fifth row in the five-row fibres.
        assert edge_set=={p for s in cliques for p in combinations(s,2)}
        if len(fs)==7:assert len(cliques)==2 and len(set(cliques[0])&set(cliques[1]))==1
        graph_type={4:'K4',5:'K4 plus isolated vertex',7:'two K4 graphs sharing one vertex'}[len(fs)]
        expected_n=5 if any(t not in S for t in theta) else (4 if (2*(theta[0]-theta[1]))%1==0 else 7)
        assert len(fs)==expected_n
        completions=[]
        for selected in cliques:
            cc=cs0+[cs[i] for i in selected];BB=Bs0+[Bs[i] for i in selected]
            cycles=verify_matrix(cc,BB);w=f4_match(cc,BB);assert w is not None,(theta,selected)
            completions.append(dict(selected_rows=[i+1 for i in selected],f4_witness=w,cycles=cycles,
                                    circle_specializations=specializations[frozenset(cs[i] for i in selected)]))
        finite.append(dict(theta=[str(t) for t in theta],rows=[[str(v) for v in row] for row in cs],
                           edges=[[i+1,j+1] for i,j in edges],graph_type=graph_type,completions=completions))
        print('F2 seed',tuple(str(t) for t in theta),'rows',len(fs),'edges',len(edges),'squares',len(cliques),'all F4',flush=True)
    out=dict(circles=circles,exceptional_seeds=finite,
             summary=dict(circle_completions=len(circles),exceptional_seed_count=len(finite),
                          finite_fibre_histogram={str(k):v for k,v in Counter((len(r['rows']),len(r['completions'])) for r in finite).items()},
                          finite_square_count=sum(len(r['completions']) for r in finite),all_completions_have_f4=True))
    (ROOT/'results/F2_completions.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['summary'],flush=True)
if __name__=='__main__':main()
