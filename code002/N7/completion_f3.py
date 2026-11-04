"""The two remaining F3 fibres: complete compatibility by binomial Gram equations."""
from common import *
from torus_solver import solve_torus
from minor_solver import solve as minor_solve

def mul(x,y):
    a,b=x;c,d=y
    return a*c-b*d,a*d+b*c+b*d

def main():
    results=[]
    for q in [Q(1,6),Q(5,6)]:
        fs=fibre('F3',(q,q));assert len(fs)==4 and all(len(c[2])==1 for c in fs)
        rows=[];self_products=[]
        for i,(c,B,p) in enumerate(fs):
            assert c[:2]==(q,q) and B[:2]==((0,),(0,))
            cc=c[2:];bb=B[2:];assert sum(row[0] for row in bb)==2 and all(row[0] in [0,1] for row in bb)
            # Distinct members of one circle cannot be orthogonal: 5+2t/s.
            left=tuple((row[0],0) for row in bb);right=tuple((0,row[0]) for row in bb)
            selfgram=gram(cc,left,cc,right);assert selfgram=={(0,0):(5,0),(1,-1):(2,0)}
            self_products.append([dict(degree=k,coefficient=v) for k,v in selfgram.items()])
            rows.append((cc,tuple(tuple(row[0] if j==i else 0 for j in range(4)) for row in bb)))
        # Thus a square must choose one row from each circle, with four phases.
        equations=[];rhs=[];pair_products=[]
        for i,j in combinations(range(4),2):
            g=gram(*rows[i],*rows[j]);assert len(g)==2
            (u,A),(v,B)=list(g.items())
            roots=[k for k in range(6) if mul(A,POWERS[k])==tuple(-a for a in B)]
            assert len(roots)==1
            row=tuple(a-b for a,b in zip(u,v));b=Q(roots[0],6)
            equations.append(row);rhs.append(b)
            pair_products.append(dict(rows=[i+1,j+1],terms=[dict(degree=k,coefficient=v) for k,v in g.items()],
                                      equation=row,rhs=str(b)))
        A=tuple(equations);b=tuple(rhs)
        sols=maximal({canon(c,B) for c,B in solve_torus(A,b,4)})
        other=maximal({canon(c,B) for c,B in minor_solve(A,b,4)})
        assert sols==other and len(sols)==1
        sol=sols[0];expected=canonical((Q(0),Q(1,2),(q+Q(1,2))%1,q),((1,),(1,),(-1,),(-1,)))
        assert sol==expected
        pc,PB,p=sol;assert len(p)==1
        cs,Bs=seed_rows('F3',(q,q),((0,),(0,)))
        for c,B in rows:
            cs.append(tuple((a+sum(x*y for x,y in zip(row,pc)))%1 for a,row in zip(c,B)))
            Bs.append(tuple(tuple(sum(row[h]*PB[h][j] for h in range(4)) for j in range(len(p))) for row in B))
        cycles=verify_matrix(cs,Bs);w=f4_match(cs,Bs);assert w is not None
        results.append(dict(theta=str(q),row_circles=[serialize(c) for c in fs],same_circle_gram=self_products,
                            pair_grams=pair_products,solutions=[serialize(s) for s in sols],
                            independent_solver_match=True,matrix=serialize_matrix(cs,Bs),cycles=cycles,f4_witness=w))
        print('F3',q,'unique completion circle',serialize(sol),'F4 witness',w,flush=True)
    (ROOT/'results/F3_completions.json').write_text(json.dumps(results,indent=2)+'\n')
    print('Both F3 seeds: all orthogonal square completions form one circle each; all are regular and contain F4.')
if __name__=='__main__':main()
