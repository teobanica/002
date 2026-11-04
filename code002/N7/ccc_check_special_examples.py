"""Direct verification of the special configurations (x) and (y).

Standalone standard-library checker. No import of the classification solvers,
support enumeration, canonicalization or family-containment code.
All computations are in Z/24Z. All 3!*7! row/column permutations are visited.
"""
from itertools import permutations,product
from pathlib import Path
import json,time

ROOT=Path(__file__).parent
M=24
def phase_rows(u,v):
    return [[0]*7,[x%M for x in u],[x%M for x in v]]
def dephased_flat(H):
    return tuple((H[r][j]-H[0][j]-H[r][0]+H[0][0])%M for r in (1,2) for j in range(7))

def old_cases():
    for w in (8,16):
        for xi in range(24):
            if 4*xi%24!=(12+w)%24:continue
            for e in (0,1):
                z=(1+e)*w+2*xi
                H=phase_rows([0,w,2*w,xi,12+z,12+xi,z],
                             [0,12+2*w-xi,z,w,2*w,2*w-xi,12+z])
                yield {'case':'x','w_exponent':w,'xi_exponent':xi,'e':e},H
        for s,t in product((1,-1),repeat=2):
            eta=6+2*w+(12 if s==-1 else 0)
            rho=6+w+(12 if t==-1 else 0)
            H=phase_rows([0,w,2*w,w+eta,2*w+rho,12+w+eta,12+2*w+rho],
                         [0,eta,rho,w,2*w,12+eta,12+rho])
            yield {'case':'y','w_exponent':w,'eta_sign':s,'rho_sign':t},H

def family(name,p):
    w=8
    if name=='F1':
        a,b,c=p
        return phase_rows([0,w,2*w,a,12+a,a,12+a],[0,2*w,w,b,c,12+b,12+c])
    if name=='F2':
        a,b=p
        return phase_rows([12,w,2*w,0,0,a,12+a],[12,2*w,w,b,12+b,0,0])
    if name=='F3':
        a,b=p
        return phase_rows([12+a+b,w,2*w,0,12+a,a,a+b],[12+a+b,2*w,w,b,a+b,0,12+b])
    a=p[0]
    if name=='F4':
        return phase_rows([0,w,2*w,w+a,12+w+a,a,12+a],[0,w,2*w,a,w+a,12+a,12+w+a])
    if name=='F5':
        return phase_rows([0,w,2*w,w-a,12+w-a,2*w+a,12+2*w+a],[0,12+w+a,a,w,2*w,w+a,12+a])
    if name=='F6':
        return phase_rows([0,w,2*w,a,12+a,-a,12-a],[0,12+2*w-a,w+a,w,12+w+a,2*w-a,2*w])
    raise ValueError(name)

def exceptional():
    A=[0,8,16,1,13,2,14];B=[0,8,16,1,13,10,22];C=[0,8,16,2,14,6,18]
    inputs=[
    (A,[0,2,11,8,23,14,16]),(A,[0,2,23,11,8,14,16]),
    (A,[0,6,5,17,14,12,22]),(A,[0,6,17,14,5,12,22]),
    (B,[0,10,7,19,16,22,8]),(B,[0,10,19,16,7,22,8]),
    (B,[0,22,9,6,21,12,14]),(B,[0,22,21,9,6,12,14]),
    (C,[0,4,2,10,16,18,12]),(C,[0,16,2,4,10,18,12]),
    (C,[0,16,14,22,4,12,6])]
    return {dephased_flat(phase_rows(u,v)):f'E{i}' for i,(u,v) in enumerate(inputs,1)}

def recover(y):
    u=y[:7];v=y[7:]
    # Recover the original parameters of the six displayed families from entries.
    return [('F1',(u[3],v[3],v[4])),
            ('F2',((u[5]+12)%24,(v[3]+12)%24)),
            ('F3',((12-v[3])%24,-u[4]%24)),
            ('F4',((u[3]-8)%24,)),
            ('F5',((8-u[3])%24,)),
            ('F6',(u[3],))]

def transformed(H,rows,cols):
    return [[H[r][j] for j in cols] for r in rows]

def check():
    E=exceptional()
    cases=list(old_cases())
    assert len(cases)==24
    records=[]
    for info,H in cases:
        matches={};visited=0
        for rows in permutations(range(3)):
            base,a,b=rows
            for anchor in range(7):
                U=[(H[a][j]-H[base][j]-H[a][anchor]+H[base][anchor])%24 for j in range(7)]
                V=[(H[b][j]-H[base][j]-H[b][anchor]+H[base][anchor])%24 for j in range(7)]
                for tail in permutations([j for j in range(7) if j!=anchor]):
                    visited+=1
                    # Every literal target has u_3/u_2=w, with w=exp(2*pi*i/3).
                    if (U[tail[1]]-U[tail[0]])%24!=8:continue
                    cols=(anchor,)+tail
                    y=tuple(U[j] for j in cols)+tuple(V[j] for j in cols)
                    if y in E:
                        name=E[y]
                        matches.setdefault(name,{'target':name,'rows':[r+1 for r in rows],
                                                 'columns':[j+1 for j in cols]})
                    for name,params in recover(y):
                        target=dephased_flat(family(name,params))
                        if y!=target:continue
                        cert={'target':name,'parameter_exponents':params,
                              'rows':[r+1 for r in rows],'columns':[j+1 for j in cols]}
                        matches.setdefault(name,cert)
        assert visited==30240 and matches
        for cert in matches.values():
            rows=[r-1 for r in cert['rows']];cols=[j-1 for j in cert['columns']]
            image=dephased_flat(transformed(H,rows,cols))
            if cert['target'].startswith('F'):
                target=dephased_flat(family(cert['target'],cert['parameter_exponents']))
            else:
                target=next(v for v,k in E.items() if k==cert['target'])
            assert image==target
        names=set(matches)
        if info['case']=='x':
            assert len(names)==1 and next(iter(names)) in {f'E{i}' for i in range(1,9)}
        else:
            assert names==({'F6'} if info['eta_sign']==info['rho_sign'] else {'F5'})
        rec={**info,'permutations_checked':visited,'matches':list(matches.values())}
        records.append(rec)
        print(info,'->',rec['matches'],flush=True)
    xs=[r for r in records if r['case']=='x']
    assert {r['matches'][0]['target'] for r in xs}=={f'E{i}' for i in range(1,9)}
    assert len(xs)==16
    # Human-readable identities, checked for every sign without permutations of rows.
    for s in (1,-1):
        same=next(H for info,H in cases if info=={'case':'y','w_exponent':8,'eta_sign':s,'rho_sign':s})
        opposite=next(H for info,H in cases if info=={'case':'y','w_exponent':8,'eta_sign':s,'rho_sign':-s})
        a=6 if s==1 else 18
        assert transformed(same,(0,1,2),(0,1,2,3,6,5,4))==family('F6',(a,))
        assert opposite==family('F5',((a+20)%24,))
    result={'status':'passed','source_cases':24,'x_cases':16,'y_cases':8,
            'permutations_checked':sum(r['permutations_checked'] for r in records),
            'records':records}
    (ROOT/'special_example_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: all 24 choices covered; x -> E1,...,E8; y same signs -> F6, opposite signs -> F5.')
    return result
if __name__=='__main__':
    start=time.perf_counter();check();print('Elapsed',round(time.perf_counter()-start,2),'seconds.')

