"""Smith-split every square locus and verify all pairwise prime cycles exactly."""
from pathlib import Path
import sys,json
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'n8_aaa_complete/code'))
from phase_algebra import *
from functools import lru_cache
from math import lcm
HOME_ROOT=R
R=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HOME_ROOT
strata=json.loads((R/'outputs/strata.json').read_text())['strata']
@lru_cache(None)
def checked_profile(q):
    p=profiles(q);assert p and 'C' not in p
    return p

MOD=json.loads((R/'outputs/strata.json').read_text())['modulus']
def decode60(row):return tuple((F(e[0],MOD),)+tuple(map(F,e[1:])) for e in row)
def solve60(eq,n):return solutions([e[:-1]+[S.Rational(30*e[-1],MOD)] for e in eq],n,0)
seen={};out=[];counts=[];modulus=MOD
for s in strata:
    p=R/f"outputs/cliques_stratum_{s['id']:02d}.json"
    if not p.exists():
        print('Not yet finished:',s['id'],flush=True);continue
    loci=json.loads(p.read_text())[0]['squares'];fs=s['families'];count=0
    for li,loc in enumerate(loci):
        chosen=[fs[i] for i in loc['families']];nf=sum(f['free'] for f in chosen);n=nf+s['base']
        for bi,values in enumerate(solve60(loc['eq'],n)):
            d=len(values[0])-1 if values else 0
            base_values=values[nf:]
            H=[substitute(decode60(row),0,base_values,d) for row in s['rows']]
            pos=0
            for f in chosen:
                H.append(substitute(decode60(f['row']),0,base_values+values[pos:pos+f['free']],d));pos+=f['free']
            assert len(H)==8 and all(e[0]%1==0 and not any(e[1:]) for e in H[0])
            assert all(not any(row[0]) for row in H)
            key=tuple(H[:3])+tuple(sorted(H[3:]));origin=dict(stratum=s['id'],locus=li,torsion_branch=bi)
            count+=1
            if key in seen:out[seen[key]]['origins'].append(origin);continue
            assert d<=5,(s['id'],li,d)
            ps=[]
            for i,j in combinations(range(8),2):
                q=tuple(tuple(e) for e in quotient(H[i],H[j]));ps.append(checked_profile(q))
            for row in H:
                for e in row:modulus=lcm(modulus,e[0].denominator)
            seen[key]=len(out);out.append(dict(square=enc(H),dimension=d,origins=[origin],profiles=ps))
    counts.append(dict(stratum=s['id'],loci=len(loci),smith_branches=count))
    print(s['id'],'loci',len(loci),'branches',count,'cumulative unique',len(out),flush=True)
(R/'outputs/completion_squares.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
(R/'outputs/square_extraction_summary.json').write_text(json.dumps(dict(modulus=modulus,counts=counts,unique=len(out),dimensions=dict(Counter(f['dimension'] for f in out)),all_strata_finished=len(counts)==len(strata)),indent=2)+'\n')
print('Unique',len(out),'modulus',modulus,'dimensions',dict(Counter(f['dimension'] for f in out)),flush=True)
targets=json.loads((HOME_ROOT.parent/'n8_aaa_complete/outputs/ten_representatives.json').read_text())
from prepare_affine import write_input
write_input(R/'inputs/completion_embeddings.txt',targets,out,modulus)
