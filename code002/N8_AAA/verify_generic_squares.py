"""Independently parameterize and verify every generic square locus."""
from phase_algebra import *

ds=read('generic_reduced.json');loci=read('generic_square_loci.json')
cs=json.loads((R/'inputs/components.json').read_text());allout=[];summary=[]
for ci,(d,l,c) in enumerate(zip(ds,loci,cs)):
    assert l['component_index']==ci
    fs=[d['families'][i] for i in d['maximal_family_indices']];base=d['dim'];templates=[]
    counts=Counter();profiles_seen=Counter()
    for li,x in enumerate(l['squares']):
        selected=[fs[i] for i in x['families']];nf=sum(f['free'] for f in selected)
        for branch,values in enumerate(solutions(x['eq'],nf,base)):
            extra=(len(values[0])-base-1) if nf else 0
            H=[pad(decode(r),extra) for r in c['rows']];pos=0
            for fam in selected:
                f=fam['free'];H.append(substitute(decode(fam['row']),base,values[pos:pos+f],extra));pos+=f
            ps=[profiles(quotient(a,b)) for a,b in combinations(H,2)]
            assert all(ps),(ci,li,branch)
            assert all(all(r[0][j]==0 for j in range(len(r[0]))) for r in H)
            dim=S.Matrix([e[1:] for r in H for e in r]).rank()
            assert dim==base+extra,(ci,li,dim,base,extra)
            counts[dim]+=1;profiles_seen[tuple(sorted(Counter(ps).items()))]+=1
            templates.append({'locus_index':li,'torsion_branch':branch,'dimension':dim,'square':H,'row_profiles':dict(Counter(ps))})
    out={'class_id':c['class_id'],'branch':c['branch'],'component_index':ci,'base_dimension':base,'templates':templates}
    (R/f'outputs/generic_squares_{ci:02d}.json').write_text(json.dumps(enc(out),separators=(',',':'))+'\n')
    summary.append({'class_id':c['class_id'],'branch':c['branch'],'loci':len(l['squares']),'templates':len(templates),'dimensions':dict(counts),'row_profile_counts':[{'counts':dict(k),'templates':v} for k,v in profiles_seen.items()]})
    print(c['class_id'],c['branch'],'verified templates',len(templates),'dimensions',dict(counts),'profiles',dict(profiles_seen),flush=True)
(R/'outputs/generic_square_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS: every square locus split into all torsion branches and every quotient checked exactly.')
