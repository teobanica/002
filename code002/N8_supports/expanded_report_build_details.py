"""Generate readable exact tables and certificates, checking every printed square."""
import json,sys
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import defaultdict,Counter
R=Path(__file__).resolve().parents[1]
P=R.parent if (R.parent/'AAA_class21').exists() else R.parent/'n8_bundle/proofs'
O=R/'report';O.mkdir(exist_ok=True)
sys.path.insert(0,str(P/'consolidation/code'))
from exact_rows import profile,profiles,quotient
def read(part,name):return json.loads((P/part/'outputs'/name).read_text())
def frac(f):return str(f.numerator) if f.denominator==1 else r'\tfrac{%d}{%d}'%(f.numerator,f.denominator)
def mono(e,names):
    c=F(e[0])%1
    factor='' if c==0 else '-' if c==F(1,2) else r'\e\!\left('+frac(c)+r'\right)'
    terms=''.join(n if k==1 else n+'^{'+str(k)+'}' for n,k in zip(names,e[1:]) if k)
    return factor+(terms or '1') if factor in ['', '-'] else factor+terms
def mat(H,names):
    return r'{\small\setlength{\arraycolsep}{5pt}\['+'\n'+r'\begin{pmatrix}'+'\n'+r'\\'.join('&'.join(mono(e,names) for e in row) for row in H)+'\n'+r'\end{pmatrix}\]}'+'\n'
def decode(H,m):return [[(F(e[0],m),)+tuple(e[1:]) for e in row] for row in H]
checks=Counter()
def check(H):
    assert len(H)==8 and all(len(r)==8 for r in H)
    assert all(profile(quotient(a,b),[2]*4) for a,b in combinations(H,2))
    checks['printed_symbolic_A_squares']+=1
def out(name,lines):(O/name).write_text('\n'.join(lines)+'\n')
cs=read('AAA_class21','components.json');sc=read('AAA_class21','generic_AAA.json')
table=[r'\begin{longtable}{rrrrl}\toprule Class&Dimension&Components&Rows per branch&Edges per branch\\\midrule\endhead']
atlas=[]
for cid in range(23):
    subset=[(c,q) for c,q in zip(cs,sc) if c['class_id']==cid]
    table.append(f"{cid}&{subset[0][0]['dim']}&{len(subset)}&"+','.join(str(len(q['rows'])) for c,q in subset)+'&'+','.join(str(q['edges']) for c,q in subset)+r'\\')
    for c,q in subset:
        if not q['k5']:continue
        H=decode(c['rows']+[q['rows'][j] for j in q['k5']],30);check(H)
        names=list('abcd')[:c['dim']]
        atlas.extend([r'\Needspace{15\baselineskip}',r'\paragraph{Class %d, branch %d.}'%(cid,c['branch']),
          'The free unit parameters are $'+','.join(names)+r'$. The first three rows are the registered Smith component.',mat(H,names)])
table.append(r'\bottomrule\end{longtable}');out('aaa_scan_table.tex',table);out('uniform_certificates.tex',atlas)
for cid in [21,22]:
    part=f'AAA_class{cid}';atlas=[]
    for i,x in enumerate(read(part,f'class{cid}_square_templates.json'),1):
        H=decode(x['matrix'],x['constant_modulus']);check(H)
        o=[F(*v) for v in x['parameter_origin']];v=x['parameter_direction']
        xy=[((o[1]-o[0])%1,v[1]-v[0]),(o[1],v[1])] if cid==21 else [(o[0],v[0]),((o[1]+F(1,2))%1,v[1])]
        # Recover the registered base rows from the same raw coordinates.
        c=next(c for c in cs if c['class_id']==cid and c['branch']==x['branch'])
        base=[[( (F(e[0],30)+sum(e[k+1]*o[k] for k in range(2)))%1,sum(e[k+1]*v[k] for k in range(2))) for e in row] for row in c['rows']]
        assert base==H[:3]
        atlas.extend([r'\Needspace{11\baselineskip}',r'\paragraph{Certificate %d.%d.}'%(cid,i),
          r'$\epsilon=%d,\quad x=%s,\quad y=%s$, with $q\in\T$. Add the following five rows to $T_{%d,\epsilon}(x,y)$:'%(x['epsilon'],mono(xy[0],['q']),mono(xy[1],['q']),cid),mat(H[3:],['q'])])
    out(f'circle{cid}_certificates.tex',atlas)
    residues=read(part,'final_residue.json')
    sets=[]
    for d in residues:
        # final_residue uses raw Smith (s,t); table uses natural (a,b) for 22.
        pts=[tuple(F(*e) for e in p) for p in d['points']]
        if cid==22:pts=[(s,(t+F(1,2))%1) for s,t in pts]
        sets.append(set(pts))
    assert sets[0]==sets[1]
    M=24 if cid==21 else 12;g=defaultdict(list)
    for s,t in sorted(sets[0]):
        assert (s*M).denominator==(t*M).denominator==1
        g[int(s*M)].append(int(t*M))
    tab=[r'\begin{longtable}{rl}\toprule $i$&Values of $j$\\\midrule\endhead']
    for i,jj in g.items():tab.append(str(i)+'&'+', '.join(map(str,jj))+r'\\')
    tab.append(r'\bottomrule\end{longtable}');out(f'residue{cid}.tex',tab)
atlas=[]
for kind in ['good_planes','good_lines']:
 for i,x in enumerate(read('AAA_classes11_19',kind+'.json'),1):
    tor=x['torus'];o=[F(*e) for e in tor['o']];V=tor['V'];dim=tor['dim'];names=list('ab')[:dim]
    H=decode(x['square'],x['modulus']) if kind=='good_planes' else [[(F(*e[0]),)+tuple(e[1:]) for e in row] for row in x['square']]
    check(H)
    base=[[( (F(e[0],30)+sum(e[k+1]*o[k] for k in range(3)))%1,)+tuple(sum(e[k+1]*V[k][j] for k in range(3)) for j in range(dim)) for e in row] for row in x['H']]
    assert base==H[:3]
    xyz=[mono((o[k],)+tuple(V[k]),names) for k in range(3)]
    label='surface' if dim==2 else 'circle'
    atlas.extend([r'\Needspace{12\baselineskip}',r'\paragraph{Class %d, %s certificate %d.}'%(x['class_id'],label,i),
      r'$x=%s,\quad y=%s,\quad z=%s$, with $%s\in\T$. Add the five rows below to $T_{%d}(x,y,z)$:'%(*xyz,','.join(names),x['class_id']),mat(H[3:],names)])
out('three_parameter_certificates.tex',atlas)
# Preserve all surface equations used before the 3D converse descends to curves.
def eq(A,c):
    terms=[]
    for a,v in zip(A,['s','t','u']):
        if a:terms.append(('+' if a>0 and terms else '-' if a<0 else '')+(str(abs(a)) if abs(a)!=1 else '')+v)
    return ''.join(terms)+'='+frac(c)
goods=read('AAA_classes11_19','good_planes.json')
goodkeys={json.dumps(g['torus'],sort_keys=True) for g in goods}
for cid in [11,19]:
 tab=[r'\begin{longtable}{ll}\toprule Surface equation&Generic completion\\\midrule\endhead']
 planes=[p for p in read('AAA_classes11_19','planes.json') if p['class_id']==cid]
 scans=[s for s in read('AAA_classes11_19','planes_scan.json') if s['class_id']==cid]
 assert len(planes)==len(scans)
 for p,s in zip(planes,scans):
    t=p['torus'];assert len(t['A'])==1
    tab.append('$'+eq(t['A'][0],F(*t['c'][0]))+'$&'+('A-only square identity' if s['k5'] else 'No generic square')+r'\\')
 tab.append(r'\bottomrule\end{longtable}');out(f'planes{cid}.tex',tab)
# All nine BBB generic fourth-row formulas in their one-parameter coordinates.
atlas=[]
cs_pure=json.loads((P/'consolidation/inputs/pure_components.json').read_text())
BBB=next(c['rows'] for c in cs_pure if c['name']=='BBB')
for i,x in enumerate(read('consolidation','BBB_circle_templates.json'),1):
    t=x['torus'];o=[F(*e) for e in t['o']];V=t['V']
    def specialize(row):return [((F(e[0],60)+sum(e[k+1]*o[k] for k in range(3)))%1,sum(e[k+1]*V[k][0] for k in range(3))) for e in row]
    rows=[specialize(row) for row in x['rows_mod_60']];base=[specialize(row) for row in BBB];assert len(rows)==1
    assert all(profiles(quotient(rows[0],h)) for h in base)
    a,b=[mono((o[k],V[k][0]),['q']) for k in [0,1]]
    atlas.extend([r'\Needspace{7\baselineskip}',r'\paragraph{BBB circle %d.}'%i,
      r'$a=%s,\quad b=%s$, with $q\in\T$. A fourth row, unique away from the finite special set, is'%(a,b),mat(rows,['q'])])
out('bbb_rows.tex',atlas);checks['BBB_symbolic_rows']=9
# Reproduce one entire CCC ten-row graph and its square transformations.
old=json.loads((P/'consolidation/inputs/n8_ccc_square_results.json').read_text());cert=read('consolidation','ccc_canonical.json')
point=next(p for p in old['points'] if (p['A'],p['B'])==(2,14));rows=[r['row'] for r in point['candidate_fourth_rows']]
tab=[r'\begin{center}\small\begin{tabular}{r|rrrrrrrr|c}\toprule Row&\multicolumn{8}{c|}{Exponents modulo $30$}&Profile\\\midrule']
for i,row in enumerate(rows):
    typ=''.join(profiles([(F(a-b,30),) for a,b in zip(row,h)]) for h in point['old_rows_exp30'])
    tab.append(str(i+1)+'&'+'&'.join(map(str,row))+'&'+typ+r'\\')
tab.append(r'\bottomrule\end{tabular}\end{center}');out('ccc_ten_rows.tex',tab)
chosen=[(i,s) for i,s in enumerate(old['square_completions']) if s['point']==point['name']];assert len(chosen)==2
lines=[]
for n,(i,s) in enumerate(chosen,1):
    c=cert[i];H=s['rows_exp30'];r0,c0=c['pivot']
    assert [[(H[a][b]-H[r0][b]-H[a][c0]+H[r0][c0])%30 for b in c['columns']] for a in c['rows']]==c['canonical']
    members=', '.join(str(j+1) for j in s['candidate_clique'])
    lines.append(r'The %s component is $\{%s\}$. Its square is sent to the displayed canonical matrix by pivot $(r,c)=(%d,%d)$, row order $%s$, and column order $%s$.'%('first' if n==1 else 'second',members,r0+1,c0+1,str(tuple(j+1 for j in c['rows'])),str(tuple(j+1 for j in c['columns']))))
out('ccc_example_certificates.tex',lines);checks['CCC_example_equivalences']=2
assert checks['printed_symbolic_A_squares']==92
(R/'checks/generated_detail_checks.json').write_text(json.dumps(dict(checks),indent=2))
print('PASS:',dict(checks))
