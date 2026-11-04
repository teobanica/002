"""Typeset exact matrices, support witnesses, and census tables from verified data."""
from local_regular import *
import math
P=R/'report'
labels=['F_{5a}','F_{5b}','F_{5c}','F_4',r'F_4^{\rm T}','F_2',r'F_2^{\rm T}','K_2','I_{12}','M_1']
raw=read('ten_representatives.json');shown=list(raw)
for i in [7,9]:shown[i]=read(f'centered_family_{i}.json')
def Hdecode(t):return [[tuple(F(*x) for x in e) for e in row] for row in t['square']]
def monomial(e,names='abcde'):
    c=e[0]%1;sgn=''
    if c>=F(1,2):sgn='-';c-=F(1,2)
    const=''
    if c:
        den=c.denominator;k=c.numerator;const=r'\zeta_{'+str(den)+'}' + ('^{'+str(k)+'}' if k!=1 else '')
    pos=[];neg=[]
    for n,k in zip(names,e[1:]):
        assert k.denominator==1
        if not k:continue
        dest=pos if k>0 else neg;v=abs(k);dest.append(n if v==1 else n+'^{'+str(v)+'}')
    top=const+''.join(pos) or '1';bottom=''.join(neg)
    return sgn+(r'\frac{'+top+'}{'+bottom+'}' if bottom else top)
def matrix(H,names='abcde'):
    return r'\begin{pmatrix}'+'\n'+(((' '+chr(92)*2+'\n').join(' & '.join(monomial(e,names) for e in row) for row in H)))+'\n'+r'\end{pmatrix}'
def display_matrix(H,name=None,names='abcde'):
    body=(name+'='+matrix(H,names)) if name else matrix(H,names)
    return '\n'+r'\begin{center}\small\setlength{\arraycolsep}{3.4pt}'+'\n'+r'\begin{adjustbox}{max width=\textwidth}$\displaystyle '+body+r'$\end{adjustbox}'+'\n'+r'\end{center}'+'\n'
texts={
0:r'This five-parameter family comes from the first free-row construction over maximal component $0$. It is ordinarily equivalent to its transpose. The common regular variation matrix has rank $44$, and the five parameter directions span its kernel.',
1:r'This is the second independent family from maximal component $0$. A third presentation obtained there is equivalent to it by a unimodular parameter substitution. Its transpose is equivalent to $F_{5c}$, while the two generic types are distinct under ordinary monomial equivalence.',
2:r'This five-parameter family comes from the free-row construction over maximal component $1$. The first three rows retain the prescribed triple and the fifth phase comes from an independent component rotation of an added row. Its ordered row and column invariants are those of $F_{5b}$ with the directions exchanged.',
3:r'This is the additional four-parameter completion exposed by the critical compatibility values in maximal component $1$. Its generic row quotients have no four-against-four opposite value split. Hence it cannot be a specialization of any of the former five-parameter templates. The common variation matrix has rank $45$.',
4:r'The literal transpose supplies the second four-parameter type. Its row and column invariants are reversed from those of $F_4$ and therefore distinguish it generically. The generic completion outputs from class $17$ and several lower-dimensional components embed in this family.',
5:r'This two-parameter family occurs as the unique generic square parameterization over one component of class $16$. Its generic row matching distribution includes one value $24$, while its column distribution has none. The common regular variation equations have rank $47$.',
6:r'The transpose gives a second two-parameter type. It is a valid regular family with an AAA triple at every point, but it is not needed by the generic prescribed-triple cover. Its existence is an explicit warning that changing the AAA frame does not automatically move a family into that generic calculation.',
7:r'This family is extracted from exact local branch calculations at catalogue matrices missed by the first seven-family cover. The constants have been shifted to a real center by an exact parameter translation. All thirty-eight positive-dimensional A branch presentations found at the missed catalogue points embed here. Its common regular variation rank is $47$.',
9:r'Here $\zeta_3=\e(1/3)$ and $\zeta_6=\e(1/6)$; both notations may occur when a negative sign is separated from a root. The centered formula uses only sixth-root constants. Its generic profile is $24$ A-only quotients and $4$ B-only quotients in both row and column directions. Its regular local dimension is one.'}
out=[];cert=[]
for i in [0,1,2,3,4,5,6,7,9]:
    H=Hdecode(shown[i]);d=len(H[0][0])-1
    assert all(profiles(quotient(a,b)) for K in (H,list(zip(*H))) for a,b in combinations(K,2))
    out += [r'\Needspace{22\baselineskip}\subsection{The family \texorpdfstring{$'+labels[i]+r'$}{Family '+str(i)+r'}}', texts[i], display_matrix(H,labels[i]+'('+','.join('abcde'[:d])+')'),r'The parameters are independent: the dephased exponent matrix has rank $'+str(d)+r'$. A complete symbolic row-pair cycle certificate is printed in Appendix~\ref{app:cycles}.']
    pairs=[]
    for a,b in combinations(range(8),2):
        q=quotient(H[a],H[b]);tag=profiles(q)[0];ps={'A':[2]*4,'B':[3,3,2],'C':[3,5]}[tag];part=cycle_partitions(q,ps)[0]
        pairs.append((f'{a+1},{b+1}',r'\;'.join(''.join(str(c+1) for c in block) for block in part)))
    cert += [r'\subsection{Cycle certificate for \texorpdfstring{$'+labels[i]+r'$}{Family '+str(i)+r'}}',r'\begin{center}\small\begin{tabular}{cl@{\qquad}cl}\toprule Pair&Cycle blocks&Pair&Cycle blocks\\\midrule']
    for a,b in zip(pairs[:14],pairs[14:]):cert.append('$'+a[0]+'$ & $'+a[1]+'$ & $'+b[0]+'$ & $'+b[1]+r'$ \\')
    cert.append(r'\bottomrule\end{tabular}\end{center}')
(P/'families.tex').write_text('\n\n'.join(out)+'\n');(P/'cycle_certificates.tex').write_text('\n'.join(cert)+'\n')
# Original natural coordinates for the corrected F4.
H=Hdecode(raw[3]);K=[]
for row in H:
    rr=[]
    for e in row:
        c,a,b,z,d=e;rr.append(((c+d/2)%1,a,b,-z,d-a-b))
    K.append(rr)
assert all(profiles(quotient(a,b))=='A' for a,b in combinations(K,2))
assert [monomial(e,'xyzw') for e in K[1]]==['1','-1','-x','x','y','-y','1','-1']
assert [monomial(e,'xyzw') for e in K[2]]==['1','-z','-w','-w','w','w','-1','z']
(P/'natural_F4.tex').write_text(display_matrix(K[3:],names='xyzw'))
# Numeric exponent matrices.
for filename,E in [('isolated_matrix.tex',read('mixed_canonical.json')[0]['canonical']),('BH18_matrix.tex',read('mixed_circle_analysis.json')['BH18_matrix'])]:
    tex='\\[\n\\begin{pmatrix}\n'+(' '+chr(92)*2+'\n').join(' & '.join(map(str,row)) for row in E)+'\n\\end{pmatrix}.\n\\]\n';(P/filename).write_text(tex)
# Generic component census.
rows=[r'\subsection{Component-by-component census}',r'The branch labels are those of the earlier AAA registry. ``Loci'' counts the retained integer relation systems before Smith splitting; ``presentations'' counts their explicit torsion-split square parameterizations.',r'\begin{longtable}{rrr r r l}\toprule Class&Branch&Base dim.&Loci&Presentations&Square dimensions\\\midrule\endfirsthead',r'\toprule Class&Branch&Base dim.&Loci&Presentations&Square dimensions\\\midrule\endhead']
cs=json.loads((R/'inputs/components.json').read_text())
for c,x in zip(cs,read('generic_square_summary.json')):
    dim=len(c['rows'][0][0])-1;desc=', '.join(f'{k}: {v}' for k,v in sorted(x['dimensions'].items(),key=lambda z:int(z[0]))) or '--'
    rows.append(f"{x['class_id']}&{x['branch']}&{dim}&{x['loci']}&{x['templates']}&{desc}"+r'\\')
rows.append(r'\bottomrule\end{longtable}');(P/'generic_table.tex').write_text('\n'.join(rows))
# Matching invariants.
rows=[r'\begin{center}\small\begin{tabular}{cll}\toprule Type&Row distribution&Column distribution\\\midrule']
for x in read('ten_family_invariants.json'):
    ds=[]
    for inv in x['matching_invariants']:ds.append(r',\;'.join(str(k)+'^{['+str(v)+']}' for k,v in sorted(inv.items(),key=lambda z:int(z[0]))))
    rows.append('$'+labels[x['representative']]+'$ & $'+ds[0]+'$ & $'+ds[1]+r'$\\')
rows.append(r'\bottomrule\end{tabular}\end{center}');rows.append(r'Here $m^{[r]}$ means that exactly $r$ unordered pairs have $m$ antipodal matchings. The zero entries record B-only pairs.')
(P/'invariant_table.tex').write_text('\n'.join(rows))
# Catalogue coverage.
a=read('aalto_final_summary.json');assert a['checked']==a['covered']==5469
rows=[r'\begin{center}\begin{tabular}{rrrrr}\toprule $q$&Catalogue classes&Regular&With AAA&Identified\\\midrule']
for x in a['by_q']:rows.append(f"{x['q']}&{x['catalogue_entries']}&{x['regular']}&{x['AAA']}&{x['covered']}"+r'\\')
rows += [r'\midrule Total entries&5469&5469&5469&5469\\\bottomrule\end{tabular}\end{center}',r'At $q=12$, exactly one catalogue class fails to admit A on every row pair; it is $I_{12}$. At the other seven root orders, every catalogue entry admits A on every row pair.']
(P/'catalogue_table.tex').write_text('\n'.join(rows))
# Navigable output guide.
items=[('generic_all.json','Complete generic fourth-row families, including disconnected supports.'),('generic_reduced.json','Exact row-family containment substitutions.'),('reduced_pairs.json','All generic compatibility congruence systems.'),('generic_square_loci.json','Five-row compatibility loci before torsion splitting.'),('generic_squares_00.json through generic_squares_37.json','All 1512 explicit squares, with dimensions and profiles.'),('generic_equivalences_expanded.json','Permutations and parameter substitutions for generic family coverage.'),('ten_representatives.json','The ten representative affine phase matrices.'),('family_rank_certificates.json','Exact generic variation equations and nonzero minors.'),('mixed_finite_squares.json; mixed_canonical.json','All completions at the thirty-two points and their canonical equivalences.'),('mixed_local_locus.json','The common rank-48 equations and rank-49 isolation certificates.'),('exceptional_local_*.json','Complete local regular branches at the twenty-seven missed catalogue points.'),('mixed_circle_analysis.json','Exact twelfth-root A locus and the mixed eighteenth-root example.'),('aalto_all_embeddings.json','Every one of the 5469 catalogue equivalence witnesses.'),('aalto_final_summary.json','Independent verification totals, by root order.')]
rows=[r'\begin{longtable}{p{0.43\textwidth}p{0.50\textwidth}}\toprule Output file&Role\\\midrule\endhead']
for name,desc in items:rows.append(r'\raggedright\small\texttt{'+name.replace('_',r'\_')+'} & '+desc+r'\\[3pt]')
rows.append(r'\bottomrule\end{longtable}');(P/'certificate_table.tex').write_text('\n'.join(rows))
print('Generated and checked all report matrices, cycle certificates, and tables.')
