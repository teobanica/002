import json,csv
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
R=Path(__file__).resolve().parents[1];O=R/'outputs';P=R/'report'
pilot=json.load(open(O/'pilot.json'))['triples'];aff=json.load(open(O/'affine.json'))['triples']
def dist(d):return ',\\ '.join(f'{k}:{v}' for k,v in d.items()) or '--'
raw=[r'\begin{tabular}{lrrl}\toprule Profile & Tested & Feasible & Dimension: count\\\midrule']
for t,d in pilot.items():raw.append(f"{t}&{d['total']:,}&{d['feasible']:,}&${dist(d['dimensions'])}$"+r'\\')
raw.append(r'\midrule Total&243,873&95,950&\\\bottomrule\end{tabular}')
(P/'raw_table.tex').write_text('\n'.join(raw))
tab=[r'\begin{tabular}{lrrlrr}\toprule Profile & Supports & Feasible & Dimension: count & Lifts & Components\\\midrule']
for t,d in aff.items():tab.append(f"{t}&{d['support_classes']}&{d['feasible_classes']}&${dist(d['class_dimensions'])}$&{d['normalized_label_lift_sum_over_classes']}&{sum(c['affine_components'] for c in d['classes'])}"+r'\\')
tab.append(r'\midrule Total&815&206&&4,505&5,572\\\bottomrule\end{tabular}')
(P/'canonical_table.tex').write_text('\n'.join(tab))
def support(P):return '|'.join(''.join(str(j+1) for j in b) for b in P)
tab=[r'\begin{longtable}{rrrlll}\toprule Class&$d$&Branches&$Q$&$R$&Completion locus\\\midrule\endhead']
for c in aff['AAA']['classes']:
    cid=c['id'];status={11:r'\S\ref{sec:class11}',19:r'\S\ref{sec:class19}',21:r'\S\ref{sec:class21}',22:r'\S\ref{sec:class22}'}.get(cid,'All parameters')
    tab.append(f"{cid}&{c['dim']}&{c['affine_components']}&${support(c['Q'])}$&${support(c['R'])}$&{status}"+r'\\')
tab.append(r'\bottomrule\end{longtable}');(P/'aaa_registry.tex').write_text('\n'.join(tab))
g=defaultdict(list)
for a,b,_ in json.load(open(O/'BBB_isolated_points.json')):g[int(F(*a)*120)].append(int(F(*b)*120))
tab=[r'\begin{longtable}{rl@{\qquad}rl}\toprule $i$ & Values of $j$ & $i$ & Values of $j$\\\midrule\endhead']
items=list(g.items());n=(len(items)+1)//2
for row in range(n):
    cells=[]
    for ix in [row,row+n]:
        if ix<len(items):
            i,jj=items[ix];cells.extend([str(i),', '.join(map(str,jj))])
        else:cells.extend(['',''])
    tab.append('&'.join(cells)+r'\\')
tab.append(r'\bottomrule\end{longtable}');(P/'bbb_points.tex').write_text('\n'.join(tab))
K=json.load(open(O/'ccc_canonical.json'))[0]['canonical'];(P/'ccc_matrix.tex').write_text('\\begin{pmatrix}\n'+'\\\\\n'.join('&'.join(map(str,r)) for r in K)+'\n\\end{pmatrix}')
print('Generated five exact-data tables/matrices.')
