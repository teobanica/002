from pathlib import Path
import json
R=Path(__file__).resolve().parents[1];OUT=R.parent/'report'
assert OUT.is_dir()
reps=json.loads((R/'outputs/c_representatives.json').read_text())
assert json.loads((R/'checks/final_verification.json').read_text())['all_checks_passed']
lines=[r'With $\zeta=\e(1/30)$, the matrices are $H_{ij}=\zeta^{E_{ij}}$',
       r'for the following exponent arrays. We repeat $K$ here so that the complete',
       r'C-sector can be read in one place.']
for r,label in zip(reps,['K','L_{30}',r'L_{30}^{\mathsf T}']):
 lines += [r'\[',r'E('+label+r')=\begin{pmatrix}']
 lines += ['&'.join(map(str,row))+r'\\' for row in r['matrix_mod30']]
 lines += [r'\end{pmatrix}.',r'\]']
(OUT/'c_representatives.tex').write_text('\n'.join(lines)+'\n')
d=json.loads((R/'outputs/fixed_pairs.json').read_text())[0]
fs=[d['families'][i] for i in d['keep']]
lines=[r'The eighteen circles over $P(1)$ are written as',
 r'\[r_j(u)_k=\e(c_k/30+m_ku),\qquad u\in\R/\Z.\]',
 r'The binary string in the last column is $(m_1,\ldots,m_8)$.',
 r'The displayed numbering starts at one; the machine family index is one less.',
 r'{\small\setlength{\tabcolsep}{4pt}',
 r'\begin{longtable}{rllc}\toprule',
 r'$j$&Types&$(c_1,\ldots,c_8)$&Mask\\\midrule\endhead']
for i,f in enumerate(fs):
 if f['free']:
  const=','.join(str(e[0]) for e in f['row']);mask=''.join(str(e[1]) for e in f['row'])
  lines.append(f"{i+1}&{''.join(f['profiles'])}&$({const})$&\\texttt{{{mask}}}"+r'\\')
lines += [r'\bottomrule\end{longtable}}']
(OUT/'c_circle_table.tex').write_text('\n'.join(lines)+'\n')
ind=json.loads((R/'outputs/independent_fixed_squares.json').read_text())[0]
lines=[r'\begin{center}\small',r'\begin{tabular}{rrr@{\qquad}rrr}\toprule',
 r'Circle $i$&Circle $j$&$c$&Circle $i$&Circle $j$&$c$\\\midrule']
ls=ind['free_lines']
for a,b in zip(ls[:9],ls[9:]):lines.append('&'.join(map(str,[a[0]+1,a[1]+1,a[2],b[0]+1,b[1]+1,b[2]]))+r'\\')
lines += [r'\bottomrule\end{tabular}\end{center}']
(OUT/'c_free_lines.tex').write_text('\n'.join(lines)+'\n')
lines=[r'\begin{center}',r'\begin{tabular}{lrrrrr}\toprule',
 r'Matrix&A-pairs&C-pairs&AAA triples&AAC triples&CCC triples\\\midrule']
for r,label in zip(reps,['K','L_{30}',r'L_{30}^{\mathsf T}']):
 a=r['rows'];b=r['columns'];assert a['profiles']==b['profiles'] and a['triple_census']==b['triple_census']
 vals=[a['profiles'].get('A',0),a['profiles'].get('C',0)]+[a['triple_census'].get(k,0) for k in ['AAA','AAC','CCC']]
 lines.append('$'+label+'$&'+'&'.join(map(str,vals))+r'\\')
lines += [r'\bottomrule\end{tabular}\end{center}',r'The exact C-cliques in the displayed representative labelling are:']
for r,label in zip(reps,['K','L_{30}',r'L_{30}^{\mathsf T}']):
 def blocks(d):return r'\quad '.join(r'\{'+','.join(map(str,c))+r'\}' for c in r[d]['C_cliques'])
 lines += [r'\[',r'\begin{aligned}',r'\text{Rows of }'+label+r':\;&'+blocks('rows')+r',\\',
 r'\text{Columns of }'+label+r':\;&'+blocks('columns')+r'.',r'\end{aligned}',r'\]']
(OUT/'c_census_table.tex').write_text('\n'.join(lines)+'\n')
print('Generated exact C-sector matrix, circle, free-line and census tables.')
