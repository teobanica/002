import json
from pathlib import Path
from fractions import Fraction as F
from torus import decode_line,line_param,encode
R=Path(__file__).resolve().parents[1]
L=json.loads((R/'outputs/class22_line_parameters.json').read_text());S=json.loads((R/'outputs/class22_line_scan.json').read_text())
MOD=int((R/'outputs/line_modulus.txt').read_text())
out=[];txt=['CLASS 22: EXPLICIT SQUARE TEMPLATES','All displayed entries are exact Laurent monomials in z, |z|=1.','Coordinates: a=s=t0, b=t+1/2=t1+1/2; x=exp(2*pi*i*a), y=exp(2*pi*i*b).','A circle is A*s+B*t=c mod 1. The stored parametrization is (s,t)=p+v*u, z=exp(2*pi*i*u).','']
def fmt(c,k):
    assert c in (0,MOD//2)
    term='1' if k==0 else ('z' if k==1 else f'z^({k})')
    return ('-' if c else '')+term
for l,s in zip(L,S):
    if not s['k5']:continue
    H=[[[int(F(*z[0])*MOD),z[1]] for z in r] for r in l['rows']]+[s['rows'][j] for j in s['k5']]
    p,v=line_param(decode_line(l['line']))
    obj=dict(branch=l['branch'],epsilon=1 if l['branch']==0 else -1,line=l['line'],parameter_origin=encode(p),parameter_direction=list(v),constant_modulus=MOD,entry_encoding=['root exponent modulo stated modulus','integer exponent of z'],matrix=H)
    out.append(obj)
    txt.append(f"branch={l['branch']} epsilon={obj['epsilon']} circle={l['line']} p={p} v={v}")
    for r in H:txt.append('  '+'  '.join(fmt(c,k) for c,k in r))
    txt.append('')
(R/'outputs/class22_square_templates.json').write_text(json.dumps(out,indent=2))
(R/'outputs/class22_square_templates.txt').write_text('\n'.join(txt)+'\n')
print('Exported',len(out),'exact square templates.')
