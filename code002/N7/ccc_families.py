"""Six affine families and eleven isolated CCC representatives.
All constants are exponents of exp(2*pi*i/24). Free phases range over R/Z.
"""
from classify_components import canonical, normalize_rows

def family(cu,cv,Bu,Bv):
    return normalize_rows(tuple(cu)+tuple(cv),tuple(Bu)+tuple(Bv),tuple(range(7)))

def templates():
    z=(0,0,0);a=(1,0,0);b=(0,1,0);c=(0,0,1)
    out=[('F1',family(
        [0,8,16,0,12,0,12],[0,16,8,0,0,12,12],
        [z,z,z,a,a,a,a],[z,z,z,b,c,b,c]))]
    z=(0,0);a=(1,0);b=(0,1);ab=(1,1)
    out.append(('F2',family(
        [12,8,16,0,0,0,12],[12,16,8,0,12,0,0],
        [z,z,z,z,z,a,a],[z,z,z,b,b,z,z])))
    out.append(('F3',family(
        [12,8,16,0,12,0,0],[12,16,8,0,0,0,12],
        [ab,z,z,z,a,a,ab],[ab,z,z,b,ab,z,b])))
    z=(0,);a=(1,);ai=(-1,)
    out.append(('F4',family(
        [0,8,16,8,20,0,12],[0,8,16,0,8,12,20],
        [z,z,z,a,a,a,a],[z,z,z,a,a,a,a])))
    out.append(('F5',family(
        [0,8,16,8,20,16,4],[0,20,0,8,16,8,12],
        [z,z,z,ai,ai,a,a],[z,a,a,z,z,a,a])))
    out.append(('F6',family(
        [0,8,16,0,12,0,12],[0,4,8,8,20,16,16],
        [z,z,z,a,a,ai,ai],[z,ai,a,z,a,ai,z])))
    pairs=[
    ([0,8,16,1,13,2,14],[0,2,11,8,23,14,16]),
    ([0,8,16,1,13,2,14],[0,2,23,11,8,14,16]),
    ([0,8,16,1,13,2,14],[0,6,5,17,14,12,22]),
    ([0,8,16,1,13,2,14],[0,6,17,14,5,12,22]),
    ([0,8,16,1,13,10,22],[0,10,7,19,16,22,8]),
    ([0,8,16,1,13,10,22],[0,10,19,16,7,22,8]),
    ([0,8,16,1,13,10,22],[0,22,9,6,21,12,14]),
    ([0,8,16,1,13,10,22],[0,22,21,9,6,12,14]),
    ([0,8,16,2,14,6,18],[0,4,2,10,16,18,12]),
    ([0,8,16,2,14,6,18],[0,16,2,4,10,18,12]),
    ([0,8,16,2,14,6,18],[0,16,14,22,4,12,6])]
    out.extend((f'E{k}',canonical(u+v,[()]*14)) for k,(u,v) in enumerate(pairs,1))
    return out

