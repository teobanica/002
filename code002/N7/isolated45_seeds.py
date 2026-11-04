"""The 45 isolated three-row representatives, in the reading draft's notation.
Integer entries are exponents of exp(2*pi*i/840). No unknown-row restriction
is assumed by using this common modulus for the seed data.
"""
from fractions import Fraction as Q
from itertools import product
MOD=840

def matrix(label,kind,rows,description=''):
    rows=[[int((Q(v)%1)*MOD) for v in row] for row in rows]
    rows=[[(x-row[0])%MOD for x in row] for row in rows]
    return dict(label=label,kind=kind,description=description,modulus=MOD,rows=rows)

def mixed(p):
    x,y,z,t=p;w=Q(1,3);h=Q(1,2)
    return [[0]*7,[0,w,2*w,w+z,w+z+h,2*w+t,2*w+t+h],
            [0,w-x,2*w-y,w,w-x+h,2*w,2*w-y+h]]

def seeds():
    out=[]
    aa=[[0,2,4,6,1,3,5],[0,3,6,2,5,1,4],[0,2,5,1,6,4,3],[0,2,6,5,3,1,4]]
    for i,p in enumerate(aa,1):out.append(matrix('A'+str(i),'AAA',[[0]*7,[Q(j,7) for j in range(7)],[Q(j,7) for j in p]]))
    for sign in [1,-1]:
        h=Q(sign,4)
        out.append(matrix('B+' if sign==1 else 'B-','BBB',[[0]*7,[0,Q(1,5),Q(2,5),Q(3,5),Q(4,5),h,h+Q(1,2)],
                    [h+Q(1,2),Q(-1,5),Q(-2,5),Q(-3,5),Q(-4,5),h,0]]))
    cc=[
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
    for i,(u,v) in enumerate(cc,1):out.append(matrix('E'+str(i),'CCC',[[0]*7,[Q(j,24) for j in u],[Q(j,24) for j in v]]))
    ac=[(1,3,4,6),(1,5,2,6),(1,6,2,5),(1,6,4,3),(2,3,4,5),(2,5,1,6),
        (2,5,4,3),(2,6,1,5),(4,3,1,6),(4,3,2,5),(4,5,2,3),(4,6,1,3)]
    for i,p in enumerate(ac,1):out.append(matrix('AC'+str(i),'ACC',mixed([Q(j,7) for j in p]),str(p)))
    i=0
    for (j,k),sign,form in product([(1,2),(1,3),(2,1),(3,1)],[1,-1],[1,2]):
        i+=1;h=Q(sign,4)-Q(j+k,5)
        p=(h,Q(j,5),h+Q(1,2),Q(k,5)) if form==1 else (Q(j,5),h,Q(k,5),h+Q(1,2))
        out.append(matrix('BC'+str(i),'BCC',mixed(p),f'j={j}, k={k}, epsilon={sign}, form={form}'))
    assert len(out)==45
    return out
