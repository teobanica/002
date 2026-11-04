"""Literal six families from page 8 of the consolidated reading draft.
Entries are exp(2*pi*i*c/6) times monomials in the original seed parameters.
The all-ones row is prepended. Seed row phases are not changed.
"""
def family(name):
    d={'F1':3,'F2':2,'F3':2,'F4':1,'F5':1,'F6':1}[name]
    z=(0,)*d
    a=(1,)+(0,)*(d-1)
    if d>=2:b=(0,1)+(0,)*(d-2);ab=tuple(x+y for x,y in zip(a,b))
    if d==3:c=(0,0,1)
    if name=='F1':
        E=[[0,2,4,0,3,0,3],[0,4,2,0,0,3,3]]
        B=[[z,z,z,a,a,a,a],[z,z,z,b,c,b,c]]
    elif name=='F2':
        E=[[3,2,4,0,0,0,3],[3,4,2,0,3,0,0]]
        B=[[z,z,z,z,z,a,a],[z,z,z,b,b,z,z]]
    elif name=='F3':
        E=[[3,2,4,0,3,0,0],[3,4,2,0,0,0,3]]
        B=[[ab,z,z,z,a,a,ab],[ab,z,z,b,ab,z,b]]
    elif name=='F4':
        E=[[0,2,4,2,5,0,3],[0,2,4,0,2,3,5]]
        B=[[z,z,z,a,a,a,a],[z,z,z,a,a,a,a]]
    elif name=='F5':
        ai=(-1,);E=[[0,2,4,2,5,4,1],[0,5,0,2,4,2,3]]
        B=[[z,z,z,ai,ai,a,a],[z,a,a,z,z,a,a]]
    elif name=='F6':
        ai=(-1,);E=[[0,2,4,0,3,0,3],[0,1,2,2,5,4,4]]
        B=[[z,z,z,a,a,ai,ai],[z,ai,a,z,a,ai,z]]
    return dict(name=name,parameters=d,constants=[[0]*7]+E,basis=[[z]*7]+B)
