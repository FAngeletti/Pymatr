from generation import *
from connectivity import *
from sympy.matrices import *
import sympy 

d=4
g=skelFull(d, ((d-1)*d)//2 )

print(g)

E=zeros(d,d)

for (i,j) in g:
	E[i,j]= sympy.Rational(1)


relbl = relabelling(E)
perm=permutationList(relbl)

print(relbl)

E2 = relabel(perm,E)

print(E2)
