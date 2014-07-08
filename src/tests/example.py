
import sympy as sym
from sympy.matrices import *
import matplotlib.pyplot as plt

sym.init_printing(use_unicode=True)

"""
def randomM(d , alea= random.random):
	M=Matrix(d,d, [alea() for i in range(d*d) ])
	return M



E0 = Matrix( [   [ 0, 1, 1  , 0, 0  , 0   ,  0  , 0    ], 
                    [ 1, 0, 0  , 0, 0  , 0   ,  0  , 0    ],
                    [ 0, 0, 0.4, 0, 0.4, 0.1 , 0.01, 0    ],
                    [ 0, 0, 0  , 1, 0  , 0   , 0   , 0    ],
                    [ 0, 0, 0.4, 0, 0.4, 0   , 0   , 0.1  ], 
		    [ 0, 0, 0  , 0, 0  , 0.2 , 0.8 , 0    ],	
                    [ 0, 0, 0  , 0, 0  , 0.1 , 0   , 0.9  ],
		    [ 0, 0, 0  , 0, 0  , 0.3 , 0.3 , 0.4  ]  ] )

d=E0.shape[0]



Q1_= randomM(d)
Q2_= randomM(d) +np.multiply(Q1_,Q1_)
"""

# Random generation
import Pymatr.generation

(E0, Qs)  = Pymatr.generation.system(12, 20, sBmx=2)

print(E0)
d=E0.shape[0]
A0 = ones((d,d))

# Jordan-Frobenius decomposition
from Pymatr.connectivity import *

Conn=computeConn(E0)
comps=relabelling(E0)
dims = blockDims(labels)

[E,A,Q1,Q2]=[ relabel(comps, M) for M in [E0,A0,Qs[0],Qs[1] ] ]
Bs = cloneBlocks(dims, E )


# Eigentriples computation
from Pymatr.eigentriples import *
eigens = eigentriples(Bs)
dEigen= dominantEigenvalue(eigens)
dominants =findDominants(eigens)

# Reduced statistics computation
from Pymatr.reduction import *
shadow=shadowTransition(E,dims, dominants)
lshadow=shadowLimit(shadow)
Er= reducedE(dominants,dims,eigens,lshadow)
Ar=reducedA(A,dominants,dims, eigens)
[Q1r,Q2r] = [ reducedStatistic(M,dominants, Bs, dims, eigens) for M in [Q1,Q2] ]
sigma2s= Q2r - Q1r*Q1r



# Limits laws computations
from Pymatr.limits import *
paths=maxPaths(Er)
probP = probPaths(Ar,Er,paths)

# LLN limit distribution
llnLaw=llnLimit(Q1r,paths,probP)


#sigma2s=np.array([0.99,0.01,0])

# CLT computation
xtcl= np.arange(-2,2,0.025)
tclLaw=tclLimit(sigma2s,paths, probP )
#ptcl = tclLaw(xtcl) 


print(" \
matrix E : \n {}, \n \
connection matrix Conn :\n {} \n \
weakly / strongly connected component: \n {} \n \
block dimensions: \n {} \n  \
relabelled matrix E : \n {} \n  \
blocks B(i) : \n {} \n \
list of eigenvalue, left and right eigenvectors : \n {} \n \
dominants blocks : \n {} \n \
shadow transition matrix X : \n {} \n \
limit shadow transition matrix X : \n {} \n \
\n {} \n \
reduced E : \n {} \n \
reduced A : \n {} \n \
reduced Q(1) : \n {} \n \
reduced Q(2) : \n {} \n \
maximal paths : \n {} \n \
path probabilities : \n{}\n \
variances : {} \n \
Law of large number limit distribution: \n {} \n \
  ".format(E0,Conn,labels,dims,E,Bs,eigens,dominants,shadow,lshadow,"="*79,Er,Ar,Q1r,Q2r,paths,probP,sigma2s,llnLaw.latex("s")))

#plt.plot(xtcl,ptcl)
#plt.show()
