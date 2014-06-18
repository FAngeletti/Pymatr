

import sys
sys.path.append('../..')

import sympy as sym
from sympy.matrices import *
import matplotlib.pyplot as plt

sym.init_printing(use_unicode=True)


d=2
p=0.25
q=1-p
E = Matrix( [ [ 2 , 1 ] ,  [1, 8] ] )
Qs =[  Matrix( [ [ -5  ,0 ] ,  [ 0 , 5] ] )
, Matrix( [ [ 1 ,2 ] ,  [0.5 , 0.1] ] ) ]



d=E.shape[0]
A = ones((d,d))

import Algo.model as Mod


red= Mod.reduced(A,E, Qs )  

import random
import math
def Gs(i,j):
	return lambda : random.gauss( Qs[0][i,j], math.sqrt(Qs[1][i,j]) )  
import Algo.Synthesis as Syn

from Algo.utils import numerical
L= red.dEigen
nsyn=100
Gen = Syn.MatrixRngOpt(numerical(A),numerical(E/L), Gs, nsyn)

def average():
	s= sum(Gen())
	av= s /nsyn
#	print(" \n sum {},  average: {}\n".format(r, av) )
	return av

import Algo.Histogram as H
nsample=1000
H.plot(nsample, average)
import matplotlib.pyplot as plt
mu=red.Qs[0][0].evalf()
print(mu)
plt.plot( [ mu, mu  ], [0, 1] )

#plt.xticks( [ Q1r[0].evalf() ], ["µ"] )
plt.show()	



