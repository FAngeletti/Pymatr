


import sympy as sym
from sympy.matrices import *
import matplotlib.pyplot as plt

sym.init_printing(use_unicode=True)


d=2
p=0.25
q=1-p
E = Matrix( [ [ 1, 1, 1], [0,1, 0], [0, 0, 1]  ]  )
Qs =[  Matrix( [ [ -5, 0,0 ], [0,5,0], [0,0,10] ] ) ]


d=E.shape[0]
A = ones((d,d))

import model as Mod


red= Mod.reduced(A,E, Qs )  

import random
import math
def Gs(i,j):
	return lambda : random.gauss( Qs[0][i,j], 1  )  
import Pymatr.synthesis as Syn

from utils import numerical
L= red.dEigen
nsyn=200
Gen = Syn.MatrixRngOpt(numerical(A),numerical(E/L), Gs, nsyn)

def average():
	s= sum(Gen())
	av= s /nsyn
#	print(" \n sum {},  average: {}\n".format(r, av) )
	return av

lln=red.lln()
import Pymatr.byPieces as Bp
Bp.plot(lln)


import Pymatr.histogram as H
nsample=20000
H.plot(nsample, average)
import matplotlib.pyplot as plt


#plt.xticks( [ Q1r[0].evalf() ], ["µ"] )
plt.show()	



