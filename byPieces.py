from sympy.printing import latex
import sympy as sp
X=sp.symbols('s')
class byPieces:
	def __init__(self, seps, funs):
		self.seps = sorted(seps)
		self.funs = funs
	def where (self, pos):
		if pos < self.seps[0]:
			return 0
		elif pos > self.seps[-1]:
			return len(self.seps)	
		mn, mx = 0, len(self.seps) -1
		mid= (mx+mn)//2
		while (mx-mn>1):
			mid= (mx+mn)//2
			if pos < self.seps[mid]:
				mx=mid
			else:
				mn=mid
		return mn+1

	def __str__(self):
		s= "\t | (-oo, {}) \t : {} \n".format(self.seps[0],self.funs[0])
		for (sup,inf,fun) in zip(self.seps[1:],self.seps[:-1], self.funs[1:]):
			s+= "\t | ({}, {}) \t : {} \n".format(inf, sup , fun)
		s+="\t | ({}, +oo ) \t : {} \n".format(self.seps[-1],self.funs[-1])
		return s


	def latex(self, variable):
		funs=self.funs
		seps=self.seps
		s= r"""\begin{{cases}}
 {} & {} \in (-\infty,{})\\
""".format(latex(funs[0]), variable, seps[0] )
		for (sup,inf,fun) in zip(seps[1:],seps[:-1], funs[1:]):
			s+=r"""{} & {} \in [{},{})\\
""".format(latex(fun),variable, inf, sup)
		s+= r""" {} & {} \in [{}, +\infty )
 \end{{cases}}""".format(latex(funs[-1]), variable, seps[-1] )
		return s    

	def __call__(self,x):
		return self.funs[self.where(x)]

	def __rmul__(self,scalar):
		funs = [scalar*fun for fun in self.funs]
		return byPieces ( self.seps, funs )
	def __add__(self, other ):
		seps =[ s for s in self.seps]
		seps.extend(other.seps)
		seps=sorted(set(seps))		
#		print( " {{ {} }} + {{ {} }} = {{ {} }} ".format(self.seps, other.seps, seps ))
		funs = [self.funs[0] + other.funs[0] ]
		for (sup,inf) in zip(seps[1:], seps ):
			mid = (sup+inf) /2
			s= self(mid) + other(mid)
#			print('sum {} at {}, left {}:{}, right {}:{} \n'.format(s, mid ,self(mid),self.where(mid), other(mid), other.where(mid) )  )
			funs.append( s )

		funs.append( self.funs[-1] + other.funs[-1] )
		return byPieces(seps,funs) 		

	def normalize(self):
		norm=0
		for (sup,inf,fun) in zip(self.seps[1:],self.seps[:-1], self.funs[1:]):			
			fi=sp.integrate(fun,X)
		#	print('fun {}, fi {}, Norm {}, sup {}, inf {} \n'.format(fun, fi,norm, sup, inf))
			norm+=fi.subs(X,sup)-fi.subs(X,inf)
			#print('Sum : {} \n'.format(norm))
		for (i,f) in enumerate(self.funs):
				self.funs[i]=f/norm
#				print('fun {} : {} \n\n'.format(i , self.funs[i] ) )
				


import matplotlib.pyplot as plt
import numpy as np
def plot(bpF, plot=plt.plot ,  n=20):
	npieces=len(bpF.funs)-2
	mn=min(bpF.seps)
	mx=max(bpF.seps)
	xs= np.zeros(n*npieces+4)
	ys=  [ 0 for x in xs]
	eps= (mx - mn) / n**2  
	xs[0]= mn - 1
	xs[1] = mn -eps 
	xs[-1]= mx + 1 
	xs[-2] = mx +eps
	offs=2
	for (f, inf,sup) in zip(bpF.funs[1:-1], bpF.seps[:-1], bpF.seps[1:] ) :
	#	print ('Inf {} - > Sup {} :  {}  \n'.format(inf,sup, f ) )
		step =(sup-inf)/(n-1)
		mesh = np.arange(inf.evalf(), (sup+0.5*step).evalf(), step.evalf()) 
	#	print(mesh)
		xs[offs : (offs +n) ] = mesh  
		ys[offs : (offs +n) ]= [ f.subs(X,x).evalf() for x in mesh ]
	#	print(xs)
	#	print(ys)
		offs+= n 
	plot(xs, ys)
	ym = max(ys)
	for sep in bpF.seps:
		plt.plot( [sep, sep], [-0.05, ym],  'r--')
	plt.xticks( bpF.seps,  bpF.seps ) 
	
	
	
