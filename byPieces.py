from sympy.printing import latex

class byPieces:
	def __init__(self, seps, funs):
		self.seps = sorted(seps)
		self.funs = funs
	def where (self, pos):
		if pos < self.seps[0]:
			return 0
		mn, mx = 0, len(self.seps) -1
		mid= (mx+mn)//2
		while (mx-mn>1):
			mid= (mx+mn)//2
			if pos < self.seps[mid]:
				mx=mid
			else:
				mn=mid
		return mid+1

	def __str__(self):
		s= "\t | (-oo, {}) \t : {} \n".format(self.seps[0],self.funs[0])
		for (sup,inf,fun) in zip(self.seps[1:],self.seps[:-1], self.funs[1:]):
			s+= "\t | ({}, {}) \t : {} \n".format(inf, sup , fun)
		s+="\t | ({}, +oo ) \t : {} \n".format(self.seps[-1],self.funs[-1])
		return s


	def latex(self, variable):
		funs=self.funs
		seps=self.seps
		s= r"""\begin{{case}}
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
		funs = [self.funs[0] + self.funs[1] ]
		for (sup,inf) in zip(seps[1:], seps ):
			mid = (sup+inf) /2
			funs.append(self(mid) + other(mid) )
		funs.append( self.funs[-1] + other.funs[-1] )
		return byPieces(seps,funs) 		


