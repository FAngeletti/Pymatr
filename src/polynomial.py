
import numpy

class vec:
	def __init__(self,arr):
		self.coeffs=[ x for x in arr]
	def __str__(self):
		s="( {}".format (self.coeffs[0])
		for x in self.coeffs[1:]:
			s=s + ", {}".format(x)
		s=s+")"
		return s

	def __len__(self):
		return len(self.coeffs)

	def dim(self):
		return len(self.coeffs)

	def __add__(self,other):
		return vec([ x+ y for (x,y) in zip(self,other) ])

	def __sub__(self,other):
		return vec([ x - y for (x,y) in zip(self,other) ])

	def __rmul__(self,scalar):
		return vec ([ x*scalar for x in self ])

	def __getitem__(self,key):
		return self.coeffs[key]	

	def __setitem__(self,key,val):
		self.coeffs[key]=val	

	def dot(self,other):
		return sum((x*y for (x,y) in zip(self,other)))

	def __truediv__(self,scalar):
		return vec([x/scalar for x in self])
	def __iter__(self):
		return vecIter(self)
	
	def zero(d):
		return vec([0]*d)


class vecIter:
	def __init__(self, vec):
		self.pos=0
		self.src=vec
	def __iter__(self):
		return self
	def __next__(self):
		if ( self.pos >= len(self.src) ) : raise(StopIteration)
		x=self.src[self.pos]
		self.pos+=1
		return x		



class Polynomial:
	def __init__(self,arr):
		self.coeffs = [ x for x in arr]
		self.degree=len(self.coeffs)-1
	def __str__(self):
		s=""
		if self.degree >-1:
			s="Polynomial : {}".format(self.coeffs[0])
			for (d,c) in enumerate(self.coeffs[1:]):
				s= s + " + {} X^{}".format(c,d+1)
		return s

	

	def __add__(self,other):
		if type(other)!=Polynomial:
			return other + self 
		if ( other.degree > self.degree ):
			return (other + self)
		coeffs =  [ x for x in self.coeffs]
		for (i,c) in enumerate (other.coeffs):
			coeffs[i]+=c
		return Polynomial(coeffs)

	def __radd__(self,other):
		coeffs= [c for c in self.coeffs ]
		coeffs[0]+= other 
		return Polynomial(coeffs)

	def __rmul__(self, other ):
		return Polynomial([ c*other for c in self.coeffs])

	def __sub__ (self, other):
		return self  + (-1) * other

	def __rsub__(self,other):
		coeffs= [ -c for c in self.coeffs ]
		coeffs[0] += other 
		return Polynomial(coeffs)


	def __mul__(self,other):
		if type(other)!=Polynomial:
			return NotImplemented
			# convolution . To be eventually replaced with exact fft
		nc = self.degree + other.degree + 1
		coeffs=[0] *nc
		for (i,c1) in enumerate(self.coeffs):
			for (j,c2) in enumerate(other.coeffs):
				coeffs[i+j]+= c1*c2
		return Polynomial(coeffs)

	def __truediv__(self,scalar):
		return Polynomial([c/scalar for c in self.coeffs])


	def __call__(self,x):
		s=self.coeffs[-1]
		for c in self.coeffs[-2::-1]:
			s=s*x +c
		return s
	
				
c = Polynomial([1] )
X= Polynomial([0,1])
print (c)
print(X)

P1= X - c
print(P1)
P2 = P1 *P1
print(P2)
print ( P2(2) )

v=vec([5,7])
v2=numpy.array([5,7])

P3= v2*(X*X-c) /2
print(P3[0])

	
				  
