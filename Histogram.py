
infinity= float('infinity')
minfinity= - infinity

class T:
	def __init__(self):		
		self.concrete= Leaf([], (minfinity,infinity) )
		self.n=0
	def add(self,x):
		self.n+=1
		self.concrete=self.concrete.add(x)
	def density(self):
		return self.concrete.density(self.n)




class Leaf:
	def __init__(self, l, lims,  nlim=100):
		self.els=l
		self.lims=lims
		self.n=len(l)
		self.nlimit=nlim
		
	def add(self,x):
		self.n+=1
		self.els.append(x)
		if self.n < self.nlimit:
			return self
		else:
			return self.split()

	def split(self):
		return Node(self.els, self.lims,  self.nlimit* (0.55 +  self.nlimit/200) )

	def density(self,n ):
		(inf,sup) = self.lims
		if inf == minfinity:
			inf = min(self.els)
		elif sup == infinity:
			sup = max(self.els)
		return [ (inf, sup , self.n / ( (sup-inf)*  n)) ]


class Node:
	def __init__(self, l , lims, nlim ) :
		s=sorted(l)
		half=len(s)//2
		left, right  = s[:half], s[half:]
		mid = ( left[-1] + right[0] ) /2
		self.lims = lims
		inf, sup = self.lims 
		self.left = Leaf(left, (inf,mid), nlim=nlim ) 
		self.right = Leaf(right, (mid,sup), nlim=nlim )
		self.mid = mid 
		self.n=len(l)
	def add(self,x):
		self.n+=1
		if x > self.mid:
			self.right= self.right.add(x)
		else :
			self.left= self.left.add(x)
		return self

	def density(self,n ):
		return self.left.density(n) + self.right.density(n)


import sys
def progressbar(f):
	width = 50
	adv = int( width * f) 
	todo =  width  - adv
	sys.stdout.write( "\r[{}>{}] ({}%) done".format('='*adv, ' '*todo, int(100*f)) ) 
	sys.stdout.flush()    

import matplotlib.pyplot as plt
def plot( n , rng, pl=plt.plot ):
	H=T()
	ndone= 0
	for x in range(n) :
		r= rng()
		H.add(r)		
		ndone+=1
		progressbar(ndone/n)
	xs=[]
	ds= [] 
	A=H.density()

	for (inf,sup, d)  in A:
		mid = (sup + inf) /2 
		#xs.extend([inf,inf,sup,sup])
	#	ds.extend([0,d,d,0])
		xs.append(mid)
		ds.append(d)

	pl(xs,ds)




 	
			
