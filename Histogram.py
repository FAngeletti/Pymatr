
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
	def __init__(self, l, lims,  nlim_d=5, nlim_b=5):
		self.els=l
		self.lims=lims
		(inf,sup) = lims
		self.n=len(l)
		self.nlimit_d=nlim_d
		self.nlim_b=nlim_b
		self.nlim_m = max(nlim_d/(sup-inf), nlim_b)
		
	def add(self,x):
		self.n+=1
		self.els.append(x)
		if self.n < self.nlim_m:
			return self
		else:
			return self.split()

	def split(self):
		s=sorted(self.els)
		half=len(s)//2
		left, right  = s[:half], s[half:]
		mid = ( left[-1] + right[0] ) /2
		shortest=min(mid-left[0], right[-1] - mid)
		nlimit= self.nlimit_d/shortest
		self.nlim_m=max(nlimit,self.nlimit_d, self.nlim_b)
		if self.n>self.nlim_m:
			return Node( left, right, mid, self.lims, 1.5*self.nlimit_d, self.nlim_b)
		else: 
			return self

	def density(self,n ):
		(inf,sup) = self.lims
		if inf == minfinity:
			inf = min(self.els)
		elif sup == infinity:
			sup = max(self.els)
		return [ (inf, sup , self.n / ( (sup-inf)*  n)) ]

	def width(self):
		inf, sup = lim
		return (sup-inf)


class Node:
	def __init__(self, left, right , mid, lims, nlim_d, nlim_b ) :
		self.lims = lims
		inf, sup = self.lims 
		self.left = Leaf(left, (inf,mid), nlim_d=nlim_d, nlim_b=nlim_b ) 
		self.right = Leaf(right, (mid,sup), nlim_d=nlim_d, nlim_b=nlim_b )
		self.mid = mid 
		self.n=len(left)+len(right)
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




 	
			
