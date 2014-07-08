
infinity= float('infinity')
minfinity= - infinity

class H:
	def __init__(self,inf, sup, nbins):
		self.nbins=nbins
		self.inf = inf
		self.sup = sup
		self.nbins=nbins
		self.width=sup-inf
		self.npoints=0
		self.array = [0] * nbins

	def localize(self,x):
		return int( self.nbins * (x-self.inf)/self.width )

	def add(self,x):
		self.npoints+=1
		if x  >= self.inf and  x <= self.sup :
			self.array[ self.localize(x) ] += 1

	def density(self):
		step = self.width/self.nbins
		xs  = [self.inf + k*step for k in range(self.nbins+1) ]
		xs.append(self.sup)
		ds = [ x/(self.npoints*step) for x in self.array]
		return zip( xs[:-1], xs[1:], ds ) 


import sys
def progressbar(f):
	width = 50
	adv = int( width * f) 
	todo =  width  - adv
	sys.stdout.write( "\r[{}>{}] ({}%) done".format('='*adv, ' '*todo, int(100*f)) ) 
	sys.stdout.flush()    

import matplotlib.pyplot as plt
def plot( n , inf, sup, nbins, rng, pl=plt.plot ):
	h=H(inf,sup,nbins)
	ndone= 0
	for x in range(n) :
		r= rng()
		h.add(r)		
		ndone+=1
		progressbar(ndone/n)
	xs=[]
	ds= [] 
	A=h.density()

	for (inf,sup, d)  in A:
		mid = (sup + inf) /2 
		#xs.extend([inf,inf,sup,sup])
	#	ds.extend([0,d,d,0])
		xs.append(mid)
		ds.append(d)

	pl(xs,ds)




 	
			
