import random
import Algo.Analysis

def mesh(n,f):
	for i in range(n):
		for j in range(n):
			yield ((i,j),f(i,j))

# Dichotomic search of x in a sorted list L 
def search(L,x):
	min,max=0,len(L)-1
	while(max-min>1):
		middle=(max+min)//2
		if L[middle]>x:
			max=middle
		elif x==L[middle]:
			return middle
		else:
			min=middle
	if L[min]>x:
		return min
	else:
		return max
	
from functools import reduce
class FiniteRng:
	def __init__(self,g):
#		print(g)
		self.values=[]
		self.cps=[]
		def build(c, tup  ): 
			x,p = tup
			self.values.append(x)
			self.cps.append(c+p)
			return c+p
		reduce(build,g,0)
		total=self.cps[-1]
		self.cps= [ x/total for x in self.cps ]
#		print(self.cps)
		
	def __call__(self):
		u=random.uniform(0,1)
		return self.values[search(self.cps,u)]


class MatrixRng:
	def __init__(self,A,E,Gs):
		self.dim=A.shape[0]
		self.At=A.T
		self.L=lambda M : trace(self.At*M)
		self.E=E
		#self.R=E**(-1)
		self.Gs=Gs


	def __call__(self,n):
		cT=self.E**n
		m=mesh(self.dim, lambda i,j : cT[i,j]*self.At[i,j] )
		fG=FiniteRng(m)
		(source,sink)=fG()
		state=source
		def transitions(i,cT):
			r=[ (j,self.E[i,j]*cT[j,sink]) for j in range(self.dim) ]
			return r
		for i in range(n):
			#cT=cT*self.R
			cT=self.E**(n-i-1)
			G=FiniteRng(transitions(state,cT))
			nstate=G()
	#		print ("{} → {}" .format(state,nstate) )
			law=self.Gs(state,nstate)
			r= law()
			yield r
			state=nstate



class MatrixRngOpt:
	def __init__(self,A,E,Gs,n ):
		self.n=n
		self.dim=A.shape[0]
		self.At=A.T
		self.L=lambda M : trace(self.At*M)
		self.E=E
		#self.R=E**(-1)
		self.Gs=Gs
		self.TransFin=[]
		self.init=None
		self.__init__Trans()
		self.__init__Init()

	def __init__Trans(self):
		for sink in range(self.dim):
			TransInh=[]				
			cT = self.E**0		
			for k in range(self.n):
				Trans=[]
				for i in range(self.dim):
					r=[ (j,self.E[i,j]*cT[j,sink]) for j in range(self.dim) ]
#					print( "i={} → j | f={} at time {} \n E={} T={} \n L={} \n ".format(i, sink, k, self.E, cT, r) )
#					print(" Correction : {},  E {} \n".format(self.E,cT) )
#					print("\n transition : {} \n".format(r) )
					Trans.append(FiniteRng(r))
				cT=cT*self.E
				TransInh.append(Trans)
			self.TransFin.append(TransInh)
	def __init__Init(self):
		cT=self.E**self.n
		m=mesh(self.dim, lambda i,j : cT[i,j]*self.At[i,j] )
		self.init=FiniteRng(m)
		

	def __call__(self):
		n=self.n
		(source,sink)=self.init()
		state=source
		Trans= self.TransFin[sink]
		for i in range(n):
			#cT=cT*self.R
			cT=self.E**(n-i-1)
			G= Trans[n-1-i][state]
			nstate=G()
	#		print ("{} → {}" .format(state,nstate) )
			law=self.Gs(state,nstate)
			r= law()
			yield r
			state=nstate

