from Pymatr.utils import ewp
from numpy import vectorize, sqrt, matrix, zeros

def momentFun(E,q):
	def f(p):
		v = vectorize(lambda p : p.moment(q) )
		m = v(p) 
		shape=E.shape
		res= ewp(E,m, zeros_fun=lambda k, l : zeros((k,l)) )   
		return res
	return f	

class model:
	def __init__(self,N,A,E,Ps):
		self.N = N
		self.A=A
		self.E=E
		self.Ps = Ps
		self.En= self.proj(E**N) 

	def proj(self, M):
		dims=M.shape
		acc=0
		for i in range(dims[0]):
			for j in range(dims[1]):
				acc+= (self.A[i,j])*M[i,j]
		return acc

	def projNorm(self,M):
		return self.proj(M)/self.En

	def expectation(self, fs, pos ):
		pos.append(self.N+1)
		prod = matrix(self.E ** (pos[0]-1), copy=True)
		for f,k_old,k in zip(fs,pos[:-1], pos[1:]):
			m = f(self.Ps)
			prod *= m * ( (self.E)** (k-k_old-1) )
		pos.pop()
		return self.projNorm(prod)




	def covariance(self, k, l ):
		if k>l:
			k ,l = l, k 
		elif k==l:
			return self.variance(k)
		fun = momentFun(self.E,1)
		moment2p= self.expectation( [fun,fun] , [k, l ] )
		av = self.moment(k,1) * self.moment(l,1)
		return ( moment2p - av ) 

	def correlation(self,k,l):
		var = sqrt( self.variance(k) * self.variance(l) )
		return self.covariance(k,l)/var
		

	def moment(self,k, q):
		return self.expectation( [momentFun(self.E,q)], [k] )

	def variance(self,k):
		return (self.moment(k,2) - self.moment(k,1)**2)
 



class reduced:
	def __init__(self,A0,E0,Qs):
		
		# Jordan-Frobenius decomposition
		import Pymatr.connectivity as Cn
		Conn=Cn.computeConn(E0)
		comps=Cn.relabelling(E0)
		dims =Cn.blockDims(comps)

		[E,A] =[ Cn.relabel(comps, M) for M in [E0,A0 ] ]
		Qs = [ Cn.relabel(comps,Q) for Q in Qs ] 
		Bs = Cn.cloneBlocks(dims, E )


		# Eigentriples computation
		import Pymatr.eigentriples as Ei
		eigens = Ei.eigentriples(Bs)
		dEigen= Ei.dominantEigenvalue(eigens)
		self.dEigen=dEigen
		dominants =Ei.findDominants(eigens)

		# Reduced statistics computation
		import Pymatr.reduction as Rd 
		shadow=Rd.shadowTransition(E,dims, dominants)
		lshadow=Rd.shadowLimit(shadow)
		self.E=Rd.reducedE(dominants,dims,eigens,lshadow)
		self.A=Rd.reducedA(A,dominants,dims, eigens, lshadow)
		self.Qs = [ Rd.reducedStatistic(M,dominants, Bs, dims, eigens) for M in Qs ]


		# Limits laws computations
		import Pymatr.limits as Lim
		self.paths=Lim.maxPaths(self.E)
		self.probP =Lim.probPaths(self.A,self.E,self.paths)

	def clt(self):
		from Pymatr.limits import tclLimit
		return tclLimit(self.Qs[1], self.paths, self.probP )

	def lln(self):
		from Pymatr.limits import llnLimit
		return llnLimit(self.Qs[0],self.paths, self.probP) 
	
		
