
class reduced:
	def __init__(self,A0,E0,Qs):
		
		# Jordan-Frobenius decomposition
		import Algo.connectivity as Cn
		Conn=Cn.computeConn(E0)
		comps=Cn.relabelling(E0)
		dims =Cn.blockDims(comps)

		[E,A] =[ Cn.relabel(comps, M) for M in [E0,A0 ] ]
		Qs = [ Cn.relabel(comps,Q) for Q in Qs ] 
		Bs = Cn.cloneBlocks(dims, E )


		# Eigentriples computation
		import Algo.eigentriples as Ei
		eigens = Ei.eigentriples(Bs)
		dEigen= Ei.dominantEigenvalue(eigens)
		self.dEigen=dEigen
		dominants =Ei.findDominants(eigens)

		# Reduced statistics computation
		import Algo.reduction as Rd 
		shadow=Rd.shadowTransition(E,dims, dominants)
		lshadow=Rd.shadowLimit(shadow)
		self.E=Rd.reducedE(dominants,dims,eigens,lshadow)
		self.A=Rd.reducedA(A,dominants,dims, eigens, lshadow)
		self.Qs = [ Rd.reducedStatistic(M,dominants, Bs, dims, eigens) for M in Qs ]


		# Limits laws computations
		import Algo.limits as Lim
		self.paths=Lim.maxPaths(self.E)
		self.probP =Lim.probPaths(self.A,self.E,self.paths)

	def clt(self):
		from Algo.limits import tclLimit
		return tclLimit(self.Qs[1], self.paths, self.probP )

	def lln(self):
		from Algo.limits import llnLimit
		return llnLimit(self.Qs[0],self.paths, self.probP) 
	
		
