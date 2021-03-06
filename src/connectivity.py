''' Compute the weak/strong component decompostion of a matrix.
    Relabel function for obtaining the Perron-Frobenius normal form of such matrix
    Block extraction function  

'''

import Pymatr.printing as pr
from Pymatr.utils import *
def computeConn(E):
	''' Compute the connection matrix Conn[i,j]>0 iff there is path from i to j in the graph
associated to E'''
	d=E.shape[0]
	#Computing epsilon such that I - \epsilon E is diagonal dominant 
	sl=2*max( [ sum(E.row(i)) for i in range(0,d) ] )  #sl = 1/epsilon
	#Inverting the matrix I - εE rather than  computing 1  + εE + (εE)^2 .... 
	M= MId(d) - E/sl
	Conn= M.inv()
	return Conn

def shape(M):
	"""Generate a latex representation of the shape of a  matrix """
	def sign(x):
		if x>0:
			return "+"
		if x<0:
			return "-"
		else:
			return "0"
	return pr.mat(M,pr=sign)	

def preorder(Conn):
	''' Compute the preorder (i connected to j) on the graph G(E) using the connection matrix Conn '''
	def prec(i,j):
		return Conn[i,j]>0
	return prec






def findExtrema(objs,prec):
	''' Find the extremal elements of a sets of objs for the partial ordering prec'''
	# We start with a trial decomposition [ potential extremal candidate , elements/ sucessors of the extremal candidate ]
	times=[set({}),objs.copy() ]
	# while there is still potential extremal elements
	while (len(times[1])>0):
		# we observe a new elements
		present=times[1].pop()
		# We then eliminate the sucessors of `present from both the list of candidate and the non-sucessor list 
		for (k,time) in enumerate(times):
			removal=set({})
			for i in time:
				if prec(present,i):
					removal.add(i)
			times[k]=time.difference(removal)
		# We then add the present element to the list of candidate. 
		times[0].add(present)	
		# Since we have removed the sucessors of `present from times[1], the invariant times[0]\inter\times[1]=ø is preserved.   
	return times[0]
				

def sortPreorder(objs,prec):
	''' Sort the elements of `objs accordingly to the preorder `prec
	 We sort the elements by separating then in level U_0, ... U_n
	 We want to preserve the invariant that elements of U_k have a chain of exactly k sucessors.
	 '''
	Us=[] 
	# The set V_0 is the whole sets
	current = { frozenset(x) for x in objs}
	step=0
	while(len(current)>0 ):
#		print("Current set at step {}: {}  \n".format(step, current) )
		# We find the extremal elements of the current sets
		news= findExtrema(current,prec)
#		print("Extrema {} : {} \n".format(step, news))
		# Add them to the current level
		Us.append(news)
		# And delete them to the current sets 
		current= current.difference(news)
		step+=1
	# We then compute a arbitrary ordering compatible with the Us ordering
	blocks=[]
	for current in Us:
		for block in current:
			blocks.append(block)
#	print(blocks)
	return (blocks,Us)


def stronglyConnected(Conn):
	'''  Partition the connected component in sorted strongly connected component '''
	d=Conn.shape[0]
	# prec(i,l) : is there a path from i to j?
	prec=preorder(Conn)
	# scon(i,j)  :is there a reversible path from i to j?
	def scon(i,j):
		return prec(i,j) and prec(j,i)
	# list of blocks
	blocks=[]
	# We construct the list of block indice by indice
	for i in range(d):
		found=False
		# we search a strongly connected component (a block) to which i is strongly connected 
		for b in blocks:
			if scon(i,b[0]):
				found=True
				b.append(i)
				break
		# if we have not found a block where i belongs, we create a new block [i]
		if not(found):
			blocks.append([i])

	# We want now to sort the strongly connected component
	# prec induces a partial ordering  cprec of the strongly connected component
	def cprec(b,c):
		eb=next(iter(b))
		ec=next(iter(c))
		return prec(eb,ec)
#	print("blocks: {} ".format(blocks))
	# Sorting the blocks
	blocks, Us = sortPreorder(blocks, cprec)
	return blocks

def antecedentSets(Conn, blocks):
	'''  Return the partial ordering of sets in functions of the antecedent level'''
	# We use the preorder on Conn
	prec=preorder(Conn)

	# We want now to sort the strongly connected component
	# prec induces a partial ordering  cprec of the strongly connected component
	def cprec(b,c):
		eb=next(iter(b))
		ec=next(iter(c))
		return prec(eb,ec)
	# Sorting the blocks
	blocks, Us = sortPreorder(blocks, cprec)
	return Us



def permutationList(comps):
	''' Compute the permutation induced by the partition in strongly connected component,
	 i.e. compute  a list of indice from the labels which are a set of sets of indices '''
	a=[]
#	print(labels)
	for c in comps:
		for i in c:
				a.append(i)
#	print("Permutation :{} \n".format(a) )
	return a


def relabelling(E):
	''' Compute the relabelling of E in order to obtain a block diagonal matrix where each block is irreducible'''
	d=E.shape[0]
	# We compute first the connection matrix Conn
	Conn=computeConn(E)
	# We then compute the strongly connected components
	comps=stronglyConnected(Conn)
	return comps



def blockDims(comps):
	''' Compute the size of each blocks '''
	dims=[]
	for c in comps:
				dims.append(len(c))
	return dims				

def relabel(comps, E ):
	''' Relabel the matrix E in function of the list of irreducible component'''
	d=E.shape[0]
	En=zeros(d,d)
	permutation=permutationList(comps)
	for (i,i_) in zip(permutation,range(d) ):
		for(j,j_) in zip(permutation,range(d) ):
			En[i_,j_]=E[i,j]
	return En

 
def cloneBlocks( dims, E):
	''' Extract the irreducible blocks from E'''
	start=0
	blocks=[]
	for d in dims:
		b=E[ start:(start+d), start:(start+d) ]
		blocks.append( b )
		start=start+d
	return blocks

