''' Compute the weak/strong component decompostion of a matrix.
    Relabel function for obtaining the Perron-Frobenius normal form of such matrix
    Block extraction function  

'''

import Algo.printing as pr
from Algo.utils import *
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
	times=[set({}),objs.copy() ]
	while (len(times[1])>0):
		present=times[1].pop()
		for (k,time) in enumerate(times):
			removal=set({})
			for i in time:
				if prec(present,i):
					removal.add(i)
			times[k]=time.difference(removal)				
		times[0].add(present)	
	return times[0]
				

def sortPreorder(objs,prec):
	levels=[]
	print()
	current = { frozenset(x) for x in objs}
	step=0
	while(len(current)>0 ):
#		print("Current set at step {}: {}  \n".format(step, current) )
		news= findExtrema(current,prec)
#		print("Extrema {} : {} \n".format(step, news))
		levels.append(news)
		current= current.difference(news)
		step+=1
	blocks=[]
	for current in levels:
		for block in current:
			blocks.append(block)
#	print(blocks)
	return blocks


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
	# sorting blocks
	blocks= sortPreorder(blocks, cprec)
	return blocks



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
	d=E.shape[0]
	En=zeros(d,d)
	permutation=permutationList(comps)
	for (i,i_) in zip(permutation,range(d) ):
		for(j,j_) in zip(permutation,range(d) ):
			En[i_,j_]=E[i,j]
	return En

 
def cloneBlocks( dims, E):
	start=0
	blocks=[]
	for d in dims:
		b=E[ start:(start+d), start:(start+d) ]
		blocks.append( b )
		start=start+d
	return blocks

