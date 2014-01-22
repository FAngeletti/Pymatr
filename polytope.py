import numpy as np
#from polynomial import (X,c,vec)
#from polynomial import Polynomial as P
import sympy as sym
from sympy.matrices import *
#import sympy.matrices as Mat
from Algo.byPieces import *



sym.init_printing(use_unicode=True)
X= sym.symbols('s')

def vec(arr):
	return Matrix(arr)


def canon(d):
	def base(i):
		v=vec([0]*d)
		v[i]=1
		return v
	return base

def zero(d):
	return vec([0]*d)



def choose(s):
	for x in s:
		return x

def fact(n):
	p=1
	for k in range(1,n+1):
		p*=k
	return p


class Hyperplane:
	def __init__(self, v, pos):
		norm= sym.sqrt(v.dot(v))
		self.normal=v/norm
		self.pos=pos/norm

	def relPos(self,v):
		r=self.normal.dot(v) - self.pos
		return r

	def immerge(self,v):
		nzds = [i for (i,val) in enumerate(self.normal) if val!=0 ]
		nzd= nzds[0]
		w=[]
		for (i,x) in enumerate(v):
			if  (i != nzd):
				w.append(x)
		return vec(w)
				

	def intersection(self, v, w):
		hv = self.relPos(v)
		hw= self.relPos(w)
		t= hw/(hw-hv)
		r= t*v+(1-t)*w
		return self.immerge(r)	

	def __str__(self):
		return "Hyperplane <{}, x> = {}\n".format(self.normal, self.pos)
	def __contains__(self,v):
		return( self.relPos(v) < 0 )



class LocalHyperplane:
	def __init__(self, v, posFun, posLoc):
		norm= sym.sqrt(v.dot(v))
		self.normal=v/norm
		self.posLoc=posLoc/norm
		self.posFun=posFun/norm

	def relPosFun(self,v):
		r=self.normal.dot(v) - self.posFun
		return r

	def relPosLoc(self,v):
		r=self.normal.dot(v) - self.posLoc
		return r


	def immerge(self,v):
		nzds = [i for (i,val) in enumerate(self.normal) if val!=0 ]
		nzd= nzds[0]
		w=[]
		for (i,x) in enumerate(v):
			if  (i != nzd):
				w.append(x)
		return vec(w)



	def intersection(self, v, w):
		hw= self.relPosFun(w)
		d= self.normal.dot(w-v)
		T= hw/d
		r= T*v+  (1-T)*w
		return self.immerge(r)	

	def __str__(self):
		return "Hyperplane <{}, x> = {}\n".format(self.normal, self.posFun)
	def __contains__(self,v):
		return( self.relPosLoc(v) < 0 )	

class Face:
	def __init__(self):
		self.childs=[]
		self.parents=[]
	def __str__(self):
		st= ["\t \t Inner faces : "]
		for s in self.childs:
			st.append("\t")
			st.append(str(s))
		st.append("\n \t \t Outer faces :" )
		for s in self.parents:
			st.append("\t")
			st.append(str(s))
		st.append("\n")
		return "".join(st)




class fullPolytope:
	def __init__(self,dim):
		self.vertices={}
		self.lattice=[]
		self.dim=dim
		self.VertexId=-1


	def genVertexId(self):
		self.VertexId+=1
		return (self.VertexId)
	def addVertex(self,v):
		key=self.genVertexId()
		self.vertices[key]=v
		return key

	def __str__(self):
		st=["Polytope with face lattice: \nvertices:\n  \t"]
		for (i,v) in self.vertices.items():
				st.append("{}:{} ".format(i,str(v)))
		st.append("\n lattice : \n")
		for (d,level) in enumerate(self.lattice):
			st.append("\t Level {}: \n".format(d))
			for (vertices, face ) in level.items():
				st.append("\t \t vertices: {} \n".format(vertices))
				st.append( str(face) )
				st.append("\n")
		return "".join(st)

	def copy(self):
		P=fullPolytope(self.dim)
		for v in self.vertices.values():
			P.addVertex(v)
		P.lattice= [dict(d) for d in self.lattice]
		return P
		


	def body(self):
		[(v,B)]=self.lattice[self.dim].items()
		return (v,B)

	def vec(self,vId):
		return self.vertices[vId]

	def volumeX(self,loc):
		def tessel(path, level, vertices, face ):
			v0=choose(vertices)
			if (len(vertices)==level+1):
				allV= vertices.difference({v0}).union(path)
				d=len(allV)
				verts=[ self.vec(v)-self.vec(v0) for v in allV ]
				m=Matrix( [v.T for v in verts] )
				svol=m.berkowitz_det()/ fact(d)
				sign = svol.subs(X,loc)
				sign=sign/abs(sign)
				return sign*svol
			else:
				sublatt=self.lattice[level-1]
				sfaces= [ (s,sublatt[s]) for s in face.childs if v0 not in s]
				npath=path.union({v0})
				vol=0
				for (s,f) in sfaces:
					vol+=tessel(npath, level -1,s, f)
				return vol
		(v,B)=self.body()
		return tessel(frozenset({}),self.dim, v,B) 
			
	def intersection( self, H):
		I=fullPolytope(self.dim-1)
		
		edges= self.lattice[1]
		collision={}
		for (s,f) in edges.items():
			vw= [ v for v  in s]
			v,w = self.vec(vw[0]),  self.vec(vw[1]) 
			if ( (v in H) and (w not in H ))  or ( ( w in H) and (v not in H) ) :
					n= H.intersection(v,w)
					key=I.addVertex(n)
					collision[s]= frozenset({key})
		nlattice={ s:Face() for s in collision.values() }
		I.lattice.append(nlattice)
		for level in self.lattice[2:]:
			sublattice={}
			ncollision={}
			for (s,f) in level.items():
				nf=Face()
				verts=frozenset({})
				for subfaces in f.childs:
					if subfaces in collision.keys():
						colSubF= collision[subfaces]
						nf.childs.append(colSubF)
						verts=verts.union(colSubF)
				if len(verts) > 0 :
					ncollision[s]=verts
					sublattice[verts]=nf
			I.lattice.append(sublattice)
			collision=ncollision
		I.repairParents()
		return I
					

	def repairParents(self):
		for (levelC,levelP) in zip(self.lattice[:-1], self.lattice[1:]):
			for (s,f) in levelP.items():
				for cs in f.childs:
					levelC[cs].parents.append(s)







def simplex(vertices):
	P=fullPolytope(len(vertices) - 1 )
	for v in vertices:
		P.addVertex(v)
	P.lattice=[]
	P.lattice.append({ frozenset([i]): Face() for i in P.vertices.keys()  })
	for ldim in range(1,P.dim+1):
		cLevel={}
		for vertex in P.vertices.keys():
			for (fvertices, face ) in P.lattice[ldim-1].items():
				if not(vertex in fvertices):
					nvertices=fvertices.union([vertex])
					face.parents.append(nvertices)
					if nvertices in cLevel:
						cLevel[nvertices].childs.append(fvertices)
					else:
						F=Face()
						F.childs=[fvertices]
						F.parents=[]
						cLevel[nvertices]=F
		P.lattice.append(cLevel)
	return P


def miror(s,t):
	s2=set({})
	for x in s:
		s2.add(x+t)
	return frozenset(s2)

def mirorFace(s,f,t):
	s2=miror(s,t)
	F2= Face ()
	F2.childs= [ miror(s,t) for s in f.childs  ]
	F2.parents= [ miror(s,t) for s in f.parents  ]
	return (s2,F2)

def lift(s,t):
	return s.union(miror(s,t)) 

def liftFace(s,f,t):
	s2=lift(s,t)
	f2=Face()
	f2.parents= [ lift(s,t) for s in f.parents ]
	f2.childs= [ lift(s,t) for s in f.childs ]
	f2.childs.append(s)
	f2.childs.append(miror(s,t))
	return(s2,f2)




def Point(z):
	P=fullPolytope(0)
	P.addVertex(z)
	level0= { frozenset({0}) : Face() } 
	P.lattice.append(level0)
	P.dim=0
	return P
	


def rCube(d ,base , dext):
	if d==0:
		z=zero(dext)
		return Point(z)
	else:

		Cr=rCube(d-1,base,dext)
		C=Cr.copy()
		vt=base(d-1)
		t=2**(d-1)
		for v in Cr.vertices.values():
			C.addVertex(vt+v)
		for (level,levelN) in  zip(Cr.lattice,C.lattice):
			for (s,f) in level.items():
				(s2,f2)= mirorFace(s,f,t)
				levelN[s2]=f2
		C.lattice.append({})
		for (levelD,levelU) in zip(Cr.lattice,C.lattice[1:]):
			for (s,f) in levelD.items():
				(s2,f2)=liftFace(s,f,t)
				levelU[s2]=f2
		C.dim=d
	return C





def Cube(d):
	return rCube(d,canon(d),d)

def SimplexRect(d):
	z=zero(d)
	base=canon(d)
	l= [base(i) for i in range(d) ]
	l.append(z)
	return simplex( l  )



def SimplexStd(d):
	base=canon(d+1)
	l= [base(i) for i in range(d+1) ]
	return simplex( l  )



def VolIntersectionSimplex( normal, offset ):
	d= len(normal)
#	print("Dimension: {}".format(d))
	base=canon(d)
	verts=[base(i) for i in range(0,d)]
	verts.append(zero(d))
	S=simplex(verts)
#	print("Normal: {}, vertices : {}".format(normal, verts ) )
	inters = sorted( set( [ normal.dot(v)-offset for v in verts ]) )
	locs = [ (v + w)/2 for (v,w) in zip(inters[1:],inters) ]
	def volumeLoc(loc):
		H = LocalHyperplane ( normal, X + offset, loc+offset)
		HnS= S.intersection(H)
#		print(HnS)
		V= sym.simplify( HnS.volumeX(loc) )
		return V
	vols = [sym.Rational(0)] 
	vols.extend( [volumeLoc(loc) for loc in locs])
	vols.append(0)
	bpP = byPieces (inters, vols )
	bpP.normalize()
	return bpP



