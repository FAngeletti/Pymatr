
from sympy import latex

def mat(M, pr=latex, env="pmatrix", spacing=""):
	n=M.shape[0]
	m=M.shape[1]
	s=r"\begin{{{env}}}".format(env=env)+"\n"
	for i in range(n):
		s+=pr ( M[i, 0] ) 
		for j in range(1,m):
			s+= "& {}".format(pr(M[i,j]))
		if not(i+1==n):
			s+=r"\\"
			if not( spacing==""):
				s+="[{sp}]".format(sp=spacing)
			s+="\n"
	s+="\n"+ r"\end{{{env}}}".format(env=env)
	return s


def vec(v, pr=latex):
	s=r"\p{{{}".format(latex(v[0]))
	d=len(v)
	for x in v[1:-1]:
		s+=", {}".format( pr ( x ) )
	s+=r"}"
	return s


def set(pr, s):
	res=r"\Set{"
	it=iter(s)
	try: 
		res+=pr(next(it))
	except StopIteration:
		res+="}"
		return res
	for x in it:
		res+="," + pr(x)	
	res+=r"}"
	return res

def nset(n, s):
	if n==1:
		return set(latex,s)
	else:
		return set( lambda s : nset(n-1, s) , s)

