
import matplotlib.pyplot as plt
import matplotlib as mtl


#mtl.rcParams['text.usetex']=True
#mtl.rcParams['text.latex.preamble']=[r'\usepackage{amsfonts}']
mtl.rcParams['text.usetex']=True
mtl.rcParams['font.family']="serif"

mtl.rcParams['lines.markersize']=4
mtl.rcParams['font.size']=12
mtl.rcParams['lines.linewidth']=2

def us2si(x,y):
	return (x/2.54,y/2.54)


import math
phi=(1+math.sqrt(5))/2

fig_center=us2si(phi*7,7)


 


styles=[ 'ro','gs','b^','cD' ]
colors=['r','g','b','c','k']

	
def  limitTicks(nx,ny):
	plt.gca().yaxis.set_major_locator(mtl.ticker.MaxNLocator(ny))
	plt.gca().xaxis.set_major_locator(mtl.ticker.MaxNLocator(nx))

