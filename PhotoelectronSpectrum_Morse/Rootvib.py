#!/usr/bin/env python3

#from Parameters_Splinev1  import grid
#from Parameters_Splinev1  import gshift

gshift = [[-0.11,8000],[-0.08,6000],[-0.02,4000],[-0.01,50]]     # for excited st
#gshift  = [[-0.18,8500],[-0.1, 7000], [-0.05, 6000],[-0.01,50]]  # for ground  st
grid=12000

import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

# open  Wfn file
with open('2.00576PlotData.pkl', 'rb') as file:   
#with open('1.41691PlotData.pkl', 'rb') as file:
    data = pickle.load(file) # ground
#---------------------------------------------------------------------

def Simpson(func,a,b,n):
# Simple Simpson Integrations over func array 
    n = n 
    h = (b - a) / (n - 1)
    I_simp = (h/3) * (func[0] + 2*sum(func[:n-2:2]) + 4*sum(func[1:n-1:2]) + func[n-1])
    return(I_simp)

#---------------------------------------------------------------------
# Utility : It find total number of nodes (Roots) between an interval [a,b] of a normalized wavefunction array
#
# "The method is guaranteed (by Brent) to converge, so long as the function can be evaluated within the initial interval
#   known to contain a root." Numerical Recipes in FORTRAN 90,pg-353: And the below function is NR's direct translation of the
#   FUNCTION zbrent(func,x1,x2,tol)!
#

import math
import numpy as np

# Test
N=grid
xgrid    = [i for i in range(N)]  # [0 1 2 ...N-1] grid of xvalues

def sign(a,b):
	sgn=+1.0
	if (b > 0):
		sgn=+1.0
	else:
		sgn=-1.0
	return a*sgn

def brent(a,b,sinf):
	ITMAX=100
	EPS=3.3e-08
	tol=12.0/12000 #1e-05  # ~ x2-x1 gap
#	print(a,b)
	iter1=1
	fa=sinf[a]#math.sin(a)
	fb=sinf[b] #math.sin(b)
	if ((fa > 0.0 and fb > 0.0) or (fa < 0.0 and fb < 0.0)):
	#	print('root must be bracketed for zbrent')
		return -1 
	c=b
	fc=fb
	while iter1<=ITMAX:
		if ((fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0)):
			c=a
			fc=fa #math.sin(a) -> fa
			d=b-a
			e=d;
		if (abs(fc) < abs(fb)):
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa;
		tol1=2.0*EPS*abs(b)+0.5*tol
		xm=0.5*(c-b)
		if (abs(xm) <= tol1 or fb == 0.0):
			zbrent=b
			#print(" Root Found(as arrayindex)=",b,round(b),"f[root-1,root,root+1]=",sinf[round(b)-1],sinf[round(b)],sinf[round(b)+1])
			#print(" Real Root ~ ",xgrid[round(b)])
			#false  0.0 values check for Padded Zeros in Wfn
			if(sinf[round(b)-1] != sinf[round(b)] != sinf[round(b)+1] and abs(sinf[round(b)]) > EPS  ): # wfn tolerence wrt. EPS
				print(" Real Root found: index= ",round(b)," in Wfn[0-grid] =",sinf[round(b)])
				return round(b)
			else:
			#	print("No Root/Nodes return -1")
				return -1 # means it is False Root due to padding
		if (abs(e) >= tol1 and abs(fa) > abs(fb)):
			s=fb/fa
			if (a == c):
				p=2.0*xm*s
				q=1.0-s
			else:
				q=fa/fc
				r=fb/fc
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
				q=(q-1.0)*(r-1.0)*(s-1.0)
			if (p > 0.0):
				e=d
				d=p/q
			else:
				d=xm
				e=d
		else:
			d=xm
			e=d
		a=b
		fa=fb
		#b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		temp1=0
		if(abs(d) > tol1):
			temp1=sign(tol1,xm)
		else:
			temp1=d
		b=b+temp1
		fb=sinf[round(b)] # here only we update the function f[index]


# end while
	print('zbrent: exceeded maximum iterations')
#	brent=b
	print("Root Not found: ",b)
	return -1 

# For ground state Wfn. normalization checking ... ######################################################

gshifti=0
gno=0
for i in range (len(data)): # j = exc.state v; i = grd. state v'
    ## FIND PSI BETWEEN[INEMIN TO INEMAX] INTERVAL
    yy =data[i][2] # shifted Wfn, un normalized, zeropadded, with complete Xrange
    xx =data[i][1] # shifted x
    # Energi_i = (data[i][0][0]

    # Normalized Psi should be used from here onwards: yy = unNormalized Psi
    # little bit Normalization work; N = Sqrt[ Int[ f*f]dx] so,---------------------- Ground
    l = [x * x for x in yy] # l=Psi^2 here
    NormalizedC=math.sqrt(Simpson(l,xx[0],xx[grid-1],len(l)))
    # note we eliminated far right/left Psi values for noise reduction
    yyNormed = [x/NormalizedC for x in yy]

    # Checking whether Normalization correct or not:
    Psi2 = [x*x for x in yyNormed] # Psi^2
    IntPsiNorm=Simpson(Psi2,xx[0], xx[grid-1],(int)(grid))

    print("--------------------------------- Nodes Count ------------------------------------\n\n")
    print("Check: ∫(Psi_Normalized)^2 dx ==1 ? If  vg = ",i,"Calculated : \t",IntPsiNorm,"E_i = ",data[i][0][0])

	# Checking every possible nodes; i => Particular energy
    totNode= 0    # total nodes
    nodesx = []   # approx location of nodes
    brakt  = 25 #  brakt number of subdivision from 1-grid

#----------------------------------------------------------
# choosing gshiftcut w.r.t Ei
#----------------------------------------------------------
    gcut=0; gshifti=0
    E=data[i][0][0]
    while( gcut  == 0):
        if(E > gshift[len(gshift)-1][0]):  # gcut near E~0; it need to be get say, 500 in 12000 grid
            gcut=gshift[len(gshift)-1][1]*40
            break
        if(E > gshift[gshifti][0] and E < gshift[gshifti+1][0]):
            gcut = gshift[gshifti][1]
            break
        else:
           gshifti = gshifti+1

    gshiftcut=gcut
#------------------------------------------------------

    print("gshiftcut",gshiftcut, data[i][0][0] ,gshift[gshifti][0],gshifti) # select gshift for Ei
    for x in range(0, grid-gshiftcut,brakt):
        sinf=yyNormed
        if(x<=(grid-gshiftcut)): # once gshift is reached , root search can be stopped
           noden=brent(x,x+brakt,sinf)  # ...)-2  is due to ... sinf[round(b)+1])
#           print(x)
        else:
            break
        if(noden != -1):
            totNode=totNode+1
            nodesx.append(xx[noden]) 
    print("Nodes for, v=",i," is: ",totNode," Approximate Locations = ", nodesx)

    print ("-----------------------------------------------------------------------------------")

