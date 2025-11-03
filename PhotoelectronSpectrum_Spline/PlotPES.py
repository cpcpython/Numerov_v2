# Plotting PES for H2-H2+ PES That can compare to Experimental spectrum
# last added FC factor data for: <v=0|v'=0-15>

import numpy as np
import matplotlib.pyplot as plt
import json
import math
from math import pi
from operator import add


# importing R and E(R) into these arrays
r = np.genfromtxt('PES.dat', usecols=0 )
ey= np.genfromtxt('PES.dat', usecols=1 )

leftShift  = 2.5
rightShift = 2.5

plt.bar(r, ey,width=0.05)
plt.xlabel("eV")
plt.ylabel("FC factor (arb. u.)")
plt.title("Photoelectron Spectrum")


plt.gca().invert_xaxis()
plt.show()

# Fitting Lorenzian; Parameters are below--------------------------------
y0 = 0.    # Shift 
A  = 0.05  # Amp.
w  = 0.1   # Width
dimn=20000 # datapoints

xplotdim=[0] * dimn # dim. of x of final spectrum
yplotdim=[0] * dimn # dim. of y of final spectrum
xploti=r[0]
xplotf=r[len(r)-1]
delx = (xplotf-xploti)/dimn

Xall = [] # final plot data for PES x,f[x]
Yall = []
# for convoluted sum spectrum only
YallSum = [0]*dimn
XallSum = [0]*dimn
delxf = (xplotf-xploti)/(dimn)
for i in range(dimn):
    XallSum[i] = xploti+i*delxf

lup = 0
for j in range(len(r)):
    ypes=[0]*dimn
    xpes=[0]*dimn 
    # Right direction wrt. center
    #for j in range(len(r)): # j is Energy
    for i in range(int(dimn/2)): # half xrange
        xpes[lup+i] = r[j]-delx*i
        ypes[lup+i] = y0+ey[j]*(2*A*w)/(math.pi*(4*(r[j]-i*delx - r[j])**2)+w**2)


    # Left direction wrt. center
    #for j in range(len(r)): # j is Energy
    for i in range(int(dimn/2)): # half xrange
        xpes[lup+int(dimn/2)+i] = r[j]+delx*i
        ypes[lup+int(dimn/2)+i] = y0+ey[j]*(2*A*w)/(math.pi*(4*(r[j]+i*delx - r[j])**2)+w**2)
    Xall.extend(xpes); Yall.extend(ypes);
    lup=lup


plt.plot(Xall,Yall, linestyle = "", marker = "o",markersize=1, color='k' )
plt.xlabel("eV")
plt.ylabel("FC factor (arb. u.)")
plt.title("Photoelectron Spectrum - unconvoluted")
plt.gca().invert_xaxis()
plt.show()

# Final convoluted Spectra--------------------------------------------------


#xplotdim=[0] * dimn # dim. of x of final spectrum
#yplotdim=[0] * dimn # dim. of y of final spectrum
xploti=r[0]-leftShift
xplotf=r[len(r)-1]+rightShift
delx = (xplotf-xploti)/dimn

# for convoluted sum spectrum only
YallSum = [0]*dimn
XallSum = [0]*dimn
delxf = (xplotf-xploti)/(dimn)

for i in range(dimn):
    XallSum[i] = xploti+i*delxf

for j in range(len(r)):
    ypes=[0]*dimn
    xpes=[0]*dimn
    # Right direction wrt. center
    #for j in range(len(r)): # j is Energy
    i=0
    while (i < dimn ): # half xrange
        xpes[i] = xploti+delx*i
        ypes[i] = y0+ey[j]*(2*A*w)/(math.pi*(4*(r[j]-(xploti+i*delx) )**2)+w**2)
        i=i+1

    #Xall.extend(xpes); Yall.extend(ypes);
    YallSum = list( map(add, YallSum, ypes))


plt.plot(xpes,YallSum, linestyle = "", marker = "o",markersize=1, color='k' )
plt.gca().invert_xaxis()
plt.xlabel("eV")
plt.title("Photoelectron Spectrum (Lorentzian Convoluted)")
plt.ylabel("FC factor (arb. u.)")
plt.show()


print("------------------------- END PES PLOT --------------------------")

'''
15.47541206	0.105419
15.73940563	0.197862
15.98962538	0.214717
16.2260713	0.177541
16.4487434	0.124655
16.65764168	0.078651
16.85276614	0.046152
17.03411677	0.025758
17.20169358	0.013892
17.35549657	0.007324
17.49552574	0.003809
17.62178108	0.001967
17.7342626	0.001015
17.83297029	0.000526
17.91790416	0.000275
17.98906421	0.000145
'''
