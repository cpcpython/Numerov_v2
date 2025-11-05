#!/usr/bin/env python

from Parameters import grid
from GeneralFunctions import Simpson

import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

# asymptot energy levels in ground and excites states
asymptg=-1.0
asympte=-0.5
arrayn = 11 # dflt array visualier size

# open ground/excited Wfn files
with open('1.41691PlotData.pkl', 'rb') as file:
# Load the object from the pickle file
    data = pickle.load(file) # ground
with open('2.00576PlotData.pkl','rb') as file:
    datae=pickle.load(file)   # excitedenfunctions
#---------------------------------------------------------------------

FCgrid = np.zeros(shape=(len(datae),len(data))) # FC visualize array

Plotdat = []                                    # for Plot.dat PES data
for i in range (len(data)): # j = exc.state v; i = grd. state v'
    for j in range (len(datae)):
        ## FIND PSI BETWEEN[INEMIN TO INEMAX] INTERVAL
        yy =data [i][2] # shifted WfnG
        yye=datae[j][2] # shifted Wfn
        xx =data [i][1] # shifted x
        # Energy_i = data[i][0][0]

        # Normalized Psi should be used from here onwards: yy = unNormalized Psi
        # little bit Normalization work; N = Sqrt[ Int[ f*f]dx] so,---------------------- Ground
        l = [x * x for x in yy] # l=Psi^2 here
        NormalizedC=math.sqrt(Simpson(l,xx[0],xx[grid-1],len(l)))
        # note we eliminated far right/left Psi values for noise reduction
        yyNormed = [x/NormalizedC for x in yy]

        # Normalized Psi should be used from here onwards: yy = unNormalized Psi
        # little bit Normalization work; N = Sqrt[ Int[ f*f]dx] so,---------------------- Excited
        l = [x * x for x in yye] # l=Psi^2 here
        NormalizedC=math.sqrt(Simpson(l,xx[0],xx[grid-1],len(l)))
        # note we eliminated far right/left Psi values for noise reduction
        yyeNormed = [x/NormalizedC for x in yye]

        # Checking whether Normalization correct or not:
        Psi2 = [x*x for x in yyNormed] # Psi^2
        IntPsiNorm=Simpson(Psi2,xx[0], xx[grid-1],(int)(grid))
#        print("Check: ∫(Psi_Normalized)^2 dx ==1 ? If  vg = ",i,"Calculated : \t",IntPsiNorm,"E_i = ",data[i][0][0])


    # for Excited state Normalization Checking ... #################################################

        # Checking whether Normalization correct or not:
        Psi2 = [x*x for x in yyeNormed] # Psi^2
        IntPsiNorm=Simpson(Psi2,xx[0], xx[grid-1],(int)(grid))
#        print("Check: ∫(Psi_Normalized)^2 dx ==1 ? If  ve = ",j,"Calculated : \t",IntPsiNorm,"E_i = ",datae[j][0][0])

        #Real FC Calculation is here
        Yge2 = list(map(lambda x, y: (x * y), yyNormed, yyeNormed))
        FCfactor=Simpson(Yge2,xx[0],xx[grid-1],(int)(grid))
        # now calculating FC by squaring it
        FCfactor = FCfactor**2
        # PE corrected vib level below
        print("FC FACTOR :: <ji> = <%2.2d|%2.2d>: %4.4f for DEL_E [j-i]eV = %4.4f"%(j,i,FCfactor,27.211*((asympte+datae[j][0][0])-(asymptg+data[i][0][0]))))
        
        FCgrid[j][i]=FCfactor
        if(i==0): # we need to plot only v=0 to other excited states
            Plotdat.append([27.211*((asympte+datae[j][0][0])-(asymptg+data[i][0][0])),FCfactor ] )

np.savetxt('PES.dat', Plotdat , delimiter=' ')

# show the FC matrix FC_ij
import matplotlib.pyplot as plt
plt.imshow(FCgrid, interpolation='none')
plt.show()
