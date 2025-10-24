#!/usr/bin/env python
# use: python GroundNumerov_v1.py  2>&1 > ground.log
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

#=======================================================================
# IMPORTANT PARAMETERS 
#=======================================================================

# importing generic parameters
from GeneralFunctions import *

# calculate X range for Potential Function
Xcomplete = []
x=0
while x < xp:
        Xcomplete.append(x+(xp/grid)) # corresponding x values
        x=x+(xp/grid)

#-------------------------------------------------------------------------------------
# Function to plot Potential etc 
#-------------------------------------------------------------------------------------
def PlotParabolaPlus(Energies,XAll,PsiAll):
    #---------------------------------------
    # Plot Eenergy Levels, Wavefunction etc. w.r.t Energy 
    #---------------------------------------
    y=[]
    xarray=[]
    y=[ Potential(x) for x in Xcomplete]
    xarray=Xcomplete
    
    #-----------------------------------------
    En=Energies[0]

    # --------------------------------------- Wavefunction Plot -----------------------------------
    # insert Lines of Energy
    fig, ax = plt.subplots()
    ax.plot(xarray, y, linewidth=2.0) # plot parabola

    plt.axvline(x = 0, color = 'b', label = 'axvline - full height', linewidth=8)
    plt.axvline(x = xp, color = 'b', label = 'axvline - full height', linewidth=8)

    ax.set_xlabel('length in atomic unit', fontweight ='bold')
    ax.set_ylabel('Energy/Wfn (arb. unit)', fontweight ='bold')
    for i in range(len(Energies)):
        print("Energy Level,n=",i+1," E = ", Energies[i])
        En=Energies[i]
        ax.hlines(y=En, xmin=0, xmax=xp, linewidth=1, color='y') #0 xp
        mxPsi=max(PsiAll[i]); mnPsi=min(PsiAll[i])
        Scale=mxPsi-mnPsi
        #Scale=1
        EadjustedWfn=list(map(lambda x:En+Amp*Egap*(x/Scale),PsiAll[i]) ) # note x/Scale : Confined Psi in Scale
        ax.plot(Xcomplete[0:len(EadjustedWfn)] ,EadjustedWfn)
# Note it is for Zoomed graph
        ax = plt.gca()
        ax.set_xlim(windowx)   #[xmin, xmax])
        ax.set_ylim(windowy)   #[ymin, ymax])

    pickle.dump(fig, open('FigureWfn.fig.pickle', 'wb'))              # saving interactive plot - wfn
    plt.savefig('P1DWavefunction.png')
    plt.show()
    plt.close()

    # --------------------------------------- Probability Plot -----------------------------------
    # insert Lines of Energy
    fig1, ax1 = plt.subplots()
    ax1.plot(xarray , y, linewidth=2.0) # plot parabola
    plt.axvline(x = 0, color = 'b', label = 'axvline - full height', linewidth=8)
    plt.axvline(x = xp, color = 'b', label = 'axvline - full height', linewidth=8)

    ax1.set_xlabel('length in atomic unit', fontweight ='bold')
    ax1.set_ylabel('Energy/Probability (arb. unit)', fontweight ='bold')
#   In first we make Psi*Psi2
    PsiAll2  = [[item**2 for item in sublist] for sublist in PsiAll]
#    PsiAll2  = np.square(PsiAll) #list(map(lambda x: x ** 2, PsiAll)) #
    for i in range(len(Energies)):
        En=Energies[i]
        ax1.hlines(y=En, xmin=0, xmax=xp, linewidth=1, color='y')
        mxPsi=max(PsiAll2[i]); mnPsi=min(PsiAll2[i])
        Scale=mxPsi-mnPsi
#        Scale=1

        EadjustedProb=list(map(lambda x:En+Amp*Egap*(x/Scale),PsiAll2[i])) # note x/Scale : Confined Psi^2 in Scale
        ax1.plot(Xcomplete[0:len(EadjustedProb)],EadjustedProb)
# Note it is for Zoomed graph
        ax = plt.gca()
        ax.set_xlim(windowx)   #[xmin, xmax])
        ax.set_ylim(windowy)   #[ymin, ymax])

    pickle.dump(fig1, open('FigureProb.fig.pickle', 'wb'))
    plt.savefig('P1DProbability.png')
    plt.show()
    plt.close()

#-------------------------------------------------------------------------------------
# Main Starting Input of Energies
# It first scans Energy from InEmin to InEmax then try to go up
#-------------------------------------------------------------------------------------

def ScanPsifunctions(InEmin,InEmax,itern,nodeMax,node):
    #Converged Energies, stores here
    Econverged=[]
    PsiAll=[]
    XAll=[]
    delE=0.05
    for i in range (itern): # upto range_32 is checked [ n=16 ] i upto 75 is checked.
        print("\n\n////////////////////////////////////////  I/ITERN -  ",i,"/",itern, "    /////////////////////////////////////")
        print("\n\n Searching Eigenfunction in [InEmin,InEmax] : " , InEmin,InEmax)
        ## FIND PSI BETWEEN[INEMIN TO INEMAX] INTERVAL ; Main Program Call Begins ...
        ## IMPORTANT If the ee=0 means it is not the Root/Eigenvalue So node is NOT counted
        ee,xx,yy = findEigenFunctions(InEmin,InEmax)

        if(ee != 0):
#            Econverged.append(ee)
            # Normalized Psi should be used from here onwards: yy = unNormalized Psi
            # little bit Normalization work; N = Sqrt[ Int[ f*f]dx] so,----------------------
            l = [x * x for x in yy] # l=Psi^2 here
            NormalizedC=math.sqrt(Simpson(l,xx[0],xx[grid-1],len(l)))
            # note we eliminated far right/left Psi values for noise reduction

            yyNormed = [x/NormalizedC for x in yy]
            #--------------------------------------------------------------------------------
            # Choice - 1
            # PsiAll.append(yy)       # for UnNormalized Psi :
            # Choice - 2
            # check-2 ; End chek and monotonically decrease check at end  len(xx): 1e-14 is crucial limit
            #--------------------------------------------------------------------------------
            # Checking whether Normalization correct or not:
            Psi2 = [x*x for x in yyNormed] # Psi^2
            IntPsiNorm=Simpson(Psi2,xx[0], xx[grid-1],(int)(grid-1))
            print("Normalization check: ∫(Psi_Normalized)^2 dx ==1 ? If  itern= \t",i,"Calculated Value: \t",IntPsiNorm,"")
            PsiAll.append(yyNormed) # if we want Normalized Psi
            Econverged.append(ee)
            node=node+1
            XAll=xx
        InEmin=InEmax

        # Tolerence between E levels
        delE = 0.05+delE/10 #if InEmax < -0.025 else 0.001
        InEmax=InEmax+delE  # for Rydberg levels this 0.005 should be decreased

        if(InEmax >EsearchMax or node>nodeMax  ): # means Node > 10 wfn dont need to be calculated
            print(node,nodeMax)
            break  # to stop loop

    return Econverged,XAll,PsiAll


# ********* Acual Calculation Section **************************************************************
# By giving a suitable numerical value of xp ( xp consists the bond elongation) and other important parameters, 
# now one can start the calculations. Note that the Probaility and Wavefunction plot are 'scaled' to fit into the parabola, legibly.
# ***************************************************************************************************

Econverged,XAll,PsiAll=ScanPsifunctions(InEmin,InEmax,itern,nodeMax,node)
print(Econverged)
PlotParabolaPlus(Econverged,XAll,PsiAll)

# Replotting
#figxw = pickle.load(open('FigureWfn.fig.pickle', 'rb'))  # Show the Wavefunction Figure interactively
#plt.show()
#figxp = pickle.load(open('FigureProb.fig.pickle', 'rb')) # Show the Probability Figure interactively
#plt.show()
