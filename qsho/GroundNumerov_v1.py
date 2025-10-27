import numpy as np
import matplotlib.pyplot as plt
import math
import pickle


# importing generic parameters
from GeneralFunctions import *

#-------------------------------------------------------------------------------------
# Function to plot Potential Parabola and Classical Turn points etc.
#-------------------------------------------------------------------------------------
def PlotParabolaPlus(Energies,K,XAll,PsiAll):
    #---------------------------------------
    # Plot (1/2)Kx**2 Parabola
    #---------------------------------------
    y=[]
    xarray=[]
    y=[ (1/2)*K*x**2 for x in XAll]
    xarray=XAll
    #-----------------------------------------
    En=Energies[0]
    # insert Lines of Energy into Parabola
    fig, ax = plt.subplots()
    ax.plot(xarray, y, linewidth=2.0) # plot parabola
    ax.set_xlabel('Distance from equilibrium (BL) in atomic unit', fontweight ='bold')
    ax.set_ylabel('Energy/Wavefunction (arb. unit)', fontweight ='bold')
    for i in range(len(Energies)):
        print("Energy Level,n=",i," E = ", Energies[i])
        En=Energies[i]
        ClassTP= math.sqrt(2*En/K) # Classical Turning Point Marking
        ax.hlines(y=En, xmin=ClassTP, xmax=-ClassTP, linewidth=1, color='y')
        ClassTPplt1 = plt.Circle((  -ClassTP,En ), 0.003 ); ClassTPplt2 = plt.Circle((  ClassTP,En ), 0.003 )
        ax.add_artist( ClassTPplt1)
        ax.add_artist( ClassTPplt2)
        mxPsi=max(PsiAll[i]); mnPsi=min(PsiAll[i])
        Scale=mxPsi-mnPsi
        #Egap = 0.02 # Empirical, which adjust height of Psi in energy level
        EadjustedWfn=list(map(lambda x:En+Egap*(x/Scale),PsiAll[i])) # note x/Scale : Confined Psi in Scale
        ax.plot(xarray,EadjustedWfn)

    ax.plot(xarray,EadjustedWfn)
    plt.ylim(0,EsearchMax) 
    pickle.dump(fig, open('FigureWfn.fig.pickle', 'wb'))  
    plt.savefig('P1DWavefunction.png')
    plt.show()
    plt.close()
    # --------------------------------------- Probability Plot -----------------------------------
    # insert Lines of Energy
    fig1, ax1 = plt.subplots()
    # plot 1D Box
    ax1.plot(xarray, y, linewidth=2.0) # plot parabola


    ax1.set_xlabel('Distance from equilibrium (BL) in atomic unit', fontweight ='bold')
    ax1.set_ylabel('Energy/Probability (arb. unit)', fontweight ='bold')
#   In first we make Psi*Psi2
    PsiAll2  = np.square(PsiAll) #list(map(lambda x: x ** 2, PsiAll))
    for i in range(len(Energies)):
        En=Energies[i]
        
        ClassTP= math.sqrt(2*En/K) # Classical Turning Point Marking
        ax1.hlines(y=En, xmin=ClassTP, xmax=-ClassTP, linewidth=1, color='y')
        ClassTPplt1 = plt.Circle((  -ClassTP,En ), 0.003 ); ClassTPplt2 = plt.Circle((  ClassTP,En ), 0.003 )
        ax1.add_artist( ClassTPplt1)
        ax1.add_artist( ClassTPplt2)
        
        ax1.hlines(y=En, xmin=-xp, xmax=xp, linewidth=1, color='y')
        mxPsi=max(PsiAll2[i]); mnPsi=min(PsiAll2[i])
        Scale=mxPsi-mnPsi
                       #0.005 # Empirical Which Determines Psi Height
        EadjustedProb=list(map(lambda x:En+1.0*Egap*(x/Scale),PsiAll2[i])) # note x/Scale : Confined Psi^2 in Scale
        ax1.plot(xarray,EadjustedProb)

    ax1.plot(xarray,EadjustedProb)
    plt.ylim(0,EsearchMax)
    pickle.dump(fig1, open('FigureProb.fig.pickle', 'wb')) 
    plt.savefig('P1DProbability.png')
    plt.show()
    plt.close()
#-------------------------------------------------------------------------------------
# Main Starting Input of Energies
# It first scans Energy from InEmin to InEmax then try to go up
#-------------------------------------------------------------------------------------

def ScanHermitefunctions(InEmin,InEmax,itern):
    #InEmin=0.00 ; InEmax = 0.01
    #Converged Energies, stores here
    Econverged=[]
    PsiAll=[]
    XAll=[]
    for i in range (itern): # upto range_32 is checked [ n=16 ] i upto 75 is checked.
        print("\n\n////////////////////////////////////////  I/ITERN -  ",i,"/",itern, "    /////////////////////////////////////")
        print("\n\n Searching Eigenfunction in [InEmin,InEmax] : " , InEmin,InEmax)
        ## FIND PSI BETWEEN[INEMIN TO INEMAX] INTERVAL ; Main Program Call Begins ...
        ee,xx,yy = findHermiteFunctions(InEmin,InEmax)
        if(ee != 0):
            Econverged.append(ee)
            # Normalized Psi should be used from here onwards: yy = unNormalized Psi
            # little bit Normalization work; N = Sqrt[ Int[ f*f]dx] so,----------------------
            l = [x * x for x in yy] # l=Psi^2 here
            NormalizedC=math.sqrt(Simpson(l,xx[0],xx[grid-2*gridShift-1],len(l)))
            # note we eliminated far right/left Psi values for noise reduction

            yyNormed = [x/NormalizedC for x in yy]
            #--------------------------------------------------------------------------------
            # Choice - 1
            # PsiAll.append(yy)       # for UnNormalized Psi :
            # Choice - 2
            PsiAll.append(yyNormed) # if we want Normalized Psi
            #--------------------------------------------------------------------------------
            # Checking whether Normalization correct or not:
            Psi2 = [x*x for x in yyNormed] # Psi^2
            IntPsiNorm=Simpson(Psi2,xx[0],xx[grid-2*gridShift-1],grid-2*gridShift)
            print("Normalization check: ∫(Psi_Normalized)^2 dx ==1 ? If  n= \t",i,"Calculated Value: \t",IntPsiNorm,"")

            XAll=xx
        InEmin=InEmax
        InEmax=InEmax+0.01
    return Econverged,XAll,PsiAll


# **Acual Calculation Section**
# By giving a suitable numerical value of xp (-xp to xp consists the bond elongation/compression with respect to the equilibrium bond length, BL, See the Paramter section), and other important parameters, now one can start the calculations. Note that the Probaility and Wavefunction plot are 'scaled' to fit into the parabola, legibly.

# In[ ]:


K=0.33                           # force constant of H2 molecule
EsearchMax =    0.25             # for plotting, y's limit
InEmin     =    0.0              # Starting E Minimum. it will automatically updated
InEmax     =    0.01             # Starting E Maximum, it will automatically updated
itern      =    20               # after this iteration program exit; if you want more levels increase it


Econverged,XAll,PsiAll=ScanHermitefunctions(InEmin,InEmax,itern)
PlotParabolaPlus(Econverged,K,XAll,PsiAll)


# In[ ]:


# Replotting
#figxw = pickle.load(open('FigureWfn.fig.pickle', 'rb'))  # Show the Wavefunction Figure interactively
#figxp = pickle.load(open('FigureProb.fig.pickle', 'rb')) # Show the Probability Figure interactively

