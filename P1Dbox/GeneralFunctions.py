# Genereic Numerov Functions (independent of Potential data)
# version 1 (05-08-2025): gpkm
#
# out data: 
# LARGE_Emidsho_-<floatn.>.png        - figures of unnormalized wavefunctions
# P1DWavefunction.png P1DProbability.png        - potential and wfn.prob. 
# <Req0>PlotData.pkl                            - Stores [E,X,Gridshifted_Wavefunction(X)_unnormalized_zeroPadded]
# FigureWfn.fig.pickle/FigureProb.fig.pickle    - Stores Wfn/Prob figures in pickle format for later use (interactive support)


import math 
from   Parameters import *
import matplotlib.pyplot as plt
import pickle


#-------------------------------------------------------------
# main function to get value of wavefunction on the array/grid
#-------------------------------------------------------------
def callNumerov(E,xp,n,ynm1,yn):
    h = (xp)/n         # grid size from -xp to xp
    m = 1.0             # mass of electron in atomic unit
# starting points
    y_n_minus_1= ynm1 # first y
    y_n        = yn   # second y
    y          =[]    # third y, the calculated point
    xarray     =[]    # x values in f(x
    x=.0 # before, for SHO only x=-xp             # x is a starting point, say from the far left point

    g_n_minus_1    =  2*m*(E-Morse(x    ))
    g_n            =  2*m*(E-Morse(x+h  ))
    g_n_plus_1     =  2*m*(E-Morse(x+2*h))
    y_n_plus_1     = (2*y_n*(1 - (5*h**2/12)*g_n) - y_n_minus_1*(1 + (h**2/12)*g_n_minus_1)) / (1 + (h**2/12)*g_n_plus_1)
    y.append(y_n_plus_1) # Solution array, unnormalized
    xarray.append(x+2*h) # corresponding x values
    x=x+h
    y_n_minus_1=y_n ;          g_n_minus_1 = g_n
    y_n        =y_n_plus_1;    g_n = g_n_plus_1
#    print("ll",y_n_plus_1)
    Ni=3
    while  x < xp :     # Numerov method loop
        g_n_minus_1 ;
        g_n         ;
        g_n_plus_1     =  2*m*(E-Morse(x+2*h))
        y_n_plus_1     = (2*y_n*(1 - (5*h**2/12)*g_n) - y_n_minus_1*(1 + (h**2/12)*g_n_minus_1)) / (1 + (h**2/12)*g_n_plus_1)
        y.append(y_n_plus_1) # Solution array, unnormalized
        xarray.append(x+2*h) # corresponding x values
        x=x+h
        y_n_minus_1=y_n
        y_n        =y_n_plus_1
        g_n_minus_1 = g_n
        g_n = g_n_plus_1
        Ni=Ni+1
    return xarray,y



#---------------------------------------------------------------------
#  Simpson 1/3 rule for Normalizing Eigenfunctions
#---------------------------------------------------------------------
def Simpson(func,a,b,n):
# Simple Simpson Integrations over func array 
    n = n 
    h = (b - a) / (n - 1)
    I_simp = (h/3) * (func[0] + 2*sum(func[:n-2:2]) + 4*sum(func[1:n-1:2]) + func[n-1])
    return(I_simp)


#-----------------------------------------------------------------------
# INITIAL BRACKETING OF ENERGY [Emin < delE < Emax] where delE is the Energy and we find Wavefunction of this delE
# Outdata: 
#-----------------------------------------------------------------------
def findEigenFunctions(Emin, Emax):

    InEmin=Emin; InEmax=Emax # Just storing the original values
    delE=Emin; xo,yo=callNumerov(delE,xp,grid,ynm1,yn)
    Eendymin=yo[(int)(grid)-1]

    delE=Emax; xo,yo=callNumerov(delE,xp,grid,ynm1,yn)
# error catch at the end point rising/lowering
    Eendymax=yo[(int)(grid)-1]

    Emid=((Emin+Emax)/2.);delE=Emid ; xo,yo=callNumerov(delE,xp,grid,ynm1,yn)
    Eendymid=yo[(int)(grid)-1]

    # BISECTION TO FIND APPROX. LOCATION OF delE in which y(+xp) ~ 0
    bisecN=50
    print("================ <<  FIRST BISECTION BEGIN >> ================")
    for egy in range(1,bisecN,1):

        delEmin=Emin; xo,yo=callNumerov(delEmin,xp,grid,ynm1,yn); Eendymin=yo[(int)(grid)-1]
        delEmax=Emax; xo,yo=callNumerov(delEmax,xp,grid,ynm1,yn); Eendymax=yo[(int)(grid)-1]
        Emid=((Emin+Emax)/2.);
        delEmid=Emid ; xo,yo=callNumerov(delEmid,xp,grid,ynm1,yn); Eendymid=yo[(int)(grid)-1]
        if(egy==1):
            print("First bisection - Loop-N, delEmin, delEmid, delEmax::\t", egy,delEmin,delEmid,delEmax,Eendymin,Eendymid,Eendymax)

        # if it finds ideal biseaction region which contains a definte Root, it will break the loop over here
        if(Eendymid*Eendymin < 0 ):
            print("Exit in First bisection - Values of delEmin,delEmid,delEmax::\t", delEmin,delEmid,delEmax)
            print("Ideal bisective region found ...")
            break;
        if(Eendymid*Eendymax < 0 ): # ideal bracket bw. mid-mx
            print("Exit in First bisection - Values of delEmin,delEmid,delEmax::\t", delEmin,delEmid,delEmax)
            print("Ideal bisective region found ...")
            break;
# Bisection method - Two possibilities [1] monotonically decreasing but all positive or negative values for Emid,Emin;
# or [2] region with real root. For Possibility [1] Lowest decreasing part is found for root search, if f(x) are having same sign
        if(abs(Eendymid) < abs(Eendymin)):
            Emin=delEmid;

    print("================ <<< FIRST BISECTION END >>> ================")

    # ------------------------------------------------------------------------------------
    # Exiting if Emin and Emax didnt change at all
    # sometime it wont give sufficient result, so it can exited from the below for loop
    # egy != 1 need since sometime a single bisection loop find optimum bisection bracket
    if((Emin==InEmin) and (Emax==InEmax) and egy != 1):
        print("No solutions in this Energy Interval,[",Emin,",",Emax, "]. Exiting ...")
        zeroP = [0] * gridShift
        return 0,xo[0:grid],yo[0:grid-gridShift]+zeroP

    # -------------------------------------------------------------------------------------
    # important: If first section gives Emin~Emid~Emax we dont go furthur and should be returned Null
    tolSec=1e-12
    if(abs(Emid-Emin)<tolSec and abs(Emid-Emax)< tolSec and abs(Emin-Emax)<tolSec):
        print("No solutions in this Energy Interval where Emin~Emid~Emax: Exiting from the Second Bisection ...")
        zeroP = [0] * gridShift
        return 0,xo[0:grid],yo[0:grid-gridShift]+zeroP # exiting ...

    bs1=Emin;bs2=Emid;bs3=Emax
    print("*** *** *** <<< SECOND BISECTION BEGIN >>> *** *** *** ",Emin,Emid,Emax)
    secondbs=1; tolSec=0.0000001    # Second bisection for finding approximated Eigenfunction
    while(secondbs < 100):          # Hopefully 50 bisection is enough !

        delEmid=(delEmax+delEmin)/2
        if((Eendymin)*(Eendymid) < 0):
            Emax=delEmid; print("Root=================1 N Emin Emax:\t",secondbs, Emin,Emax)
        if ((Eendymax)*(Eendymid) <0 ):
            Emin=delEmid; print("Root=================2 N Emin Emax:\t",secondbs, Emin,Emax)

        delEmin=Emin; xo,yo=callNumerov(delEmin,xp,grid,ynm1,yn); Eendymin=yo[(int)(grid)-1]
        delEmax=Emax; xo,yo=callNumerov(delEmax,xp,grid,ynm1,yn); Eendymax=yo[(int)(grid)-1]
        Emid=((Emin+Emax)/2.);
        delEmid=Emid ; xo,yo=callNumerov(delEmid,xp,grid,ynm1,yn); Eendymid=yo[(int)(grid)-1]

        secondbs=secondbs+1
        if(abs(Emid-InEmax) <=2e-15):# sometime Emid tends to the limit of InEmax, that give error, so it should be avoided
            print("Emid-InEmax are very small, Exiting....")
            zeroP = [0] * gridShift
            return 0,xo[0:grid],yo[0:grid-gridShift]+zeroP

        if(abs(Emid-Emin)< 2e-15): # Crucial step ~ Machine Precision
            print("*** Break in Second Bisection ***")
            # Right-End Points contains noises so it should be "symmetrically" removed like xo[100:1500] yo[100:1500]
# if we want to save PNG files
            plt.plot(xo[0:grid-gridShift],yo[0:grid-gridShift], linestyle='--', marker='o', color='g')
            plt.savefig('LARGE_Emidsho_'+str(Emid)+'H2.png')
            # PlotData contains the Full gridso Padding of zeros applied only in yy (due to its deviated values at ends)
            zeroP = [0] * gridShift
            PlotData.append([[Emid],xo[0:grid],yo[0:grid-gridShift]+zeroP ])
            plt.close()
            print("*************  Convergence Achieved ************")

            # save this FULL grid based x,Wavefunction(un normalized) for FC factor info.
            with open("PlotData.pkl", "wb") as file:
                pickle.dump(PlotData, file)

            # Main Return
            global node
            node=node+1
            return Emid,xo[0:grid],yo[0:grid-gridShift]+zeroP #break

    # Main Results, generally it wont come below; added just for an 'in case ...'
    zeroP = [0] * gridShift
    return 0,xo[0:grid],yo[0:grid-gridShift]+zeroP 
