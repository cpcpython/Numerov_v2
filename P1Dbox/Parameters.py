#======================================================================
# IMPORTANT PARAMETERS (General)
#=======================================================================

import pickle
import math

########################################################################
# Numerov Parameters
########################################################################
grid      = 12000        # N division from 0 to xp ie. grid number
global      gridShift    # An Important Parameter in Forward Numerov Method, see Manual

itern     = 200          # after this iteration program exit; if you want more levels increase it

ynm1      = 1e-16        # far left point value
yn        = 2e-16        # ynm1+1th point value "representably small number" ; should be higher than |machine precision|of f32
xp        = 8.0          # range of oscillator as [0...xp] in atomic unit
Egap      = 0.01         # an empirical one to fit WFN at energy levels
#global node
node       =    0
EsearchMax =    10.      # for plotting, y's limit
Amp        =    10.0
# for Display X Y Limt
windowx    =    [0.0 ,xp]                   # xmin-xmax
windowy    =    [0.,EsearchMax]             # ymin-ymax

global Xcomplete                            # for complete range in x

print("Step size:h = ",xp/grid)             # it was 2*xp/grid since potential is symmetric wrt x=0 , now x=[0,xp]

# Important:
# PlotData- to save E,R,unNormalized Wfn upto grid shifted point, and then padded with zeros for reusing
# Storesdata as: PlotData[j][2] # (gridshifted,padded) WfnUnnormalized(Energy) ; PlotData[i][1] = x ;  Energy_i = PlotData[i][0][0]
global PlotData
PlotData = []


####################################################################################
# Ground state Parameters
####################################################################################
InEmin=0.0                # Starting E Minimum. it will automatically updated
InEmax=0.01                # Starting E Maximum, it will automatically updated

gridShift = 50

##################################################################################

# define Potential function Here Zero-shifted asymptot is used
def Potential(x):
# here V=0 inside the box.    
    return 0.0

#Finalalizing important Parameters
EsearchMax =     10.      # for plotting, y's limit
itern      =     200         # after this iteration program exit; if you want more levels increase it
node       =     0
nodeMax    =     20         # Limit the node of wfn to stop the run

