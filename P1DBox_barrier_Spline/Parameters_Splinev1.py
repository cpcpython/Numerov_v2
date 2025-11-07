#======================================================================
# IMPORTANT PARAMETERS (General)
#=======================================================================

import pickle
import math

########################################################################
# Numerov Parameters
########################################################################
grid      = 12000        # N division from 0 to xp ie. grid number
global      gridshift    # An Important Parameter in Forward Numerov Method,

itern     = 50           # after this iteration program exit; if you want more levels increase it

ynm1      = 1e-16        # far left point value
yn        = 2e-16        # ynm1+1th point value "representably small number" ; should be higher than |machine precision|of f32
xp        = 8.0          # range of oscillator as [0...xp] in atomic unit
Egap      = 0.01         # an empirical one to fit WFN at energy levels
#global node
node       =    0
nodeMax    =    10
EsearchMax =    5.       # for plotting, y's limit
Amp        =    20.
# for Display X Y Limt
windowx    = [0.0 ,xp]           # xmin-xmax
windowy    = [0.,EsearchMax]

global Xcomplete                 # for complete range in x

print("Step size:h = ",xp/grid)  # it was 2*xp/grid since potential is symmetric wrt x=0 , now x=[0,xp]

# Important:
# PlotData- to save E,R,unNormalized Wfn upto grid shifted point, and then padded with zeros for reusing
# Storesdata as: PlotData[j][2] # (gridshifted,padded) WfnUnnormalized(Energy) ; PlotData[i][1] = x ;  Energy_i = PlotData[i][0][0]
global PlotData
PlotData = []

###############################################################
# trying array stored values for x and Energy(x)
# The below block is Applicable only for SPLINE based Numerov method
# Un comment this if you dont use SPLINE-Numerov
###############################################################

global fittedX
global fittedE
# Select Ground/Excited State fitted data [xfitg/yfitg] or [xfite/yfite] (for Excited state)
with open('xfitg.pkl', 'rb') as file:
    fittedX = pickle.load(file)
with open('yfitg.pkl', 'rb') as file:
    fittedE = pickle.load(file)


# for Starting Bisection 
InEmin=0.                 # Starting E Minimum. it will automatically updated
InEmax=0.01               # Starting E Maximum, it will automatically updated

# gridShift for all energy range; if needed increase it
gridShift = 50


