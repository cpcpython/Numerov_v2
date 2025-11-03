#======================================================================
# IMPORTANT PARAMETERS (General)
#=======================================================================

import pickle
import math

########################################################################
# Numerov Parameters
########################################################################
grid      = 12000        # N division from 0 to xp ie. grid number
global      gridshift         # An Important Parameter in Forward Numerov Method,

itern     = 100               # after this iteration program exit; if you want more levels increase it

ynm1      = 1e-16        # far left point value
yn        = 2e-16        # ynm1+1th point value "representably small number" ; should be higher than |machine precision|of f32
xp        = 12.0         # range of oscillator as [0...xp] in atomic unit
Egap      = 0.01         # an empirical one to fit WFN at energy levels
#global node
node       =    0
nodeMax    =    20
EsearchMax =    0.       # for plotting, y's limit
Amp        =    1.0
# for Display X Y Limt
windowx    = [0.0 ,xp+0.1]            # xmin-xmax
windowy    = [-0.2,0.20]            # ymin-ymax

global Xcomplete                 # for complete range in x

print("Step size:h = ",xp/grid)  # it was 2*xp/grid since potential is symmetric wrt x=0 , now x=[0,xp]

# Important:
# PlotData- to save E,R,unNormalized Wfn upto grid shifted point, and then padded with zeros for reusing
# Storesdata as: PlotData[j][2] # (gridshifted,padded) WfnUnnormalized(Energy) ; PlotData[i][1] = x ;  Energy_i = PlotData[i][0][0]
global PlotData
PlotData = []

##############################################################
#  Use Either Ground state ot Excited State Block
##############################################################

####################################################################################
# Ground state Parameters Block
####################################################################################
'''
Req0 =   1.41691
# for Ground state only
InEmin     =    -0.18        # Starting E Minimum. it will automatically updated
InEmax     =    -0.16        # Starting E Maximum, it will automatically updated
# for Ground state :  gridshift =b when (InEmin to  a) in [[a,b],...]; a=gshift[0][0] b=gshift[0][1]     
gridShift = 8500    # default gridShift; need to cut end point's higher values; for lower vib.quantum numbers we need to increase it
gshift = [[InEmin,8500],[-0.1, 7000], [-0.05, 6000],[-0.01,50]]

global fittedX # trying array stored values for x and Energy(x)
global fittedE
# Select Ground/Excited State fitted data [xfitg/yfitg] or [xfite/yfite] (for Excited state)
with open('xfitg.pkl', 'rb') as file:  # use xfite.pkl for Excited state data
    fittedX = pickle.load(file)
with open('yfitg.pkl', 'rb') as file:  # use yfite.pkl for Excited state data
    fittedE = pickle.load(file)
print(len(fittedX))
'''
##################################################################################
# Excited State Parameters Block
##################################################################################

Req0 = 2.00576

# for Excited state only
InEmin=-0.11                # Starting E Minimum. it will automatically updated
InEmax=-0.09                # Starting E Maximum, it will automatically updated

# for Excited state
gridShift = 8000
gshift = [[InEmin,8000],[-0.08,6000],[-0.02,4000],[-0.01,50]]


global fittedX # trying array stored values for x and Energy(x)
global fittedE
# Select Ground/Excited State fitted data [xfitg/yfitg] or [xfite/yfite] (for Excited state)
with open('xfite.pkl', 'rb') as file:  # use xfite.pkl for Excited state data
    fittedX = pickle.load(file)
with open('yfite.pkl', 'rb') as file:  # use yfite.pkl for Excited state data
    fittedE = pickle.load(file)
print(len(fittedX))



