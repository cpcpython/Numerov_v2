#======================================================================
# IMPORTANT PARAMETERS (General)
#=======================================================================

import pickle
import math

########################################################################
# Numerov Parameters
########################################################################
grid      = 12000        # N division from 0 to xp ie. grid number
global      gridShift         # An Important Parameter in Forward Numerov Method,

itern     = 10               # after this iteration program exit; if you want more levels increase it

ynm1      = 1e-16        # far left point value
yn        = 2e-16        # ynm1+1th point value "representably small number" ; should be higher than |machine precision|of f32
xp        = 12.0         # range of oscillator as [0...xp] in atomic unit
Egap      = 0.01         # an empirical one to fit WFN at energy levels
#global node
node       =    0
EsearchMax =    0.       # for plotting, y's limit
Amp        =    1.0
# for Display X Y Limt
windowx    = [0.0 ,12.1]            # xmin-xmax
windowy    = [-0.2,0.10]            # ymin-ymax

global Xcomplete                 # for complete range in x

print("Step size:h = ",xp/grid)  # it was 2*xp/grid since potential is symmetric wrt x=0 , now x=[0,xp]

# Important:
# PlotData- to save E,R,unNormalized Wfn upto grid shifted point, and then padded with zeros for reusing
# Storesdata as: PlotData[j][2] # (gridshifted,padded) WfnUnnormalized(Energy) ; PlotData[i][1] = x ;  Energy_i = PlotData[i][0][0]
global PlotData
PlotData = []

###################################################################################
# IMPORTANT : PLEASE USE EITHER GROUND STATE PARAMETERS OR EXCITED STATE PARAMETERS
#             COMMENT THE FOLLOWING BLOCK IF YOU NEEDED
####################################################################################
# Ground state Parameters
####################################################################################

# Ground state Morse Parameters
D0   =   0.17675
a0   =   1.0494
Req0 =   1.41691
# for Ground state only
InEmin     =    -0.18        # Starting E Minimum. it will automatically updated
InEmax     =    -0.16        # Starting E Maximum, it will automatically updated
# for Ground state :  gridshift =b when (InEmin to  a) in [[a,b],...]; a=gshift[0][0] b=gshift[0][1]     
gridShift = 8500         # default gridShift; need to cut end point's higher values; for lower vib.quantum numbers we need to increase it
gshift = [[InEmin,8500],[-0.1, 7000], [-0.05, 6000],[-0.01,50]]

##################################################################################
# Excited State Parameters
##################################################################################
'''
# Excited state Morse parameters, un comment if needed
D0   = 0.102928
a0   = 0.681859
Req0 = 2.00576

# for Excited state only
InEmin=-0.11                # Starting E Minimum. it will automatically updated
InEmax=-0.09                # Starting E Maximum, it will automatically updated

# for Excited state
gridShift = 8000
gshift = [[InEmin,8000],[-0.08,6000],[-0.02,4000],[-0.01,50]]
'''
##################################################################################

# define Morse function Here Zero-shifted asymptot is used
def Morse(x,D0,a0,Req0):
    return ( D0*(1-math.exp(-a0*(x - Req0)))**2 - D0)

