#=======================================================================
# IMPORTANT PARAMETERS  
#=======================================================================

grid      = 12000 #3200        # N division from -xp to xp ie. grid number
gridShift = 1000         # need to cut end point higher values
ynm1      = 1e-15       # far left point value
yn        = 2e-15       # ynm1+1th point value "representably small number"
xp        = 2.0         # range of oscillator as [-xp...0...xp] in atomic unit
Egap      = 0.01        # an empirical one to fit WFN at energy levels
node      = 0
EsearchMax =    0.1     # for plotting, y's limit

InEmin=0.0              # Starting E Minimum. it will automatically updated
InEmax=0.01             # Starting E Maximum, it will automatically updated
itern =10               # after this iteration program exit; if you want more levels increase it

# Important:
# PlotData- to save E,R,unNormalized Wfn upto grid shifted point, and then padded with zeros for reusing
# Storesdata as: PlotData[j][2] # (gridshifted,padded) WfnUnnormalized(Energy) ; PlotData[i][1] = x ;  Energy_i = PlotData[i][0][0]
global PlotData
PlotData = []

# define Potential function Here Zero-shifted asymptot is used

K=0.33
def Potential(K,x):
# here V=0 inside the box.
    return K*(x)*(x)/2
