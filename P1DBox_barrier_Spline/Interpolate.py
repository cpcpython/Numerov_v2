# this make a Potential file (V vs x) for the Spline Numerov method
# i/p: linear.dat which contains piecewise informatoin of the potential
# o/p: xfitg.pkl and yfitg.pkl contains interpolated grid points of x and V respectively
#      that is used by Numerov's code

import numpy as np
import re
import matplotlib.pyplot as plt
from numpy import loadtxt
import pickle

#from Parameters_Splinev1 import grid as grid
#from Parameters_Splinev1 import xp as xp
grid=12000
xp=8.0

lines = loadtxt("linear.dat", comments='#', delimiter="\t", unpack=False) # using whitespace delimiter


divX=len(lines)
print(divX)
subX = grid/divX

unitx=xp/grid

inn=0
AllX=[]
AllV=[]
for item in lines:
    xi=item[0];xf=item[1];  v=item[2]
    xp = [xi,xf]; fp = [v,v]
    x_new = np.arange(xi,xf,unitx ).tolist() #
    linearV=np.interp(x_new, xp, fp)
    AllX.extend(x_new)
    AllV.extend(linearV)

    inn=inn+1
# add 3 points extra
x_new = np.linspace(0., 8., grid+3) # 3 
# we need to add 3 extra points for Numerov method for the splines
AllX.append(AllX[len(AllX)-1]+(AllX[1]-AllX[0])*1); AllV.append(AllV[len(AllV)-1])
AllX.append(AllX[len(AllX)-1]+(AllX[1]-AllX[0])*2); AllV.append(AllV[len(AllV)-1])
AllX.append(AllX[len(AllX)-1]+(AllX[1]-AllX[0])*3); AllV.append(AllV[len(AllV)-1])

# saving to Pickles
f = open("test.csv", "w")

for i in range(len(AllX)):
    f.write("{:.8f} {:.18f}\n".format(AllX[i],AllV[i]))

f.close()

# Plotting the Potential
plt.figure(figsize = (10,8))
plt.plot(AllX, AllV, 'b')
plt.title('Spline Potential Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

rounded_AllX = [round(x, 20) for x in AllX]
rounded_AllV = [round(x, 20) for x in AllV]

# save fitted x,F(x)
# Open the file in binary read mode and load the array
try:
    with open('xfitg.pkl', 'wb') as file:
        pickle.dump(rounded_AllX, file)
    with open('yfitg.pkl','wb') as file:
        pickle.dump(rounded_AllV,file)
except Exception as e:
    print(f"An error occurred: {e}")

