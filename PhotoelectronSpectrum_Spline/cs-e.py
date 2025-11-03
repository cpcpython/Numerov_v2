# this interpolate (x y) dataset from the text file, e.txt (excited state data) to xfite.pkl and yfite.pkl pickle files
# Note: asymptotically zero adjusted data used (ie when r->large; y-> 0)
#

import pickle
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import json

# importing R and E(R) into these arrays
r = np.genfromtxt('e.txt', usecols=0 )
ey= np.genfromtxt('e.txt', usecols=1 )

x =  r
y = ey

grid=12000
# use bc_type = 'natural' adds the constraints as we described above
f = CubicSpline(x, y, bc_type='natural')
x_new = np.linspace(0, 12, grid+3) # 3 is needed for grid adjustments
y_new = f(x_new)

plt.figure(figsize = (10,8))
plt.plot(x_new, y_new, 'b')
plt.plot(x, y, 'ro')
plt.title('Cubic Spline Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# save fitted x,F(x)
# Open the file in binary read mode and load the array
try:
    with open('xfite.pkl', 'wb') as file:
        pickle.dump(x_new, file)
    with open('yfite.pkl','wb') as file:
        pickle.dump(y_new,file)
except Exception as e:
    print(f"An error occurred: {e}")
