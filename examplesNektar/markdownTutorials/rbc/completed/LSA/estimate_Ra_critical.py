import numpy as np
from scipy.interpolate import splrep, splev

data = np.loadtxt('data.d')
Ra = data[:,0]
eigen_value = data[:,1]





############ Spline ###########################

A = splrep(eigen_value,Ra,k=1)
Ra_c = splev(0,A)
print("Ra critical = ",Ra_c)
print(np.interp(0, eigen_value, Ra))
