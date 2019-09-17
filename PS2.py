import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
import pandas as pd



def modeldefs(Xm, Xn, A, *mparams):
    # unpack sets
    k = Xm                            # k(t)
    kp = Xn                           # k(t+1)
    a = A                             # z(t)
    # find variables
    y = a*k**alpha                   # y(t)  
    r = alpha*a*k**(alpha-1)                
    c = y + (1-delta)*k - kp         # c(t)
    i = y - c                        # i(t)
    
    if c<=0:
        c = max(c,0.00000001)
    if sigma == 1.0:                 # u(t)
        u = np.log(c)
    else:
        u = c**(1-sigma)/(1-sigma)
    if np.isnan(u) or np.isinf(u):
        u = -1.0E+99
    return y,r,c,i,u

alpha = .35
beta = .99
sigma = 2
delta = .025

mparams = (alpha, beta, delta, sigma)

# find the steady state
abar = 1
rbar = 1/beta + delta - 1
kbar = (alpha*abar/(rbar))**(1/(1-alpha))
#ybar = abar*kbar**alpha
#cbar = ybar - delta*kbar
#ubar = (cbar**(1-sigma)-1)/(1-sigma)

# set up grid for k

klow = 0.01*kbar
khigh = 1.99*kbar
knpts = 21 
kgrid = np.linspace(klow, khigh,knpts)

# set up grid for A
anpts = 2 
agrid = [1.1,0.678] # state space of A
Pimat = np.array([[0.977,1-0.977],[1-0.926,0.926]]) # hh = 0.977


# initialize arrays
VF = np.zeros((knpts, anpts))
VFnew = np.zeros((knpts, anpts))
PF = np.zeros((knpts, anpts))

# set parameter values for the iterative prodcedure.
ccrit = 1.0E-6 #  definition of "close" of VF and VFnew 
maxit = 2000  # maximum number of iterations allowed
damp = 1.      # weight we put on the new value function relative to the old one
dist = 1.0E+99 # initial distance measure and needs to be larger than ccrit.
iters = 0      # nitial count for the number of iterations peformed

# starting A = high

while (dist > ccrit) and (iters < maxit):
    VFnew = np.zeros((knpts, anpts))
    iters = iters + 1
    for i in range (0, knpts):
        for j in range(0, anpts):
            maxval = -1.0E+98
            
            for m in range(0, knpts):
                # get current period utility
                yout, rat, con, inv, u =  \
                    modeldefs(kgrid[i], kgrid[m], agrid[j], *mparams)
                # get expected value
                val = 0.
              
                val = val + u + beta*(np.dot(Pimat[j],VF[m]))
                
                    # if this exceeds previous maximum do replacements
                if val > maxval:
                    maxval = val
                    VFnew[i, j] = val
                    PF[i, j] = kgrid[m]
                    
    # dist = np.mean(np.abs(VF - VFnew))
    dist = np.amax(np.abs(VF - VFnew))
    print('iteration: ', iters, 'distance: ', dist)
    VF = VFnew 
    
PFdata = pd.DataFrame(PF)
plt.plot(kgrid,PFdata[1])
plt.plot(kgrid,PFdata[0])

VFdata = pd.DataFrame(VF)    
a = VFdata[1]
b = VFdata[0]
plt.plot(kgrid, a, 'r') # plotting t, a separately 
plt.plot(kgrid, b, 'b') # plotting t, b separately 
plt.show()
