# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 23:21:40 2017

@author: mexicore
"""

import numpy as np
import random as r
import math as ma
from pylab import * 
from operator import *
import random as r
from scipy.stats import norm



def brownian(x0, n, dt, delta, out=None):
    """\
    Generate an instance of Brownian motion (i.e. the Wiener process):

        X(t) = X(0) + N(0, delta**2 * t; 0, t)

    where N(a,b; t0, t1) is a normally distributed random variable with mean a and
    variance b.  The parameters t0 and t1 make explicit the statistical
    independence of N on different time intervals; that is, if [t0, t1) and
    [t2, t3) are disjoint intervals, then N(a, b; t0, t1) and N(a, b; t2, t3)
    are independent.
    
    Written as an iteration scheme,

        X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)


    If `x0` is an array (or array-like), each value in `x0` is treated as
    an initial condition, and the value returned is a numpy array with one
    more dimension than `x0`.

    Arguments
    ---------
    x0 : float or numpy array (or something that can be converted to a numpy array
         using numpy.asarray(x0)).
        The initial condition(s) (i.e. position(s)) of the Brownian motion.
    n : int
        The number of steps to take.
    dt : float
        The time step.
    delta : float
        delta determines the "speed" of the Brownian motion.  The random variable
        of the position at time t, X(t), has a normal distribution whose mean is
        the position at time t=0 and whose variance is delta**2*t.
    out : numpy array or None
        If `out` is not None, it specifies the array in which to put the
        result.  If `out` is None, a new numpy array is created and returned.

    Returns
    -------
    A numpy array of floats with shape `x0.shape + (n,)`.
    
    Note that the initial value `x0` is not included in the returned array.
    """

    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    
    r = norm.rvs(size=x0.shape + (n,), scale=delta*ma.sqrt(dt))


    # If `out` was not given, create an output array.
    if out is None:
        out = np.empty(r.shape)

    # This computes the Brownian motion by forming the cumulative sum of
    # the random samples. 
    np.cumsum(r, axis=-1, out=out)

    # Add the initial condition.
    out += np.expand_dims(x0, axis=-1)

    return out

def main(xmin,xmax,ymin,ymax):

    # The Wiener process parameter.
    delta = 0.25
    # Total time.
    T = 20.0
    # Number of steps.
    N = 1
    # Time step size
    dt = T/N
    TFrr = np.empty((200,200))
    #perimetro celular
    xcell = [-30,-30,30,30,-30]
    ycell = [-10,10,10,-10,-10]
    #area celular
    def insideCell(x,y):
        if x>-30 and x<30 and y>-10 and y<10: return True
        else: return False
    #elementos

    for i in range(200):

	# Initial values of x.

        x = np.empty((2,N+1))
        sol = np.empty((2,201))
        for j in range(len(x)):
            x[:,0] = 0.0

        sol[0][0] = 0.0
        sol[1][0] = 0.0
	TFr=[]
        for k in range(len(sol[0])-1):
            brownian(x[:,0], N, dt, delta, out=x[:,1:])
            while not(insideCell(x[0][1],x[1][1])):
                brownian(x[:,0], N, dt, delta, out=x[:,1:])
            sol[0][k+1] = x[0][1]
            sol[1][k+1] = x[1][1]
            x = np.empty((2,N+1))
            x[0][0] = sol[0][k+1]
            x[1][0] = sol[1][k+1]
            
    
            TF=0
	    # sumo uno si la trayectoria esta en el perimetro	    
	    if x[0][len(x)-1]>xmin and x[0][len(x)-1]<xmax and x[1][len(x)-1]>ymin and x[1][len(x)-1]<ymax:
       	    	TF += 1
            # junto cada trayectoria en un append   
            TFr.append (TF)	
        #creo una matriz con todas las trayectorias
        TFrr[i]=TFr
    TFr=[]
    #sumo cada TF para cada tiempo (cada columna de la matriz) 
    for i in xrange(200):
		    TF= 0
	            
		    for j in TFrr:
		       	TF+=j[i]
                    
                    TFr.append(TF)   
    	

    return (TFr)
    
    
def gill(K1, K_1,K2, K3, K4, K5, K6, TF, Pm, Pm_A, m, p, t, tmax, xmin,xmax,ymin,ymax):

    trange = []
    prange = []    
    TFr= main(xmin,xmax,ymin,ymax)

    tau=0         

    while t < tmax:
    	if int(t)>(t-tau) and int(t)%1==0:
        	TF = TFr[int(t)]     

        a =  [K1*Pm*TF, K_1*Pm_A, K2*Pm_A, K4*m, K3*m, K5*p, K6*Pm] 
        a0 = sum(a)	     
        r1 = r.random() 
        tau = -ma.log(r1)/a0    
        t = t+tau 
        r2 = r.random() 
        acumsum = np.cumsum(a)/a0   
        chosen_reaction = min([i for i in range(len(a)) if acumsum[i] >= r2]) 
 
        if chosen_reaction == 0:
               TF-=1; Pm-=1; Pm_A+=1;
        if chosen_reaction == 1:
               TF+=1; Pm+=1; Pm_A-=1;
        if chosen_reaction == 2:
               m+=1;
        if chosen_reaction == 3:
               m-=1;
        if chosen_reaction == 4:
               p+=1;
        if chosen_reaction == 5:
               p-=1;
        if chosen_reaction == 6:
               m+=1
        trange.append(t) 
        prange.append(p)
       
	

    return (trange,prange)
    