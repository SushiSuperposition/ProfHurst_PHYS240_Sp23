# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:28:00 2023

@author: Matthew
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from lmfit import Model, Parameters


def poly(x, a, b):
    """
    Function returns a line with independent variable x
    :param x: independent variable
    :param a: slope
    :param b: intercept
    :return: a*x+b
    """
    return a*x+b*x**2

def misfit(speed, drag):
    # Euler-Cromer Method


    # Set initial position and velocity of the comet.
    r0 = 0.3
    v0 = 5.8
    h = r0 * v0
    r = np.array([r0, 0.0])
    v = np.array([0.0, v0])


    # Set physical parameters
    GM = 3.986  # Gravitational constant * Mass of Earth [*10^5 km^3/s^2]
    mass = 1.0  # Mass of missile (reference mass)
    time = 0.001
    
    if drag == 0:
        rho = 0

    # Loop over the desired number of steps using the specified numerical method.
    nStep = 560
    tau = .01

    xplot = np.empty(nStep)
    yplot = np.empty(nStep)
    laststep = nStep



    for iStep in tqdm(range(nStep)):
        
        # equations of motion, energy, and position
        E = (np.linalg.norm(v))**2 / 2 - GM / np.linalg.norm(r) 
        e = np.sqrt(1 + (2 * E * h**2) / GM**2)
        a = mass / 2*E
        x1 = time / a
        z1 = (x1)**2 / a
        C1 = (1 - np.cos(np.sqrt(z1))) / z1
        S1 = (np.sqrt(z1) - np.sin(np.sqrt(z1))) / np.sqrt(z1**3)
        t1 = ((np.linalg.norm(r) * np.linalg.norm(v)) / np.sqrt(GM) * x1**2 * C1 + (1 - np.linalg.norm(r) / a) * x1**3 * S1 + np.linalg.norm(r) * x1) / np.sqrt(GM)
        dt1 = ( x1**2 * C1 + (np.linalg.norm(r) * np.linalg.norm(v)) / np.sqrt(GM) * x1 * (1 - z1 * S1) ) / np.sqrt(GM)
        x2 = x1 + (time - t1) / dt1
        f = 1 - a / np.linalg.norm(r) * (1 - np.cos(x1 / np.sqrt(a)))
        g = t1 - x1**3 / (np.sqrt(GM) * S1)
        fdot = np.sqrt(GM) / (np.linalg.norm(r) * np.linalg.norm(r)) * x1 * (z1 * S1 - 1)
        gdot = 1 - a / np.linalg.norm(r) + a / np.linalg.norm(r) * np.cos(x1 / np.sqrt(a))
        
        # Record position for plotting
        xplot[iStep] = f * np.linalg.norm(r)
        yplot[iStep] = g * np.linalg.norm(v)

        # Calculate new position and velocity using the desired method
        accel = -GM * r / np.linalg.norm(r) ** 3
        r += tau * v  # Euler Step
        v += tau * accel
        
        time += tau
        
        
        
        
    xdata = xplot[:laststep]
    ydata = yplot[:laststep]

    
    model = Model(poly, independent_vars = ['x'])
    bestfit = model.fit(ydata, x=xdata, a = .01, b = .5)
    
    
    fig, ax = plt.subplots(figsize = (7, 5))
    ax.set_title('Trajectory of Missile')
    ax.plot(xplot[:laststep+1], yplot[:laststep+1], '.', label='Euler-Cromer method', markersize = 4)
    ax.plot(xplot[:len(bestfit.best_fit)], bestfit.best_fit, '-', label = 'Best fit', markersize = 10)
    ax.set_xlabel('Range (* 10^3 km)')
    ax.set_ylabel('Height (* 10^3 km)')
    plt.show()
    
    
    
traj = np.arange(5.8)
for t in traj:
    misfit(t, 0)

