# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:09:33 2023

@author: Matthew
"""

# orbit - Program to compute the orbit of a comet.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm



# Runge-Kutta method


def rk4(x, t, tau, derivsRK, param):
    """
    Runge-Kutta integrator (4th order)
    Input arguments
    :param x: current value of dependent variable
    :param t: independent variable (usually time)
    :param tau: step size (usually time step)
    :param derivsRK: right hand side of the ODE; derivsRK is the name of the function which returns dx/dt
    Calling format derivsRK (x, t, param).
    :param param: estra parameters passed to derivsRK
    :return:
    xout: new value of x after a step of size tau
    """

    half_tau = 0.5*tau
    F1 = derivsRK(x, t, param)
    t_half = t + half_tau
    xtemp = x + half_tau*F1
    F2 = derivsRK(xtemp, t_half, param)
    xtemp = x + half_tau*F2
    F3 = derivsRK(xtemp, t_half, param)
    t_full = t + tau
    xtemp = x + tau*F3
    F4 = derivsRK(xtemp, t_full, param)
    xout = x + tau/6.0 * (F1 + F4 + 2.0*(F2+F3))

    return xout

# Define gravrk function used by the Runge-Kutta routines
def gravrk(s, t, GM):
    """
    Returns the right-hand side of the Kepler ODE; used by Runge-Kutta routines
    :param s: State vector [r(0), r(1), v(0), v(1)]
    :param t: Time (not used here, included to match derivsRK input)
    :param GM: Parameter G*M - gravitational constant * solar mass Units [AU^3/yr^2]
    :return: deriv: Derivatives [dr(0/dt), dr(1)/dt, dv(0)/dt, dv(1)/dt]
    """

    # Compute acceleration
    r = s[:2]  # Unravel the vector s into position and velocity
    v = s[2:]
    accel = -GM * r / 0.5**3  # Gravitational acceleration

    # Return derivatives
    deriv = np.array([v[0], v[1], accel[0], accel[1]])

    return deriv


# Set initial position and velocity of the comet.
r0 = .1
v0 = 5.8
r = np.array([r0, 0])
v = np.array([0.0, v0])



state = np.array([r[0], r[1], v[0], v[1]])  # State used by R-K routines

# Set physical parameters
GM = 3.98600  # Gravitational constant * Mass of Earth [*10^5 km^3/s^2]
time = 0.0

# Loop over the desired number of steps using the specified numerical method.
nStep = 555
tau = .001

rplot = np.empty(nStep)
thplot = np.empty(nStep)
tplot = np.empty(nStep)
stoplist = []


for iStep in tqdm(range(nStep)):

    # Record position and energy for plotting
    
    

    rplot[iStep] = np.linalg.norm(r)  # Record radial position and angle for polar plot
    
    
    # print(rplot[iStep])
    
    tplot[iStep] = time
    
    state = rk4(state, time, tau, gravrk, GM)
    r = state[:2]  # 4th Order Runge-Kutta
    v = state[2:]
    time += tau
    


fig = plt.figure(figsize=(10.0, 5.25))
ax = fig.add_subplot(121)
ax.plot(tplot, rplot, '+',)
ax.set_title('Trajectory of Missile')
ax.set_xlabel('Range (* 10^3 km)')
ax.set_ylabel('Height (* 10^3 km)')
ax.grid(True)
fig.tight_layout(pad=5.0)


plt.show()
