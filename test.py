import numpy as np
import mathplotlib.pyplot as plt


# Constants
gravity = 9.8
drag_coeff = 0.53
rho_air = 1.2


# Parameters
dt = 1e-5

mass = 58e-3
radius = 3.4e-2

t0 = 0.0
rx0 = 0.0
ry0 = 50

v0 = 20.0
theta = 0.0

vx0 = v0 * np.cos(theta)
vy0 = v0 * np.sin(theta)

ax0 = 0.0
ay0 = -g


# Pre-compute
cross_section = np.pi * radius**2


# Compute loop
tt = [t0]

rxrx = [rx0]
ryry = [ry0]

vxvx = [vx0]
vyvy = [vy0]

axax = [ax0]
ayay = [ay0]


while ryry[-1] > 0.0:

    # Update time
    tt += [tt[-1] + dt]

    # Drag force
    vsquared = vxvx[-1]**2 + vyvy[-1]**2
    vtheta = np.abs(np.atan(vyvy[-1] / vxvx[-1]))

    fdrag = 0.5 * drag_coeff * rho_air * cross_section * vsquared
    fdragx = -fdrag * np.cos(vtheta)
    fdragy = +fdrag * np.sin(vtheta)

    # Update acceleration
    axax += [fdragx / mass]
    ayay += [fdragy / mass - g]

    # Update velocity
    vxvx += [vxvx[-1] + dt * axax[-1]]
    vyvy += [vyvy[-1] + dt * ayay[-1]]

    # Update position
    rxrx += [rxrx[-1] + dt * vxvx[-1]]
    ryry += [ryry[-1] + dt * vyvy[-1]]



print(f"\nFinal values: xf = {rxrx[-1]}, yf = {ryry[-1]}, tf = {tt[-1]}")

