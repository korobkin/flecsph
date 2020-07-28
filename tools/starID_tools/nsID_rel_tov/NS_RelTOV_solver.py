#!/usr/bin/env python3

import numpy as np
import h5py
import sys
import multiprocessing as multi
import scipy.optimize as opt

from matplotlib import pylab
from scipy.constants import pi, G, c, hbar, m_n, m_e     # Physical constants
from scipy.interpolate import interp1d
from random import uniform

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

# User input parameters in CGS
rho_center     = float(sys.argv[1]) #1e18                          # Central density [g/cm^3],
dr             = float(sys.argv[2]) #10                           # Radial step [cm]

# Convert physical constants from MKS into CGS
Gr   = G*1.e3                                 # Grav. constant from m^3/(kg*s^2) to cm^3/(g*s^2)
hb   = hbar*1.e7                              # from kg*m^2/s to g*cm^2/s
me   = m_e*1.e3                               # from kg to g
mn   = m_n*1.e3                               # from kg to g            
cs   = c*1.e2                                 # from m/s to cm/s
Msun = 1.98892e33                             # Solar mass [g]

# Proton fraction
Ye = 0.5                                      # Proton fraction

# Piecewise polytrope equations of state
Gamma0 = 5.0/3.0
Gamma1 = 2.5
rho1 = 1e17

K0 = (3.0*pi**2)**(2.0/3.0)*hb**2/(5.0*mn**(8.0/3.0))
KR = (3.0*pi**2)**(1.0/3.0)*hb*cs/(4.0*mn**(4.0/3.0))
P1 = K0*rho1**Gamma0
K1 = P1/rho1**Gamma1

# Interpolation of the mass density and the pressure for a 
# SLy4 EOS in taubulated form. 
TabEos = True # Read eos from table 
if TabEos == True:
  p_conv_fact = 1.6022e33 #[Mev/fm^3] to [dyn/cm^2] set as one if it is already give in [dyn/cm^2]
  eos = "eos_data/sly4d_eos.dat"
  with open(eos, "rb") as eos_table:
      density_col = eos_table.read().split()[1::3]
  with open(eos, "rb") as eos_table:
      pressure_col = eos_table.read().split()[2::3]
  for i in range(len(density_col)):
      density_col[i] = float(density_col[i])
      pressure_col[i] = p_conv_fact*float(pressure_col[i])
  pressure_interp = interp1d(density_col, pressure_col, fill_value="extrapolate")
  density_interp = interp1d(pressure_col, density_col, fill_value="extrapolate")
  
  # Getting density 
  def get_density_tab(rho):
    return density_interp(rho)

  # Getting pressure
  def get_pressure_tab(P):
    return pressure_interp(P)

############ START ##############
print("---- TOV Equation for Neutron Star ID ----")
print("")
print("Using:")
print("Gravitational constant: {} cm^3/(g*s^2)".format(Gr))
print("Neutron mass:           {} g".format(mn))
print("Electron mass:          {} g".format(me))
print("Speed of light:         {} cm/s".format(cs))
print("Planck constant:        {} g*cm^2/s".format(hb))
print("Polytorpe values in cgs:{}, {}, {}".format(K0,KR,K1))
print("")

# pressure (difference) for given mass density 
def eos(rho):
     if rho < rho1:
         return K0*rho**Gamma0
     else:
         return K1*rho**Gamma1

def inv_eos(P):
    if P < P1:
        return (P/K0)**(1.0/Gamma0)
    else:
        return (P/K1)**(1.0/Gamma1)


# 4th order Runge-Kutta routines for integration
def rk4(f,y,x,h,rho):
	k1 = f(y,x,rho)*h
	k2 = f(y + 0.5*k1, x + 0.5*h, rho)*h
	k3 = f(y + 0.5*k2, x + 0.5*h, rho)*h
	k4 = f(y + k3, x + h, rho)*h
	return y + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0

# Define TOV equation
def tov(y,r,rho):
    if TabEos == True:
        P = get_pressure_tab(y[0]) # Pressure from table
        rho = get_pressure_tab(rho) # Density from table
        m = y[1] 
    else:
	    P, m = y[0], y[1]                             # pressure and mass
	    rho  = inv_eos(P)                             # Density from inverse of EOS
    dPdr = -Gr*rho*m/(cs**2*r**2)*(1.0+P/(rho*cs**2))*(1.+4.*pi*r**3*P/(m*cs*cs))*1./(1.-2.*Gr*m/(cs*cs*r))  # P eqn
    dmdr = 4.0*pi*rho*r**2                               # M eqn
    return pylab.array([dPdr, dmdr])

# For calculation purpose, define derivative of quantities with respect to r
def get_dfdr(f,dr):
	dfdr       = np.empty_like(f)
	dfdr[1:-1] = (f[2:] - f[:-2])/(2*dr)
	dfdr[0]    = (f[2] - f[0])/dr
	dfdr[-1]   = (f[-1] - f[-2])/dr
	return dfdr

# Solve TOV equation within a given central density
def sol_tov(rho_center):
    rmin = dr
    rmax = 20000.0

    r = pylab.arange(rmin,rmax+dr,dr)
    m = pylab.zeros_like(r)
    P = pylab.zeros_like(r)
    rho = pylab.zeros_like(r)
    i = 0
    rho[i] = rho_center
    m[i]   = (4.0/3.0)*pi*rho_center*r[0]**3
    P[i]   = eos(rho_center)
    y = pylab.array([P[i],m[i]])    
    while P[i]>0.0 and i<len(r)-1:
        y = rk4(tov,y,r[i],dr,rho[i])
        m[i+1] = y[1]
        P[i+1] = y[0]
        i = i+1
    rho = pylab.array(list(map(lambda p: inv_eos(p),P)))
    if P[i]<0.0:
        P[i] = 0.0
        rho[i] = 0.0
    m,r,rho,P = m[:i],r[:i],rho[:i],P[:i]

    return m,m[-1]/Msun,r,rho,P

# Solving the equations
mball, M, r, rho, P = sol_tov(rho_center)

j=0
drhodr = get_dfdr(rho,dr)
i = len(r)-1
print_every = int(round(i/40))
Rtot = r[i]
Mtot = mball[i]

print("Resulting Star:")
print("Mtot: {} g; Rtot {} cm".format(Mtot, Rtot))
print("rho_central: {} g/cm^3".format(rho[0]))
print("P_central: {} g/cm^3".format(P[0]))
print(" ")
print("c1: r/Rtot, c2: rho*Rtot^3/Mtot, c3: m/Mtot, c4: (drho/dr)*Rtot^4/Mtot")
print('{:<18.11E}{:<18.11E}{:<18.11E}{:<18.11E}'.format(0.0, rho[0]*Rtot**3/Mtot, 0.0, 0.0))
for j in range(0, i, print_every):
 print('{:<18.11E}{:<18.11E}{:<18.11E}{:<18.11E}'.format(r[j]/Rtot, rho[j]*Rtot**3/Mtot, mball[j]/Mtot, drhodr[j]*Rtot**4/Mtot))
print('{:<18.11E}{:<18.11E}{:<18.11E}{:<18.11E}'.format(r[i]/Rtot, rho[i]*Rtot**3/Mtot, mball[i]/Mtot, drhodr[i]*Rtot**4/Mtot))
print("Done")

