from scipy import constants
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

pi = constants.pi # Pi
ε0 = constants.epsilon_0 # Permittivity of free space (F/m)
hbar = constants.hbar # Reduced Planck constant (Js)
c = 299792458 # Speed of light (m/s)


fid = int(1e4) # Program fidelity

dc= 2e-3 # Control laser beam diameter (m)
dp= 2e-3 # Control laser beam diameter (m)
Ac = pi*(dc/2)**2 # Laser spot size of control (m^2)
Ap = pi*(dp/2)**2 # Laser spot size of control (m^2)
n = 1 # Average index of refraction
Pc = 1.0e-3 # Control laser power (W)
Pp = 20.0e-6 # Probe laser power (W)
E0c = (2*Pc/(c*ε0*Ac*n))**0.5 # Amplitude of control field (V/m)
E0p = (2*Pp/(c*ε0*Ap*n))**0.5 # Amplitude of probe field (V/m)

dip = -2.537e-29 # D1 transition dipole moment (Cm)
dip1 = (1/4)**0.5*dip # <1|μ|3> Dipole moment (Cm)
dip2 = (1/12)**0.5*dip # <2|μ|3> Dipole moment (Cm)



λ0p =  795.0e-9 # Vacuum wavelength of the probe field (m)
λ0c =  795.0e-9 # Vacuum wavelength of the control field (m)
k0p = 2*pi/λ0p # Vacuum wave number of tuned probe field (1/m)
k0c = 2*pi/λ0c # Vacuum wave number of tuned probe field (1/m)


Ω0c = -dip2*E0c/hbar # Control Rabi frequency (rad Hz)
Ω0p = -dip1*E0p/hbar # Probe Rabi frequency (rad Hz)


Γ3 = 36.1e6 # Spontaneous emission rate out of state |3> (Hz) 
Γ32 = Γ3/2 # Spontaneous emission rate |3>->|2> (Hz)
Γ31 = Γ3/2 # Spontaneous emission rate |3>->|1> (Hz)
γ3d = 0.1*Γ3 # Dephasing rate of state |3> (Hz)
# γ2d = 1e3 # Dephasing rate of state |2> (Hz)
γ2d = 1e6 # Dephasing rate of state |2> (Hz)


Δ1 = 0 # Detuning of probe field (Hz)
Δ2 = 0 # Detuning of control field (Hz)


T=298.15 # Temperature of the Rb-87 vapour sample (K)
N = 9.715e15 # Atomic number density of the Rb-87 sample (1/m^3)


L = 5e-2 # Length cell (m)


σ = 3*λ0p**2/(2*pi) # Resonant cross section (m^2)
optical_depth_per_L = N*σ # Optical depth/density per length (1/m)

# print("Ω0p/Ω0c = %.3f" % (Ω0p/Ω0c))
# print("E0p/E0c = %.3f" % (E0p/E0c))

# print(Pp/Pc,Pp,Pc,E0p,E0c,Ω0p,Ω0c)
# print(Pp/Pc)

print('%.3E' % Decimal(Ω0c**2))
print('%.3E' % Decimal((Γ3+γ3d)*γ2d))