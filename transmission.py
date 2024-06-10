import numpy as np
from density_solver import *
from tqdm import tqdm

from parameters import *

def co_T(Ω0p,Ω0c,E0p,E0c):

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ω0c; Ec = E0c # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0c*c+Δ2) # of probe and control field
    d = L/fid # Distance between consecutive propagation points (m)

    for xi in np.linspace(0,L,fid):
        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2)
        ρsol = ρ.ρsol(Δ=Δ1) # Get solution for particular probe detuning Δ1
        ρ31=ρsol[6]; ρ32=ρsol[7]

        χp = N*dip1*ρ31/(ε0*Ep) # Susceptibility at probe frequency
        χc = N*dip2*ρ32/(ε0*Ec) # Susceptibility at probe frequency 

        αp = 4*pi*np.imag(np.sqrt(1+χp))/λp # Absorption coefficient for probe field
        αc = 4*pi*np.imag(np.sqrt(1+χc))/λc # Absorption coefficient for probe field

        rp = np.exp(-αp*d/2) # Field dissipation ratio
        rc = np.exp(-αc*d/2) # Field dissipation ratio
        Ωp*=rp; Ep*=rp # Dissipation of field amplitude and Rabi frequency
        Ωc*=rc; Ec*=rc # Dissipation of field amplitude and Rabi frequency

    Tp = (Ep/E0p)**2
    Tc = (Ec/E0c)**2

    return Tp,Tc


def counter_prop(Ω0p,Ω0c,E0p,E0c):

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ω0c; Ec = E0c # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0c*c+Δ2) # of probe and control field
    d = L/fid # Distance between consecutive propagation points (m)

    for xi in np.linspace(0,L,fid):

        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2)
        ρsol = ρ.ρsol(Δ=Δ1) # Get solution for particular probe detuning Δ1
        ρ31=ρsol[6]; ρ32=ρsol[7]


        χp = N*dip1*ρ31/(ε0*Ep) # Susceptibility at probe frequency 
        χc = N*dip2*ρ32/(ε0*Ec) # Susceptibility at probe frequency 

        αp = 4*pi*np.imag(np.sqrt(1+χp))/λp # Absorption coefficient for probe field
        αc = 4*pi*np.imag(np.sqrt(1+χc))/λc # Absorption coefficient for probe field

        rp = np.exp(-αp*d/2) # Field dissipation ratio
        rc = np.exp(αc*d/2) # Field dissipation ratio, counter-propagation
        Ωp*=rp; Ep*=rp # Dissipation of field amplitude and Rabi frequency
        Ωc*=rc; Ec*=rc # Dissipation of field amplitude and Rabi frequency

    return Ωp,Ep,Ωc,Ec # If lower limit of Ep is not reached, return it


def counter_T(Ω0p,Ω0c,E0p,E0c):

    Ec_guess = E0c
    Ωc_guess = Ω0c

    Ωp,Ep,Ωc,Ec = counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)
    # Eclast = 1e10

    while abs(Ec-E0c)/E0c>1e-7:

        # print(abs(Ec-E0c)/E0c)

        # if Ec-E0c>Eclast: print(Ec-E0c)
        
        Ec_guess += (E0c-Ec)
        Ωc_guess = -dip2*Ec_guess/hbar

        # Eclast=Ec-E0c

        Ωp,Ep,Ωc,Ec = counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)

    # print(abs(Ec-E0c)/E0c)


    Tp = (Ep/E0p)**2 # Transmission is the square of the ratio of output and input field amplitudes
    Tc = (Ec_guess/E0c)**2 # Transmission is the square of the ratio of output and input field amplitudes

    return Tp,Tc

# print("Co-propagation : (probe,control) transmission : (%f,%f)" % (co_T(Ω0p,Ω0c,E0p,E0c)))
# print("Counter-propagation : (probe,control) transmission : (%f,%f)" % (counter_T(Ω0p,Ω0c,E0p,E0c)))
