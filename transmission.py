import numpy as np
# Import the class from density_solver.py
from density_solver import *
from tqdm import tqdm

# Importing all parameters
from parameters import *


# Co-propagation transmission

def co_T(Ω0p,Ω0c,E0p,E0c):

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ω0c; Ec = E0c # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0c*c+Δ2) # of probe and control field
    d = L/fid # Distance between consecutive finite step points (m)

    for xi in np.linspace(0,L,fid):

        # Instantiating class at every step is necessary because Rabi frequencies evolve
        # which changes the solution of the density matrix elements
        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2)
        ρsol = ρ.ρsol(Δ=Δ1) # Get solution for particular probe detuning Δ1
        ρ31=ρsol[6]; ρ32=ρsol[7] # Choosing solutions of elements ρ31 and ρ32

        χp = N*dip1*ρ31/(ε0*Ep) # Susceptibility at probe frequency
        χc = N*dip2*ρ32/(ε0*Ec) # Susceptibility at control frequency 

        αp = 4*pi*np.imag(np.sqrt(1+χp))/λp # Absorption coefficient for probe field
        αc = 4*pi*np.imag(np.sqrt(1+χc))/λc # Absorption coefficient for control field

        rp = np.exp(-αp*d/2) # Probe field dissipation ratio
        rc = np.exp(-αc*d/2) # Control field dissipation ratio
        Ωp*=rp; Ep*=rp # Dissipation of probe field amplitude and Rabi frequency
        Ωc*=rc; Ec*=rc # Dissipation of control field amplitude and Rabi frequency

    # Transmission is the ratio of final and initial intensities
    # or equivalently, the square of the ratio of final and initial field amplitudes

    Tp = (Ep/E0p)**2
    Tc = (Ec/E0c)**2

    return Tp,Tc # Return probe and control transmissions

# Counter-propagation

def counter_prop(Ω0p,Ω0c,E0p,E0c):

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ω0c; Ec = E0c # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0c*c+Δ2) # of probe and control field
    d = L/fid # Distance between consecutive finite step points (m)

    for xi in np.linspace(0,L,fid):

        # Instantiating class at every step is necessary because Rabi frequencies evolve
        # which changes the solution of the density matrix elements
        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2)
        ρsol = ρ.ρsol(Δ=Δ1) # Get solution for particular probe detuning Δ1
        ρ31=ρsol[6]; ρ32=ρsol[7] # Choosing solutions of elements ρ31 and ρ32


        χp = N*dip1*ρ31/(ε0*Ep) # Susceptibility at probe frequency 
        χc = N*dip2*ρ32/(ε0*Ec) # Susceptibility at control frequency 

        αp = 4*pi*np.imag(np.sqrt(1+χp))/λp # Absorption coefficient for probe field
        αc = 4*pi*np.imag(np.sqrt(1+χc))/λc # Absorption coefficient for control field

        rp = np.exp(-αp*d/2) # Probe field dissipation ratio
        rc = np.exp(αc*d/2) # Control field amplification ratio
        Ωp*=rp; Ep*=rp # Dissipation of probe field amplitude and Rabi frequency 
        Ωc*=rc; Ec*=rc # Amplification of control field amplitude and Rabi frequency

    return Ωp,Ep,Ωc,Ec


def counter_T(Ω0p,Ω0c,E0p,E0c):

    # The first guess for the final field amplitude and Rabi frequency
    # is the same as the initial values, which is obviously wrong but
    # in the EIT regime (Ωc>>Ωp) the control field experiences small absorption
    # This first guess lowers convergence time 
    Ec_guess = E0c
    Ωc_guess = Ω0c

    # Control field amplitude and Rabi frequency are extrapolated from the guess
    # by amplifying the control field instead of dissipating it and comparing the
    # predicted value with the known
    Ωp,Ep,Ωc,Ec = counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)

    # While the predicted and knwon control field amplitude have a relative difference of
    # 10^(-7), adjust the guess
    while abs(Ec-E0c)/E0c>1e-7:

        # If the program takes too long to execute, the first place the user should
        # look for debugging is here, because sometimes instead of converging,
        # the guess diverges or oscillates from the known

        # print(abs(Ec-E0c)/E0c)
        
        Ec_guess += (E0c-Ec) # Aggressive guess adjustment 
        Ωc_guess = -dip2*Ec_guess/hbar # Getting control Rabi frequency from field amplitude

        Ωp,Ep,Ωc,Ec = counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)


    # Transmission is the ratio of final and initial intensities
    # or equivalently, the square of the ratio of final and initial field amplitudes
    Tp = (Ep/E0p)**2 
    Tc = (Ec_guess/E0c)**2

    return Tp,Tc # Return probe and control transmissions


# Printing co- and counter-propagation probe and control transmissions

print("Co-propagation : (probe,control) transmission : (%f,%f)" % (co_T(Ω0p,Ω0c,E0p,E0c)))
print("Counter-propagation : (probe,control) transmission : (%f,%f)" % (counter_T(Ω0p,Ω0c,E0p,E0c)))
