import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from density_solver import *
from tqdm import tqdm


def co_prop_save(Ω0p,Ω0c,E0p,E0c):

    Es = []

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ω0c; Ec = E0c # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0p*c+Δ2) # of probe and control field
    d = L/fid # Distance between consecutive propagation points (m)

    for xi in tqdm(np.linspace(0,L,fid)):

        Es.append([Ep,Ec])

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


    return np.array(Es)

def counter_prop(Ω0p,Ω0c_guess,E0p,E0c_guess):

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ω0c_guess; Ec = E0c_guess # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0p*c+Δ2) # of probe and control field
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
        rc = np.exp(αc*d/2) # Field dissipation ratio
        Ωp*=rp; Ep*=rp # Dissipation of field amplitude and Rabi frequency
        Ωc*=rc; Ec*=rc # Dissipation of field amplitude and Rabi frequency

    return Ωc,Ec
    

def counter_prop_save(Ω0p,Ω0c,E0p,E0c):

    Ec_guess = E0c
    Ωc_guess = Ω0c

    Ωc,Ec = counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)

    while abs(Ec-E0c)/E0c>1e-7:

        # print(abs(Ec-E0c)/E0c)

        Ec_guess += (E0c-Ec)
        Ωc_guess = -dip2*Ec_guess/hbar

        Ωc,Ec = counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)
    

    Es = []

    Ωp = Ω0p; Ep = E0p # Setting initial Rabi frequencies 
    Ωc = Ωc_guess; Ec = Ec_guess # and field amplitudes
    λp = 2*pi*c/(k0p*c+Δ1) # Accounting for detuning
    λc = 2*pi*c/(k0c*c+Δ2) # of probe and control field    
    d = L/fid # Distance between consecutive propagation points (m)

    for xi in tqdm(np.linspace(0,L,fid)):

        Es.append([Ep,Ec])

        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2)
        ρsol = ρ.ρsol(Δ=Δ1) # Get solution for particular probe detuning Δ1
        ρ31=ρsol[6]; ρ32=ρsol[7]

        χp = N*dip1*ρ31/(ε0*Ep) # Susceptibility at probe frequency
        χc = N*dip2*ρ32/(ε0*Ec) # Susceptibility at probe frequency 

        αp = 4*pi*np.imag(np.sqrt(1+χp))/λp # Absorption coefficient for probe field
        αc = 4*pi*np.imag(np.sqrt(1+χc))/λc # Absorption coefficient for probe field

        rp = np.exp(-αp*d/2) # Field dissipation ratio
        rc = np.exp(αc*d/2) # Field dissipation ratio
        Ωp*=rp; Ep*=rp # Dissipation of field amplitude and Rabi frequency
        Ωc*=rc; Ec*=rc # Dissipation of field amplitude and Rabi frequency


    return np.array(Es)



from parameters import *

Esco = co_prop_save(Ω0p,Ω0c,E0p,E0c)
Escn = counter_prop_save(Ω0p,Ω0c,E0p,E0c)

x = np.linspace(0,L,fid)
fig, (ax1, ax2) = plt.subplots(1, 2)

fig.set_figheight(5)
fig.set_figwidth(13)
fig.set_dpi(30)

ax1.ticklabel_format(scilimits=(0,0))
ax2.ticklabel_format(scilimits=(0,0))

ax1.plot(x,Esco[:,0],c='orangered',linewidth=7,label=r'Co')
ax1.plot(x,Escn[:,0],c='steelblue',linewidth=3,label=r'Counter')
ax1.set_xlabel(r'Distance along the cell (m)')
ax1.set_ylabel(r'$E_p$ (V/m)')
ax1.set_title(r'Probe')
ax1.legend()

ax2.plot(x,Esco[:,1],c='orangered',linewidth=3,label=r'Co')
ax2.plot(x,Escn[:,1],c='steelblue',linewidth=3,label=r'Counter')
ax2.set_xlabel(r'Distance along the cell (m)')
ax2.set_ylabel(r'$E_c$ (V/m)')
ax2.set_title(r'Control')
ax2.legend()

ax1.grid()
ax2.grid()


fig.suptitle(r'Propagation: $L=%.1f$ cm, $P_p=%.1f$ \textmu W, $P_c=%.1f$ mW' % (L*100,Pp*1e6,Pc*1e3))
fig.savefig("propagation/L=%.1fcm_Pp=%.1fμW_Pc=%.1fmW.png" % (L*100,Pp*1e6,Pc*1e3))


plt.show()

