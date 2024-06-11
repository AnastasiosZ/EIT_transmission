# Importing all parameters
from parameters import *
from tqdm import tqdm
# Import the class from density_solver.py
from density_solver import *

def plot_elements():

    # For each of the nine plot elements
    for i in tqdm(range(9)):
        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ω0p,Ωc=Ω0c,Δ2=Δ2) # Instantiate a ρcl class
        ρ.sweep_Δ1(i,fid) # Sweep over single photon detuning Δ1
        ρ.plotρ(True,False) # and plot. In this case, the program saves the plots but doesn't show them

    # This below is optical and was uesd during the development of the tool to verify that an increased
    # control power lowered the EIT dip at two-photon resonance 
    print("Imaginary part of ρ31 at two-photon resonance for Pp=%.2fμW, Pc=%.2fmW is %.2E" % (Pp*1e6,Pc*1e3,Decimal(np.imag(ρ.ρsol(0)[6]))))
    print("Real part of ρ31 at two-photon resonance for Pp=%.2fμW, Pc=%.2fmW is %.2E" % (Pp*1e6,Pc*1e3,Decimal(np.real(ρ.ρsol(0,)[6]))))

# Calling the function plots the nine density matrix elements
plot_elements()