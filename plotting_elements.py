from parameters import *
from tqdm import tqdm
from density_solver import *

def plot_elements():

    for i in tqdm(range(9)):
        ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ω0p,Ωc=Ω0c,Δ2=Δ2)
        ρ.sweep_Δ1(i,fid)
        ρ.plotρ(True,False)

    print("Imaginary part of ρ31 at two-photon resonance for Pp=%.2fμW, Pc=%.2fmW is %.2E" % (Pp*1e6,Pc*1e3,Decimal(np.imag(ρ.ρsol(0)[6]))))
    print("Real part of ρ31 at two-photon resonance for Pp=%.2fμW, Pc=%.2fmW is %.2E" % (Pp*1e6,Pc*1e3,Decimal(np.real(ρ.ρsol(0,)[6]))))


plot_elements()