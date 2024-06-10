from parameters import *
from transmission import co_T,counter_T
from tqdm import tqdm

ns=200

def plot_T_varPc():


    filename = "var_Pc"


    Ppratios = np.logspace(-5,0.1,ns)


    Tpco = []; Tpcn = []
    Tcco = []; Tccn = []

    try:
        with open('npy/Tpco%s.npy' % filename, 'rb') as f: Tpco=np.load(f)
        with open('npy/Tpcn%s.npy' % filename, 'rb') as f: Tpcn=np.load(f)
        with open('npy/Tcco%s.npy' % filename, 'rb') as f: Tcco=np.load(f)
        with open('npy/Tccn%s.npy' % filename, 'rb') as f: Tccn=np.load(f)
    except (FileNotFoundError,NameError) as e:
        for rat in tqdm(Ppratios):
            
            Pci = Pp/rat
            E0ci = (2*Pci/(c*ε0*Ac*n))**0.5 # Amplitude of probe field (V/m)
            Ω0ci = -dip2*E0ci/hbar # Probe Rabi frequency (rad Hz)

            temp = co_T(Ω0p,Ω0ci,E0p,E0ci)
            Tpco.append(temp[0])
            Tcco.append(temp[1])
            temp = counter_T(Ω0p,Ω0ci,E0p,E0ci)
            Tpcn.append(temp[0])
            Tccn.append(temp[1]) 


        with open('npy/Tpco%s.npy' % filename, 'wb') as f: np.save(f,Tpco)
        with open('npy/Tpcn%s.npy' % filename, 'wb') as f: np.save(f,Tpcn)
        with open('npy/Tcco%s.npy' % filename, 'wb') as f: np.save(f,Tcco)
        with open('npy/Tccn%s.npy' % filename, 'wb') as f: np.save(f,Tccn)


    plt.ticklabel_format(scilimits=(0,0))
    plt.plot(Ppratios,Tpco,c='orangered', linewidth=7,label=r'Probe: Co')
    plt.plot(Ppratios,Tpcn,c='steelblue', linewidth=3,label=r'Probe: Counter')
    plt.plot(Ppratios,Tcco,c='darkorchid',linewidth=7,label=r'Control: Co')
    plt.plot(Ppratios,Tcco,c='forestgreen',linewidth=3,label=r'Control: Counter')
    plt.grid()
    plt.legend(fontsize=17)
    plt.title(r'Transmission, varying $P_c$: $L=%.1f$ cm, $P_p=%.1f$ \textmu W' % (L*100,Pp*1e6),fontsize=20)
    plt.xlabel(r'$P_p/P_c$',fontsize=20)
    plt.xscale('log')
    plt.ylabel(r'Transmission',fontsize=20)

    plt.savefig("transmission/varPc_L=%.1fcm_Pp=%.1fmW.png" % (L*100,Pp*1e6))
    plt.show()



def plot_T_varPp():

    filename = "var_Pp"


    Ppratios = np.logspace(-5,4,ns)


    Tpco = []; Tpcn = []
    Tcco = []; Tccn = []
    
    try:
        with open('npy/Tpco%s.npy' % filename, 'rb') as f: Tpco=np.load(f)
        with open('npy/Tpcn%s.npy' % filename, 'rb') as f: Tpcn=np.load(f)
        with open('npy/Tcco%s.npy' % filename, 'rb') as f: Tcco=np.load(f)
        with open('npy/Tccn%s.npy' % filename, 'rb') as f: Tccn=np.load(f)
    except (FileNotFoundError,NameError) as e:
        for rat in tqdm(Ppratios):
            
            Ppi = Pc*rat
            E0pi = (2*Ppi/(c*ε0*Ap*n))**0.5 # Amplitude of probe field (V/m)
            Ω0pi = -dip1*E0pi/hbar # Probe Rabi frequency (rad Hz)

            temp = co_T(Ω0pi,Ω0c,E0pi,E0c)
            Tpco.append(temp[0])
            Tcco.append(temp[1])
            temp = counter_T(Ω0pi,Ω0c,E0pi,E0c)
            Tpcn.append(temp[0])
            Tccn.append(temp[1]) 


        with open('npy/Tpco%s.npy' % filename, 'wb') as f: np.save(f,Tpco)
        with open('npy/Tpcn%s.npy' % filename, 'wb') as f: np.save(f,Tpcn)
        with open('npy/Tcco%s.npy' % filename, 'wb') as f: np.save(f,Tcco)
        with open('npy/Tccn%s.npy' % filename, 'wb') as f: np.save(f,Tccn)




    plt.ticklabel_format(scilimits=(0,0))
    plt.plot(Ppratios,Tpco,c='orangered', linewidth=7,label=r'Probe: Co')
    plt.plot(Ppratios,Tpcn,c='steelblue', linewidth=3,label=r'Probe: Counter')
    plt.plot(Ppratios,Tcco,c='darkorchid',linewidth=7,label=r'Control: Co')
    plt.plot(Ppratios,Tcco,c='forestgreen',linewidth=3,label=r'Control: Counter')
    plt.grid()
    plt.legend(fontsize=17)
    plt.title(r'Transmission, varying $P_p$: $L=%.1f$ cm, $P_c=%.1f$ mW' % (L*100,Pc*1e3),fontsize=20)
    plt.xlabel(r'$P_p/P_c$',fontsize=20)
    plt.xscale('log')
    plt.ylabel(r'Transmission',fontsize=20)

    plt.savefig("transmission/varPp_L=%.1fcm_Pc=%.1fmW.png" % (L*100,Pc*1e3))
    plt.show()



plot_T_varPc()
# plot_T_varPp()
