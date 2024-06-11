# Importing all parameters
from parameters import *
from transmission import co_T,counter_T # Importing functions transmission.py
from tqdm import tqdm

# Number of probe/control power ratios to be plotted
ns=200 # Overhead is polynomial

# Plotting probe and control transmissions, varying the control power

def plot_T_varPc():

    filename = "var_Pc"

    Pratios = np.logspace(-5,0.1,ns) # Power ratios are generated in log10 space


    Tpco = []; Tpcn = [] # Probe transmissions for co- and counter-propagation
    Tcco = []; Tccn = [] # Control transmissions for co- and counter-propagation

    # If npy files exist with the same names, then load them to avoid unnecessarily
    # determining transmissions with the same parameters

    # If the user wants to determine transmissions with a different set of parameters,
    # but files already exist, simply delete them 
    try:
        with open('npy/Tpco%s.npy' % filename, 'rb') as f: Tpco=np.load(f)
        with open('npy/Tpcn%s.npy' % filename, 'rb') as f: Tpcn=np.load(f)
        with open('npy/Tcco%s.npy' % filename, 'rb') as f: Tcco=np.load(f)
        with open('npy/Tccn%s.npy' % filename, 'rb') as f: Tccn=np.load(f)
    except (FileNotFoundError,NameError) as e:
        
        # If no such files exist, then determine transmissions

        for rat in tqdm(Pratios):
            
            Pci = Pp/rat
            E0ci = (2*Pci/(c*ε0*Ac*n))**0.5 # Amplitude of control field (V/m)
            Ω0ci = -dip2*E0ci/hbar # Control Rabi frequency (rad Hz)

            # Appending to transmission lists
            temp = co_T(Ω0p,Ω0ci,E0p,E0ci)
            Tpco.append(temp[0])
            Tcco.append(temp[1])
            temp = counter_T(Ω0p,Ω0ci,E0p,E0ci)
            Tpcn.append(temp[0])
            Tccn.append(temp[1]) 

        # Saving the transmission lists into npy files in npy folder        
        with open('npy/Tpco%s.npy' % filename, 'wb') as f: np.save(f,Tpco)
        with open('npy/Tpcn%s.npy' % filename, 'wb') as f: np.save(f,Tpcn)
        with open('npy/Tcco%s.npy' % filename, 'wb') as f: np.save(f,Tcco)
        with open('npy/Tccn%s.npy' % filename, 'wb') as f: np.save(f,Tccn)

    # LaTeX symbols are used for labels

    plt.ticklabel_format(scilimits=(0,0))
    # The colours below were chosen in consideration of colourblindness
    plt.plot(Pratios,Tpco,c='orangered', linewidth=7,label=r'Probe: Co')
    plt.plot(Pratios,Tpcn,c='steelblue', linewidth=3,label=r'Probe: Counter')
    plt.plot(Pratios,Tcco,c='darkorchid',linewidth=7,label=r'Control: Co')
    plt.plot(Pratios,Tcco,c='forestgreen',linewidth=3,label=r'Control: Counter')
    plt.grid()
    plt.legend(fontsize=17)
    plt.title(r'Transmission, varying $P_c$: $L=%.1f$ cm, $P_p=%.1f$ \textmu W' % (L*100,Pp*1e6),fontsize=20)
    plt.xlabel(r'$P_p/P_c$',fontsize=20)
    plt.xscale('log')
    plt.ylabel(r'Transmission',fontsize=20)

    # Plot is saved in transmission folder
    plt.savefig("transmission/varPc_L=%.1fcm_Pp=%.1fmW.png" % (L*100,Pp*1e6))
    plt.show()

# Plotting probe and control transmissions, varying the probe power

def plot_T_varPp():

    filename = "var_Pp"


    Pratios = np.logspace(-5,4,ns) # Power ratios are generated in log10 space


    Tpco = []; Tpcn = [] # Probe transmissions for co- and counter-propagation
    Tcco = []; Tccn = [] # Control transmissions for co- and counter-propagation

    # If npy files exist with the same names, then load them to avoid unnecessarily
    # determining transmissions with the same parameters

    # If the user wants to determine transmissions with a different set of parameters,
    # but files already exist, simply delete them   
    try:
        with open('npy/Tpco%s.npy' % filename, 'rb') as f: Tpco=np.load(f)
        with open('npy/Tpcn%s.npy' % filename, 'rb') as f: Tpcn=np.load(f)
        with open('npy/Tcco%s.npy' % filename, 'rb') as f: Tcco=np.load(f)
        with open('npy/Tccn%s.npy' % filename, 'rb') as f: Tccn=np.load(f)
    except (FileNotFoundError,NameError) as e:
        
        # If no such files exist, then determine transmissions

        for rat in tqdm(Pratios):
            
            Ppi = Pc*rat
            E0pi = (2*Ppi/(c*ε0*Ap*n))**0.5 # Amplitude of probe field (V/m)
            Ω0pi = -dip1*E0pi/hbar # Probe Rabi frequency (rad Hz)

            # Appending to transmission lists
            temp = co_T(Ω0pi,Ω0c,E0pi,E0c)
            Tpco.append(temp[0])
            Tcco.append(temp[1])
            temp = counter_T(Ω0pi,Ω0c,E0pi,E0c)
            Tpcn.append(temp[0])
            Tccn.append(temp[1]) 

        # Saving the transmission lists into npy files in npy folder    
        with open('npy/Tpco%s.npy' % filename, 'wb') as f: np.save(f,Tpco)
        with open('npy/Tpcn%s.npy' % filename, 'wb') as f: np.save(f,Tpcn)
        with open('npy/Tcco%s.npy' % filename, 'wb') as f: np.save(f,Tcco)
        with open('npy/Tccn%s.npy' % filename, 'wb') as f: np.save(f,Tccn)


    # LaTeX symbols are used for labels

    plt.ticklabel_format(scilimits=(0,0))
    # The colours below were chosen in consideration of colourblindness
    plt.plot(Pratios,Tpco,c='orangered', linewidth=7,label=r'Probe: Co')
    plt.plot(Pratios,Tpcn,c='steelblue', linewidth=3,label=r'Probe: Counter')
    plt.plot(Pratios,Tcco,c='darkorchid',linewidth=7,label=r'Control: Co')
    plt.plot(Pratios,Tcco,c='forestgreen',linewidth=3,label=r'Control: Counter')
    plt.grid()
    plt.legend(fontsize=17)
    plt.title(r'Transmission, varying $P_p$: $L=%.1f$ cm, $P_c=%.1f$ mW' % (L*100,Pc*1e3),fontsize=20)
    plt.xlabel(r'$P_p/P_c$',fontsize=20)
    plt.xscale('log')
    plt.ylabel(r'Transmission',fontsize=20)

    # Plot is saved in transmission folder
    plt.savefig("transmission/varPp_L=%.1fcm_Pc=%.1fmW.png" % (L*100,Pc*1e3))
    plt.show()


# The user can call the functions above to plot the transmission as a function of
# probe/control power ratio, varying probe or control power

plot_T_varPc()
# plot_T_varPp()
