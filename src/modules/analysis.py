import numpy as np
import scipy.integrate as integrate

N_A = 6.022e23  # Avogadro number - atoms/mol
e = 1.602176634e-19  # elementary charge - C


ug_to_gram = 1e-6
barn_to_cm2 = 1e-24


def thickess_atoms (thickness_ug, n_atoms, mol_weight):
        # thickness is in ug/cm2, n_atoms is number of active atoms per molecule, mol_weight is the mass number of the molecule
        # thick_m is in atoms/cm2 
    thick_m = n_atoms * thickness_ug * ug_to_gram * N_A / mol_weight
    return thick_m


# Infinitesimal of the yield integral
def infinitesimal( energy, cross_section, stop_power ):
    # cross_graph wants energies in MeV, stop_graph in keV
    cross = np.interp( energy / 1e3, cross_section[:,0], cross_section[:,1] )
    stop = np.interp( energy, stop_power[:,0], stop_power[:,1] )
    return barn_to_cm2 * cross / stop * 1e15 * 1e3 # 1e15 atoms/cm^2, 1e3 eV

# Calculate the yield
def calculate_yield( energy, delta_e, cross_section, stop_power ):
    # Integrate the infinitesimal between [energy, energy - delta_e]
    return integrate.quad( lambda energy: infinitesimal( energy, cross_section, stop_power ), energy - delta_e, energy )[0]
        #infinitesimal, energy - delta_e, energy )[0]
