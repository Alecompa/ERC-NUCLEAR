import numpy as np
import ROOT
import modules.analysis as analysis


import matplotlib.pyplot as plt

stop_power = np.loadtxt("data/stopping/helium_in_boron.dat")
cross_section = np.loadtxt("data/cross_sections/10b_an.cross")

#print( stop_power )
#print( cross_section )

energies = np.linspace( 300, 400, 100 ) 
delta_e = 33 # keV 

yields = []
for energy in energies:
    yields.append( analysis.calculate_yield( energy, delta_e, cross_section, stop_power ) )

#print( yields )

yield_array = np.array( yields )
#set log y
plt.yscale('log')
plt.plot( energies, yield_array )

plt.show()