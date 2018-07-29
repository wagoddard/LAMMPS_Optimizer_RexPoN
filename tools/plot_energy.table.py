# Input: history

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def main(input):
    dat = np.genfromtxt(input, delimiter=None, skip_header=1)      # Read from CSV (fort.99)

    name_dat    = dat[:,0]
    etotal_dat  = dat[:,1]
    pe_dat      = dat[:,2]
    ke_dat      = dat[:,3]
    evdw_dat    = dat[:,4]
    electro_dat = dat[:,5]
    ebond_dat   = dat[:,6]
    eangle_dat  = dat[:,7]
    edihed_dat  = dat[:,8]

    # Convert pandas to np
    name, etotal, pe, ke, evdw, electro, ebond, eangle, edihed = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
    for i in range(len(name_dat)):
        name    = np.append(name, name_dat[i])
        etotal  = np.append(etotal , etotal_dat[i])
        pe      = np.append(pe, pe_dat[i])
        ke      = np.append(ke, ke_dat[i])
        evdw    = np.append(evdw, evdw_dat[i])
        electro = np.append(electro, electro_dat[i])
        ebond   = np.append(ebond, ebond_dat[i])
        eangle  = np.append(eangle, eangle_dat[i])
        edihed  = np.append(edihed, edihed_dat[i])


    # TODO: Figure out what the X-axis should be
    X = np.array([0, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 10, 20, 30, 40, 50, 60, 70, 80])

    # Make scatter plots
    plt.scatter(X, etotal,s=20,  c='red')
    plt.scatter(X, evdw,s=20,    c='green')
    plt.scatter(X, electro,s=20, c='blue')
    plt.scatter(X, ebond, s=20,  c='orange')
    plt.scatter(X, eangle,s=20,  c='purple')
    plt.scatter(X, edihed, s=20, c='yellow')

    plt.legend(['etot', 'evdw', 'electr', 'ebond', 'eang', 'edih'])

    plt.xlabel('Degrees')
    plt.ylabel('Energy (kcal/mol)')

    plt.xlim(0, 180)

    plt.show()

if __name__ == "__main__":

    helpstring = """ Input: energy_table.txt Output: *png """

    if len(sys.argv) != 2:
        print helpstring
        sys.exit()
  
    input = sys.argv[1]

    main(input)

