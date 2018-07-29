# Input: history

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def main(input, name):
    f = open(input, 'r')
    lines = f.readlines()
    f.close()

    file = []
    reorder = []
    for line in lines:
        if line[0] != '#' and line.split()[3].split('.')[0] not in file:
            file.append(line.split()[3].split('.')[0])
            reorder.append([])
            reorder[-1].append(line)
        elif line[0] != '#':
            reorder[-1].append(line)

    for c,a in enumerate(reorder):
        x   = []
        lmp = []
        qm  = []

        for b in a:
            x.append(float(b.split()[3].split('-')[-1].split('dp')[0]))
            lmp.append(float(b.split()[6]))
            qm.append(float(b.split()[7]))

        plt.plot(x, lmp, '.-r', label='PQEq')
        plt.plot(x, qm, '.-b', label='QM')
        #plt.show()
        plt.xlabel('Distance (A)')
        plt.ylabel('Energy (kcal/mol)')
        plt.legend(['LAMMPS','QM'])
        #plt.xlim(0, 180)
        #plt.show()
        plt.savefig('%s-%s.png' % (name, file[c]))
        plt.clf()

    return

    dat = np.genfromtxt(input, delimiter=None, skip_header=1)      # Read from CSV (fort.99)
    lmp_dat = dat[:,6]
    qm_dat  = dat[:,7]
    lmp, qm = np.array([]), np.array([])                         # Convert pandas to np
    for i in range(len(lmp_dat)):
        lmp = np.append(lmp, lmp_dat[i])
        qm  = np.append(qm , qm_dat[i])
    X = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180])
    #X = np.array([80,90,100,110,120,130,140])
    #X = np.array([88,92,96,100,104,108,112,116,120,124,128,132,136,140])
    # Make a scatter plot
    #LMP = plt.scatter(X, lmp, s=20, c='red', label='LAMMPS')
    #QM  = plt.scatter(X,  qm, s=20, c='blue', label='QM')
    #LMP = plt.plot(X, lmp, '.-r', label='LAMMPS')
    #QM = plt.plot(X, qm, '.-b', label='QM')
    #plt.legend(['LAMMPS','QM'])
    #plt.xlabel('Degrees')
    #plt.ylabel('Energy (kcal/mol)')
    #plt.xlim(0, 180)
    #plt.show()
    #plt.savefig('%s.png' % name)

if __name__ == "__main__":

    helpstring = """ Input: fort.99, name Output: *png """

    if len(sys.argv) != 3:
        print helpstring
        sys.exit()
  
    input = sys.argv[1]
    name = sys.argv[2]

    main(input,name)


