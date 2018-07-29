#!usr/bin/env python 

"""
Convert xyz file to lammps_data file for RexPoN force field 
Usage: python xyz2rexpon.py input.xyz ffield output.data
       then insert the box info in the next line:
       a b c alpha beta gamma
"""

import os, sys
from sys import argv
from math import sin, cos, sqrt, radians

if len(argv) != 4:
	print "Usage: python xyz2rexpon.py input.xyz ffield output.data"
	sys.exit()

fin = open(os.path.join(os.getcwd(), argv[1]),'r')
ffield = open(os.path.join(os.getcwd(), argv[2]),'r')
fout = open(os.path.join(os.getcwd(), argv[3]),'w')
#element = {}

mass = {'H':1.00794,'He':4.002602,'Li':6.941,  'Be':9.012182,\
        'B':10.811,  'C':12.0107, 'N':14.00674,'O':15.9994, \
        'F':18.9984032,'Ne':20.1797,'Na':22.98977,'Mg':24.3050,\
        'Al':26.981538,'Si':28.0855,'P':30.973761,'S':32.066, \
        'Cl':35.4527,'Ar':39.948,'K':39.0983,'Ca':40.078, \
	'Pd':106.42,'Mo':95.94, 'Cu':63.546, 'Ga':69.723, 'Ge':72.64, \
        'Sn':118.71,'Pb':207.2,'W':183.84,'Hw':1.00794,'Ow':15.9994, 'Fe':55.845,\
        'Ns':14.00674,'Nf':14.00674 }

def abc2vector_alongx(a, b, c, alpha, beta, gamma):
	"""
	Convter the a/b/c/alpha/beta/gamma to lattice vector of a long x and b in xy plane.
	Calculated from lammps manual
	by Saber Naserifar
	"""

	r_alpha = radians(alpha)
	r_beta = radians(beta)
	r_gamma = radians(gamma)

	lx = a
	ly = b*sin(r_gamma)
	lz = c*sqrt(
				sin(r_beta)**2 - \
				(cos(r_alpha)-cos(r_beta)*cos(r_gamma))**2/(sin(r_gamma)**2)
				)
	xy = b*cos(r_gamma)
	xz = c*cos(r_beta)
	yz = c*(cos(r_alpha)-cos(r_beta)*cos(r_gamma))/sqrt(1-cos(r_gamma)**2)

	return lx, ly, lz, xy, xz, yz

def readxyz(file):
	list = []
	file.next()
	file.next()
	for i in file:
		list.append(i.split())
	return list

def readffield(file):
        atomdict = {}
        cnt = 0 
        lines = file.readlines()
        ntypes = int(lines[41].split()[0])
        line = 45 
        for i in range(ntypes):
            atom = lines[line].split()[0]
            cnt += 1 
            atomdict[atom] = cnt 
            line += 5  
        return atomdict


#Remove the duplicate element in the list
def thinList(list):
    for i in list:
        num = list.count(i)
        if num > 1:
            t=1
            while t < num:
                list.remove(i)
                t=t+1
    list.sort()
    return list

xyzlist = readxyz(fin)
atomdict = readffield(ffield)

#atomlist = thinList([i[0] for i in xyzlist])

a, b, c, alpha, beta, gamma = raw_input("Input the cell info: a b c alpha beta gamma\n").split()
lx, ly, lz, xy, xz, yz = abc2vector_alongx(float(a), float(b), float(c), \
										float(alpha), float(beta), float(gamma))
#print lx, ly, lz, xy, xz, yz

fout.write("# Data file generated from xyz2data_pqeq.py by Saber\n\n")
fout.write("%-7s atoms\n" %len(xyzlist))
fout.write("%-2s atom types\n\n" %len(atomdict))
fout.write("0.0000 %.4f xlo xhi\n" %lx)
fout.write("0.0000 %.4f ylo yhi\n" %ly)
fout.write("0.0000 %.4f zlo zhi\n\n" %lz)
fout.write("%.4f %.4f %.4f xy xz yz\n\nMasses\n\n" %(xy, xz, yz))
for key,value in sorted(atomdict.iteritems(), key=lambda (k,v): (v,k)):
        fout.write("%-2s%7f   # %-2s\n" %(value,mass[key],key))
fout.write("\nAtoms\n\n")
for i in xyzlist:
	fout.write("%s 1  %i %f %f %f %f %f %f %f\n" \
		    %(xyzlist.index(i)+1,atomdict[i[0]],0.0\
		    ,float(i[1]),float(i[2]),float(i[3]),0.0,0.0,0.0))
fout.write("\n\n")

fin.close()
fout.close()

