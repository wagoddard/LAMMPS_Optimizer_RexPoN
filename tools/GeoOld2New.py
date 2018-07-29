# Parse a geo file and use tpascal script to convert to data_lammps_ files

import os
import sys
import glob
from sys import argv


def readgeo(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    # Parse the geo file by BIOGRF
    out_array = []
    for line in lines:
        if line.startswith("BIOGRF"):
            out_array.append(line)
        else:
            out_array[-1] += line

    print '(2) Writing data files...'
    for i in range(len(out_array)):
        s = out_array[i]
        name = s[s.find('DESCRP')+len('DESCRP'):].split()[0]
        outbgf = '%s.tmp.bgf'%name
        h = open(outbgf, 'w')
        h.write(s)
        h.close()

    return

def writexyz(col):
    files = glob.glob('*.tmp.bgf')
    for file in files:
      name = '.'.join(file.split('.')[:-2])+'.xyz'
      f = open(file, 'r')
      g = open(name, 'w')
      lines = f.readlines()
      natoms = 0
      for line in lines:
        if line.startswith('HETATM') or line.startswith('ATOM'):
          natoms += 1
      g.write('%i \n\n'%natoms)
      for line in lines:
        if line.startswith('HETATM') or line.startswith('ATOM'): 
          s = line.split()
          type = s[2]; x = s[col-1]; y = s[col]; z = s[col+1]   
          g.write('%s %s %s %s\n'%(type,x,y,z))
    os.system('rm *.tmp.bgf')
    g.close()
    f.close()
    return 

def newgeo():
    files = glob.glob('*.xyz')
    f = open('newgeo','w')
    for file in files:
      name = '.'.join(file.split('.')[:-1])
      os.system('jaguar babel -ixyz %s.xyz -obgf %s.bgf'%(name,name)) 
    files = glob.glob('*.bgf')
    files.sort()
    for file in files:
      lines = open(file,'r').readlines()
      for line in lines:
        if line.startswith('DESCRP'):
          geoname = '.'.join(file.split('.')[:-1])
          f.write('DESCRP %s\n'%geoname)
        else:
          f.write(line)
    os.system('rm *.bgf *.xyz')


if __name__ == "__main__":

    if len(argv) != 3:
        print "Usage: python GeoOld2New.py [geo] [column]"
        print " column: the number of column that xyz coordinates start"
        sys.exit()

    ifile = argv[1]
    col = int(argv[2])


    readgeo(ifile)
    writexyz(col)
    newgeo()
