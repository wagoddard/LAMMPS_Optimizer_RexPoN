import os
import sys
import glob, re
from sys import argv


def readfile(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    energies = []
    header = ' '
    for line in lines:
        if not line.startswith('#') and len(line.split()) > 1:
            energies.append([])
            for i in range(len(line.split())):
                energies[-1].append(line.split()[i])
        elif line.startswith('#str'):
            header = line
    return energies,header

def writefile(etable,train,header):
    g = open('eng.table','w')
    g.write('%-30s Eqm %s'%('#str',header.split('#str')[1]))
    geos = {}
    nelem =  len(header.split())
    for i in range(len(etable)):
      geos[etable[i][0]] = etable[i][1:]
      temp = len(etable[i][1:])
      if temp+1 != nelem:  
        print '%s does not have correct number of cols '%etable[i][0]
        sys.exit()
    for i in range(len(train)):
        weight = float(train[i][0])
        struct = train[i][2::2]
        sign   = train[i][1:-1:2]
        tr_err =  float(train[i][-1])
        eng_line = ' '
        for k in range(nelem-1):
          eng_comp = 0.0
          for j in range(len(struct)):
            factor = -1 if sign[j] == '-' else 1
            temp = struct[j].split('/')[0]
            eng_comp += float(geos[temp][k]) * factor 
          eng_line = eng_line + ' ' + str('%0.4f'%eng_comp) + ' '
        name  = struct[0].split('/')[0]
        g.write('%-30s %0.4f %s\n'%(name,tr_err,eng_line))

if __name__ == "__main__":

    if len(argv) != 3:
        print "Usage: python analysis.py [energ_table] [trainset]"
        print " column: the number of column that xyz coordinates start"
        sys.exit()

    table = argv[1]
    trainset =  argv[2]

    etable,header = readfile(table)
    etrain,temp   = readfile(trainset)
    writefile(etable,etrain,header)  

 
