#from decimal import *
import os
import sys
import glob

helpstring = """ Convert .grd file to trainset. Input structure and reference. """

#if len(sys.argv) != 3:
#    print helpstring
#    sys.exit()

#name = sys.argv[1]
#num  = sys.argv[2]

dirs = glob.glob("*.grd")

for name in dirs:
    f = open("%s" % name, 'r')
    lines = f.readlines()
    f.close()

    name = name[:-4]
    
    print name
    num = 0
    j = 1
    min = float(1000000.0)
    for i in range(len(lines)):
        if len(lines[i].split()) == 1 and float(lines[i].split()[0]) < min:
            num = j
            min = float(lines[i].split()[0])
            j += 1
            print float(lines[i].split()[0]), j-1
        elif len(lines[i].split()) == 1:
            j += 1
 
    j = 1
    towrite = "ENERGY\n"
    for i in range(len(lines)):
        if len(lines[i].split()) == 1:
            towrite += "1.0 + %s%d/1 - %s%d/1 %s\n" % (name, j, name, num, str(float(lines[i].split()[0])/float(4.184)))
            j += 1
    towrite += "ENDENERGY"

    f = open("trainset-%s" % name, 'w')
    f.write(towrite)
    f.close()
