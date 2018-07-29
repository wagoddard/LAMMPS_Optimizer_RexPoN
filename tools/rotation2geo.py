# Two command line inputs: name of structure, number of structures

import os
import sys
import glob

helpstring = """ USAGE: python rotation2geo.py [file.01.mae]
                      file.01.mae: output 01.mae file from jaguar rigid rotation scan
                      converts to .mol2, then to .bgf, and finally to geo
             """

if len(sys.argv) != 2:
    print helpstring
    sys.exit()

name = sys.argv[1]
#num  = int(sys.argv[2])

os.system('/$SCHRODINGER/utilities/structconvert -imae %s.01.mae -omol2 %s.mol2' % (name,name))
print 'Converted %s.01.mae to %s.mol2' % (name,name)
os.system('jaguar babel -imol2 %s.mol2 -obgf %s.bgf -m' % (name, name))
print 'Converted %s.mol2 to %s.bgf[s]' % (name,name)

num = len(glob.glob("%s*.bgf" % name))
print 'Reading %d files' % num

for i in range(num):
    f = open("%s%d.bgf" % (name,i+1), 'r')
    lines = f.readlines()
    f.close()

    f = open("%s%d.bgf" % (name,i+1), 'w')
    for line in lines:
        if line.startswith('DESCRP'):
            f.write("DESCRP %s%d\n" % (name,i+1))
        else:
            f.write(line)

    f.close()

lst = "cat "
lst2 = "rm "
for i in range(num):
    lst += "%s%d.bgf " % (name,i+1)
    lst2 += "%s%d.bgf " % (name,i+1)
lst += "> geo-%s" % name

print lst
os.system(lst)

print lst2
os.system(lst2)

