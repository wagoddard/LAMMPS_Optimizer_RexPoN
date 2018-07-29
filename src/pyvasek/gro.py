""" reading and outputting gro file

    title string (free format string, optional time in ps after 't=')
    number of atoms (free format integer)
    one line for each atom (fixed format, see below)
    box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), the last 6 values may be omitted (they will be set to zero). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0. 


    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places) 

"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"

"""

import os,sys,time,copy,numpy,math
import atoms

# population functions for one atom
def readLine(line):
    """store bgf line info into atom and return the atom"""
    a = atoms.Atom()

    a.rNo   = int(line[0:5])
    a.rName = line[5:10].strip()
    a.aName = line[10:15].strip()
    a.aNo   = int(line[15:20])
    a.xyz   = 10.0*numpy.array( [float(line[20:28]), float(line[28:36]), float(line[36:44]) ])
    a.velocities = numpy.array( [float(line[44:52]), float(line[52:60]), float(line[60:68]) ])
    return a

# make output lines:
def returnLine(atom):
    """make a bgf line from atom"""
    format = "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"
    aName = atom.aName.strip()
    rName = atom.rName.strip()
    if hasattr(atom, 'velocities'):
        vel = atom.velocities
    else:
        vel = numpy.array([0.,0.,0.])
    xyzNM = atom.xyz / 10.0
    line = format % (atom.rNo,rName,aName,atom.aNo,xyzNM[0],xyzNM[1],xyzNM[2],vel[0],vel[1],vel[2])
    return line

def readLines(lines):
    """reads the bgf which is passes as list of lines, returns list of atoms, discard advanced features of bgf"""
    al = []
    n = int(lines[1])
    for i in range(n):
        al.append(readLine(lines[2+i]))
    boxline = lines[2+n].split()
    boxGro = 10.0*numpy.array([float(boxline[0]),float(boxline[1]),float(boxline[2])])
    return al,boxGro

def writeLines(AL,boxGro):
    """returns list of lines corresponding to the atom list"""
    lines = []
    getlogin = os.getenv('LOGNAME')
    getenv = os.getenv('HOSTNAME')
    gettime = time.strftime('%X %x %Z')
    lines.append('REMARK saved by %s@%s on %s\n' %(getlogin,getenv,gettime))
    if len(AL) > 99999:
        sys.exit('Error in writing .gro file: There are too many atoms: %d, max is 99999'%len(AL))
    lines.append('%5d\n'%len(AL))
    atoms.renumberAtomsFrom(AL)
    for a in AL: lines.append(returnLine(a))
    boxNM = boxGro / 10.0
    lines.append('%10.5f%10.5f%10.5f\n'%(boxNM[0],boxNM[1],boxNM[2]))
    return lines

def read(file, debug=False):
    AL,boxGro = readLines(open(file).readlines())
    if debug: print '%d atoms read from %s.' %(len(AL), file)
    return AL,boxGro

def write(outfile,AL,boxGro,debug=False):
    '''boxGro are 3 numbers with the size of the rectangular box: box_x, box_y, box_z'''
    fout = open(outfile,'w')
    fout.writelines(writeLines(AL,boxGro))
    fout.close()
    if debug: print 'File %s written.' % outfile
      

