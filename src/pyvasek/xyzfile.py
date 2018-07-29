""" reading and outputting xyz file
"""

import sys,numpy
import atoms,bgf

def write(filename,xyzlist,names=[],al=[]):
    ''' names are names of the frames, al is the list of atoms'''
    nframes = len(xyzlist)
    if nframes == 0:
        print >> sys.stderr, 'Error: number of frames to store in xyz format is zero.'
        sys.exit(1)
        
    if len(names) == 0:
        names = ['frame%d' % i for i in range(nframes)]
    elif len(names) != nframes:
        print >> sys.stderr, 'Frame names have incorrect length %d, should be %d' % (len(names),nframes)
        sys.exit(1)

    natoms = len(xyzlist[0])
    if len(al) == 0:
        atomnames = ['%d' % i for i in range(natoms)]
    elif len(al) != natoms:
        print >> sys.stderr, 'Error: atom list al has length %d but the frame has %d atoms' % (len(al),natoms)
        sys.exit(1)
    else:
        atomnames = [a.aName for a in al]
        
    with open(filename,'w') as f:
        for iframe in range(nframes):
            f.write('%d\n' % natoms)
            f.write('%s\n' % names[iframe])
            for iatom in range(natoms):
                f.write('%s %.10f %.10f %.10f\n' % (atomnames[iatom],xyzlist[iframe][iatom,0],xyzlist[iframe][iatom,1],xyzlist[iframe][iatom,2]))

def read(filename):
    """ read frames with coordinates from xyz file: returns list of coordinate xyz (numpy array) and list of strings with frame names"""
    xyzlist = []
    namelist = []
    with open(filename) as f:
        natoms = int(f.readline())
        while natoms > 0:
            namelist.append( f.readline().strip() ) # comment line
            xyz = numpy.zeros((natoms,3)) 
            for i in range(natoms):
                s = f.readline().split()
                xyz[i,0] = float(s[1])
                xyz[i,1] = float(s[2])
                xyz[i,2] = float(s[3])
            xyzlist.append(xyz)
            try:
                natoms = int(f.readline())
            except ValueError:
                natoms = 0
    return numpy.array(xyzlist), namelist

def readAtomList(filename):
    """ read single frame from xyz file """
    al = []
    with open(filename) as f:
        natoms = int(f.readline())
        f.readline() # comment line
        for i in range(natoms):
            atom = atoms.Atom()
            s = f.readline().split()
            atom.aName = s[0]
            atom.xyz[0] = float(s[1])
            atom.xyz[1] = float(s[2])
            atom.xyz[2] = float(s[3])
            al.append(atom)
    return al

def writeAtomList(filename,al):
    xyz = atoms.getXyzFromList(al)
    write(filename,numpy.array([xyz]),names=[filename],al=al)    
    

def catbgf2xyz(catbgfFile,xyzFile,firstLabel = 1,keepFrames = []):
    xyzlist = []
    names = []
    with open(catbgfFile) as f:
        lines = []
        counter = firstLabel
        
        line = f.readline()
        while line != '':
            lines.append(line)
            if line.startswith('END'):
                al = bgf.readLines(lines)
                if (keepFrames == []) or (counter in keepFrames):
                    xyz = atoms.getXyzFromList(al)
                    xyzlist.append(xyz)
                    names.append('%d'% counter)
                    print counter
                counter += 1
                lines = []
            line = f.readline()
    write(xyzFile,numpy.array(xyzlist),names=names,al=al)
    

