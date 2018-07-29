""" outputting pqr file for abps 
using amber radii

Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius
"""

import os,sys,time,copy,numpy,math
import atoms,forceField,amberFF,defaults

# make output lines:
def returnLine(atom, ff):
    """make a bgf line from atom"""
    if atom.aName.strip() == '':
        aName = '?'
    else:
        aName = atom.aName

    if atom.rName.strip() == '':
        rName = '?'
    else:
        rName = atom.rName

    if atom.chain.strip() == '':
        chain = 'c'
    elif atom.chain.isdigit():
        chain = 'c'
    else:
        chain = atom.chain
    
    if atom.ffType in ff.atomTypes:
        radius = ff.atomTypes[atom.ffType].vdwR / 2.0
    else:
        radius = -9999.9999
        print 'Error when making pqr: atom type %s is not in the forcefield' % atom.ffType

    line = '%s %d %s %s %s %d %.5f %.5f %.5f %.5f %.5f \n'
    line %= (atom.aTag, atom.aNo, aName, rName, chain, atom.rNo, atom.xyz[0], \
             atom.xyz[1], atom.xyz[2], atom.charge, radius)
    return line

def writeLines(AL, ff):
    """returns list of lines corresponding to the atom list"""
    lines = []
    getlogin = os.getenv('LOGNAME')
    getenv = os.getenv('HOSTNAME')
    gettime = time.strftime('%X %x %Z')
    lines.append('REMARK saved by %s@%s on %s\n' %(getlogin,getenv,gettime))
    for a in AL: lines.append(returnLine(a, ff))
    return lines


def write(outfile,AL,amberRadii=True):
    """uses amber radii for the radius field of the pqr file
    charges are taken from the current structure

    if amberRadii == False will use Dreiding vdw radii (using Dreiding does not work yet)
    """
   
    ff = forceField.ForceField()
    if amberRadii:
        amberFF.readAmberDat(ff,fileName=defaults.amberDat)
        amberFF.readAmberDat(ff,fileName=defaults.amberDatGaff)  
    else:
        shortpar.read(ff)
    forceField.fixAtomLabels(AL,ff)  

    fout = open(outfile,'w')
    fout.writelines(writeLines(AL, ff))
    fout.close()

def insertXyz(origPqrFile,newPqrFile,xyz):
    """inserts new coordinates xyz into the given pqr file"""
    i = 0
    lines = open(origPqrFile).readlines()
    with open(newPqrFile,'w') as f:
        for line in lines:
            s = line.split()
            if s[0] in ['ATOM','HETATM']:
                newline = '%s %s %s %s %s %s %.5f %.5f %.5f %s %s \n'
                newline %= (s[0],s[1],s[2],s[3],s[4],s[5],xyz[i,0],xyz[i,1],xyz[i,2],s[9],s[10])
                f.write(newline)
                i += 1
            else:
                f.write(line)


