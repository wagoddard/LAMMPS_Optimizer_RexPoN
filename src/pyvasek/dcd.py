""" reading dcd file
"""

import sys,numpy
from packages.pyvmd import VMD, Molecule

import periodic

def read(bgffile,dcdfile):
    m = Molecule.Molecule()
    m.load(bgffile)
    m.load(dcdfile)
    m.delFrame(first=0, last=0)
    return m

def getFrame(mol, frameid):
    '''isnumpy array of float32'''
    return VMD.vmdnumpy.timestep(int(mol),frameid)

def getFrames(mol):
    frames = []    
    for i in range(mol.numFrames()):
        frames.append( getFrame(mol,i) )
    return numpy.array( frames )

def setFrame(mol, frameid, frame):
    '''isnumpy array of float32'''
    y = numpy.array(frame,dtype=numpy.float32)
    x = VMD.vmdnumpy.timestep(int(mol), frameid)
    x[:] = y

def setFrames(mol, frames):
    for i in range(len(frames)):
        setFrame(mol,i,frames[i])

