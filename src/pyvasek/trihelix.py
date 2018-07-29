#!/usr/bin/env python

import sys, string,copy, numpy, time
from pyvasek import *


def getChain(al,chain):
    alchain = []
    indices = []
    for i,a in enumerate(al):
        if a.chain == chain:
            alchain.append(a)
            indices.append(i)
    return alchain, indices

def getBackbone(al):
    backbone = ['C','N','CA']
    albone = []
    for a in al:
        if a.aName.strip() in backbone:
            albone.append(a)
    return albone 

def tmAxis(al,v):
    xyz = atoms.getXyzFromList(getBackbone(al))
    center = numpy.mean(xyz, 0)
    xyzc = xyz - center
    inertia = -numpy.dot(xyzc.transpose(), xyzc)+numpy.diag((1,1,1))*sum(sum(numpy.square(xyzc)))-numpy.diag((1,1,1))
    eigenvalues, eigenvectors = numpy.linalg.eig(inertia)
    eigenvalue = min(eigenvalues)
    eigenvector = eigenvectors[:,list(eigenvalues).index(eigenvalue)]
    print eigenvectors
    print eigenvalues
    plotAxis(center, eigenvector, v)
    return center,eigenvector 

def plotAxis(center, axis, v):
    a1 = atoms.Atom()
    a2 = atoms.Atom()
    a1.aName = 'K'
    a2.aName = 'K'
    a1.makeBond(a2)
    a1.xyz = center + 30.* axis
    a2.xyz = center - 30.* axis
    v.show([a1,a2])    

def run(protein,plot=None):
    v = vmd.vmd()
    al = bgf.readFile(protein)
    en = energy.Energy(debug = False)
    en.load(al)

    chain = '3'
    tm,indices = getChain(al,chain)
    
    center, axis = tmAxis(tm, v)

    xyzal = atoms.getXyzFromList(al)
    xyztm = xyzal[indices,:]
    xyztmcopy = xyztm.copy()
    v.show(al)
    v.pvmd.__call__('mol modselect 0 top noh')
    v.pvmd.__call__('mol addrep top')
    v.pvmd.__call__('mol modstyle 1 top NewRibbons 0.300000 10.000000 3.000000 0')
    
    
    angles = numpy.arange(-90,+100,20) 
    for angle in angles:
        print '========================='
        print 'ANGLE   %f' %angle
        xyztm = geometry.rotate(xyztmcopy,axis,angle,origin = center, radians=False)
        xyzal[indices,:] = xyztm
        en.setXyz(xyzal)
        v.minimize(en,200,plot=plot)

if __name__ == '__main__':
    from pizza import matlab
    m=None
    #m = matlab.matlab()
    run(sys.argv[1],plot=m)
    #m.stop()
