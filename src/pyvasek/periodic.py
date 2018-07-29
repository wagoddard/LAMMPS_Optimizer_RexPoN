''' tools for periodic boxes
note that boxsize here is using the NAMD notation, there is an example in defaults.py

waterBoxSize= numpy.array([[30.0,0.0,0.0],\
                           [0.0,30.0,0.0],\
                           [0.0,0.0,30.0],\
                           [0.0,0.0,0.0]])

3 cell vectors a,b,c
1 cell origin o

a_x a_y a_z
b_x b_y b_z
c_x c_y c_z
o_x o_y o_z

in many case I use only orthogonal boxes, in which case a has to be in x direction,
b in y direction and c in z direction

'''


import numpy, copy, random, sys, math
import bgf, defaults, atoms, solvate, bonds, geometry


def copyTo(al, vec):
    nal = copy.deepcopy(al)
    xyz = atoms.getXyzFromList(nal)
    xyz = xyz + vec
    atoms.putXyzIntoList(nal,xyz)
    return nal

def getBoxCoordinates(point, boxSize):
    '''returns fractional coordinates in lattice
    vectorized, so point can be more points'''
    orig = boxSize[3,:]
    p = point - orig
    x = numpy.dot( p, numpy.linalg.inv(boxSize[0:3,:]))
    return x

def coordinatesModuloBox(point, boxSize):
    '''puts coordinates of points into the periodic box
    vectorixed, so point can be more points'''
    x = getBoxCoordinates(point, boxSize)
    intx = numpy.floor(0.5+x)
    x = x - intx
    origin = boxSize[3,:]
    xInBox = numpy.dot(x,boxSize[0:3,:])
    return origin + xInBox 

def framesModuloBox(frames, boxSize):
    for i in range(len(frames)):
        frames[i] = coordinatesModuloBox(frames[i], boxSize)

def putCoordinatesIntoBox(al,boxSize):
    for a in al:
        a.xyz = coordinatesModuloBox(a.xyz, boxSize)

def isInsideBox(point, boxSize):
    '''not vectorized yet'''
    x = getBoxCoordinates(point, boxSize)
    if -0.5<=x[0]<0.5 and -0.5<=x[1]<0.5 and  -0.5<=x[2]<0.5:
        return True
    else:
        return False   

def replicate(al,boxSize,na,nb,nc,cutoff = 2.0, maxBonds = 4,includePeriodic=True):
    '''na=nb=nc=1 would just keep the same atoms'''
    putCoordinatesIntoBox(al,boxSize)
    xyz = atoms.getXyzFromList(al).copy()
    
    cal = copy.deepcopy(al)
    
    # cell vectors are boxSize[i,:]   for i 0,1,2
    a = boxSize[0,:]
    b = boxSize[1,:]
    c = boxSize[2,:]
    
    for i in range(1,na):
        v = i*a
        cal2 = copy.deepcopy(cal)
        geometry.moveAtoms(cal2, v)
        al.extend(cal2)
    cal = copy.deepcopy(al)
    for i in range(1,nb):
        v = i*b
        cal2 = copy.deepcopy(cal)
        geometry.moveAtoms(cal2, v)
        al.extend(cal2)
    cal = copy.deepcopy(al)
    for i in range(1,nc):
        v = i*c
        cal2 = copy.deepcopy(cal)
        geometry.moveAtoms(cal2, v)
        al.extend(cal2)
    
    newbox = numpy.vstack([na*a, nb*b, nc*c, boxSize[3,:]])
    putCoordinatesIntoBox(al,newbox)
    bonds.findBonds(al, cutoff=cutoff, maxBonds=maxBonds)
    
    if includePeriodic:
        bonds.findBondsAcrossBoundary(al, newbox, cutoff=cutoff, maxBonds=maxBonds)
        
    return newbox
    
def graphene(na,nb,nc,includePeriodic = True):
    lattice = 1.42
    planar = 3.35
    a = numpy.array([lattice*math.sqrt(3),0.,0.])
    b = numpy.array([0.,3.*lattice,0.])
    c = numpy.array([0.,0.,planar])
    boxSize = numpy.vstack([a,b,c,numpy.array([0.,0.,0.])])
    # make two basic atoms
    a = atoms.Atom()
    a.aName = 'C'
    a.rName = 'GRA'
    a.rNo = 1
    a.chain = 'G'
    a.ffType = 'ca'
    b = copy.deepcopy(a)
    b.xyz = numpy.array([0.,lattice,0.])
    a.makeBond(b)
    al = [a,b]
    cal = copy.deepcopy(al)
    geometry.moveAtoms(cal,numpy.array([lattice*math.sqrt(3)/2.0,lattice * 1.5 ,0.]))
    al.extend(cal)
    
    boxSize = replicate(al,boxSize,na,nb,nc,cutoff=lattice*1.05, maxBonds=3,includePeriodic=includePeriodic)    
    return (al,boxSize)

def grapheneDelete(al,carbonsToDelete):
    
    for c in carbonsToDelete:
        c.ffType = 'ha'
        c.aName = 'H'
    for c in carbonsToDelete:
        hasOutsideNeighbor = False
        for b in c.bonds[:]:
            if b.ffType == 'ha':
                b.deleteBond(c)
            else:
                hasOutsideNeighbor = True
        if not hasOutsideNeighbor:
            al.remove(c)
    # set the charges
    #for a in al:
    #    if a.ffType == 'ha':
    #        a.charge = 0.12
    #        for b in a.bonds:
    #            b.charge = -0.12 

def grapheneHole(al,center,radius):
    todelete = []
    for a in al:
        if numpy.linalg.norm(a.xyz-center) < radius:
            todelete.append(a)
    grapheneDelete(al,todelete)
    
def grapheneFixHydrogen(al):
    '''fixes graphene hydrogen distance and charge'''
    for a in al:
        if a.ffType == 'ha':
            a.charge = 0.12
            if len(a.bonds) > 1:
                print 'the movable atom and the fixed atom are connected by more then a direct bond, so cannot move the movable atom'
                sys.exit(1)
            else:
                b = a.bonds[0]
                b.charge = -0.12
                geometry.setBondLength(b,a,1.1)
                        
def readXsc(xscfile):
    '''read namd xsc file to get coordinates '''
    lines = open(xscfile).readlines()
    for line in lines:
        if line[0:1] != '#':
            box = numpy.array([ float(a) for a in line.split() ])
    box = numpy.reshape(box[1:13],(4,3))
    return box

def convertCellInto(a,b,c,alphaDeg,betaDeg,gammaDeg):
    '''converts cell size in the crystalographic coordinates to the NAMD box size I use 
    i.e. each cell vector is specified with 3 coordinates
    alpha, beta, gamma are in degrees'''
    # convert angles to radians
    alpha = alphaDeg/180.0*math.pi
    beta  = betaDeg/180.0*math.pi
    gamma = gammaDeg/180.0*math.pi
    # convert crystal cell a,b,c,alpha,beta,gamma to the namd cell notation
    ca = math.cos(alpha)
    cb = math.cos(beta)
    cg = math.cos(gamma)
    sg = math.sin(gamma)
    box = numpy.zeros((4,3))
    box[0,0] = a
    box[1,0] = b*cg
    box[1,1] = b*sg
    v = math.sqrt(1-math.pow(ca,2)-math.pow(cb,2)-math.pow(cg,2)+2*ca*cb*cg)
    box[2,0] = c*cb
    box[2,1] = c*(ca-cb*cg)/sg
    box[2,2] = c*v/sg
    # fix the round off error
    largenumber = abs(box[0,0])+abs(box[1,1])+abs(box[2,2])
    smallnumber = abs(box[1,0])+abs(box[2,0])+abs(box[2,1])
    if smallnumber/largenumber < 1e-10:
        box[1,0] = 0.0
        box[2,0] = 0.0
        box[2,1] = 0.0
    # origin in the a,b,c axis system is 0,0,0
    # but origin in the namd notation is in the middle of the cell
    box[3,:] = (box[0,:]+box[1,:]+box[2,:])/2
    return box

def readBgf(bgfLines):
    '''read the cell size from list of lines from a bgf file
    only reads CRYSTX entry
    returns None if there is no CRYSTX entry
    
    'CRYSTX', Specification of unit cell parameters, lengths of axes (a,b,c) in A and angles 
    (alpha, beta, gamma) in degrees.    
    '''
    # find crystx line:
    crystx = []
    for line in bgfLines:
        if line[0:6]=='CRYSTX':
            crystx = line.split()
            break
    if len(crystx) == 0:
        return None
    # get crystal cell parameters
    a = float(crystx[1])
    b = float(crystx[2])
    c = float(crystx[3])
    alpha = float(crystx[4])
    beta  = float(crystx[5])
    gamma = float(crystx[6])
    return convertCellInto(a,b,c,alpha,beta,gamma)

def unwrap(xyzs,boxSizes):
    '''xyzs are multiple frames of some trajectory
    boxSizes are box sizes from those steps
    
    returned xyzs2 contains unwrapped coordinates 
    (so that protein is in one piece and not broken on periodic cell bounrady)
    works only for rectangular cells
    '''
    xyzs2 = []
    xyzs2.append(xyzs[0])
    for i in range(len(xyzs)-1):
        prevbox = boxSizes[i]
        prevxyz = xyzs2[-1]
        prevx = getBoxCoordinates(prevxyz, prevbox)

        nextbox = boxSizes[i+1]
        nextxyz = xyzs[i+1]
        nextx = getBoxCoordinates(nextxyz, nextbox)

        x = prevx-nextx
        intx = numpy.floor(0.5+x)

        #if numpy.linalg.norm( intx ) > 0:
        #    print intx

        x = nextx + intx
        xInBox = numpy.dot(x, nextbox[0:3,:])
        origin = nextbox[3,:]

        xyzs2.append(origin + xInBox)
    return numpy.array(xyzs2)




