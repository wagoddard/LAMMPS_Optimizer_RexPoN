"""
in this module there are functions for
- rotation of atoms in 3d
- alignment of one atom selection to another
- filling in missing atoms positions for various hybridizations (sp2, sp3)
"""
import math, numpy, copy, sys
import atoms, residues
# from packages import kdtree # only used in findStericClash

def rotate(xyz,axis,angle,origin = None, radians=False):
    '''returns new array
    computes vector which is produced by rotating vector A around vector B (axis), using given origin by given:
    
    C = cos(angle) * A + (1-cos(angle)) (A dot B) * B / abs(B)^2 + sin(angle) * (A cross B) / abs(B)
    '''
    if not radians:
        angle = angle/180.0*math.pi
    if origin != None:
        nxyz = xyz-origin
    else:
        nxyz = xyz.copy()
    
    naxis = axis / numpy.linalg.norm(axis)
    cos = math.cos(angle)
    sin = -math.sin(angle)  # fixing the orientation of the rotation
    
    dot = numpy.dot(nxyz,naxis)
    dots = numpy.array([dot,dot,dot]).T
    
    nxyz = cos*nxyz + (1-cos)*dots*naxis + sin*numpy.cross(nxyz,naxis)
    if origin != None:
        nxyz = nxyz+origin

    return nxyz
    
def rotateAtoms(al, axis,angle,origin = None, radians=False):
    xyz = atoms.getXyzFromList(al)
    xyz = rotate(xyz,axis,angle,origin = origin, radians=radians)
    atoms.putXyzIntoList(al, xyz)

def moveAtoms(al, vector):
    xyz = atoms.getXyzFromList(al)
    xyz = xyz + vector
    atoms.putXyzIntoList(al, xyz)
    
def setBondLength(atomFixed,atomMovable,distance):
    # get atoms which are only attached to the movable atom
    movable = [atomMovable]
    fifo = []
    for b in atomMovable.bonds:
        if b != atomFixed:
            fifo.append(b)
    while len(fifo) > 0:
        a = fifo.pop()
        if a==atomFixed:
            print 'the movable atom and the fixed atom are connected by more then a direct bond, so cannot move the movable atom'
            sys.exit(1)
        if not a in movable:
            movable.append(a)
            fifo.extend(a.bonbds)
    v0 = atomFixed.xyz
    v1 = atomMovable.xyz
    v2 = v1-v0
    v2 = v2*(distance/numpy.linalg.norm(v2) - 1)
    moveAtoms(movable,v2)
    
def findStericClash(al,distance = 3.0):
    '''find pairs of atoms which are not connected by a bond or don't have a common neighbor which are closer than the given distance '''
    from packages import kdtree
    xyz = atoms.getXyzFromList(al)
    
    kd = kdtree.KDTree(3,10) # dimension = 3, bucket size = 10
    kd.set_coords(xyz)
    
    counter = 0
    for i in range(len(al)):
        a = al[i] 
        kd.search(a.xyz, distance)
        result_indices = kd.get_indices()
        for j in result_indices:
            if i<j:
                if not a.lessThanThreeBondsFrom(al[j]):
                    print 'atoms %s (index %d) and %s (index %d) are at distance %f ' % (a, i, al[j],j,a.distance(al[j]))        
                    counter += 1
    print '  summary: in total found %d pairs of atoms (nonbonded or separated by more than two bonds) closer than %f' % (counter, distance)


def align( set1_orig, set2_orig, selection = None, selection1 = None, selection2 = None, printWarning = True):
    """
    align performs the Kabsch alignment algorithm
    aligns set2 to set1

    set1 and set2 are xyz numpy arrays with atom coordinates

    returns: coordinates of aligned set2, and 3 RMSD numbers:
                    - between the aligned points (in selection) 
                    - between the other points
                    - between all points

    optional: selection is numpy array with indices of atoms to align
                (in this case RMSD on otehr points than selection is returned)
    optional2: apply different selections to the two sets
                (in this case RMSD on otehr points than selection is NOT returned)

    adopted from Jason Vertrees: QKabsch.py. 
    """

    # get local copies to manipulate the data
    if selection != None:
        set1 = set1_orig[selection]
        set2 = set2_orig[selection]
    elif selection1 != None and selection2 != None:
        set1 = set1_orig[selection1]
        set2 = set2_orig[selection2]
    elif selection1 != None:
        set1 = set1_orig[selection1]
        set2 = set2_orig.copy()
    elif selection2 != None:
        set1 = set1_orig.copy()
        set2 = set2_orig[selection2]
    else:
        set1 = set1_orig.copy()
        set2 = set2_orig.copy()

    # check for consistency
    assert len(set1) == len(set2)
    L = len(set1)
    assert L > 0

    # must alway center the two proteins to avoid
    # affine transformations.  Center the two proteins
    # to their selections.

    COM1 = numpy.sum(set1,axis=0) / float(L)
    COM2 = numpy.sum(set2,axis=0) / float(L)
    set1 -= COM1
    set2 -= COM2

    # Initial residual, see Kabsch.
    E0 = numpy.sum( numpy.sum(set1 * set1,axis=0),axis=0) + numpy.sum( numpy.sum(set2 * set2,axis=0),axis=0)

    #
    # This beautiful step provides the answer.  V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(set2), set1))

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]
        if printWarning: print 'Warning in residues.align(): The best match is reflection. Need to correct that.'

    RMSD = E0 - (2.0 * sum(S))
    RMSD = numpy.sqrt(abs(RMSD / L))

    #U is simply V*Wt
    U = numpy.dot(V, Wt)
    # this might come handy: angle of rotation of the 3d rotation
    # rotationMatrix = U
    # trace = numpy.trace(U)  # trace = 1+2*cos(angle)
    # angleDeg = numpy.arccos( (trace-1.)/2. )*180./numpy.pi


    # rotate and translate the molecule
    set2 = numpy.dot(set2_orig - COM2, U) + COM1
    # center the molecule

    # compue other RMSDs if the alignment was done only on selection
    if selection != None:
        n = len(set1_orig)
        nsel = len(selection)

        RMSDallsn = numpy.sum(numpy.power(set2-set1_orig,2))
        RMSDall = numpy.sqrt(RMSDallsn / n)

        if n-nsel > 0:
            tempnumber = ( RMSDallsn-nsel*numpy.power(RMSD,2) )/(n-nsel)
            if tempnumber > 1e-10:
                RMSDnonsel = numpy.sqrt( ( RMSDallsn-nsel*numpy.power(RMSD,2) )/(n-nsel) )    
            else:
                RMSDnonsel = 0.0
        else:
            RMSDnonsel = 0.0
    else:
        RMSDnonsel = 0.0
        RMSDall = RMSD 

    # print "RMSD of the Kabsch allignment =%f" % RMSD
    return set2, RMSD, RMSDnonsel, RMSDall

def alignBackboneAtoms(al1,al2):
    """aligns al2 to the backbone of al1
    the two provided atoms list do not need to have the same number of atoms
    but they have to have the same number of residues 
    returns al2 list"""
    backbone = ['N','CA','C','O']
    xyz2 = atoms.getXyzFromList(al2)
    xyz1 = xyz2.copy()
    rl1 = residues.makeList(al1)
    rl2 = residues.makeList(al2)
    selection = []
    if len(rl1) != len(rl2):
        sys.exit('Error in alignBackboneAtoms: the two list have different numebr of residues: %d and %d' % (len(rl1), len(rl2)))
    for i in range(len(rl1)):
        if rl1[i].rName == 'RES':
            continue
        r1 = rl1[i]
        r2 = rl2[i]
        for aName in backbone:
            a1 = r1.atomWithName(aName)
            a2 = r2.atomWithName(aName)
            i2 = al2.index(a2)    
            selection.append(i2)
            xyz1[i2] = a1.xyz
    xyz3, RMSD, RMSDnonsel, RMSDall = align( xyz1, xyz2, selection = selection)
    atoms.putXyzIntoList(al2,xyz3)
    return al2 

def alignAtoms(al1,al2,sel1=None,sel2=None):
    """ another interface to the align function
    will align atoms in al2 to atoms in al1
    optional: sel1: list of indices af atoms in al1 (does not have to be numpy.array)
    optional: sel2: list of indices af atoms in al2 (does not have to be numpy.array)
    lengths of these two selections have to be the same
    """
    if sel1 == None:
        sel1 = range(len(al1))
    if sel2 == None:
        sel2 = range(len(al2))
    assert len(sel1) == len(sel2), 'Error: lenghts of the two selections in geometry.alignAtomsAngle is not the same: %d and %d' % (len(sel1), len(sel2))

    # create a list of dumme atoms with the target orientations that we are aligning to
    bl2 = []
    for i in range(len(al2)):
        b = atoms.Atom()
        bl2.append(b)

    for i in range(len(sel1)):
        b = bl2[sel2[i]]
        b.xyz = al1[sel1[i]].xyz.copy()

    set1_orig = atoms.getXyzFromList(bl2)
    set2_orig = atoms.getXyzFromList(al2)    
    
    set2, RMSD, RMSDnonsel, RMSDall = align( set1_orig, set2_orig, selection = numpy.array(sel2))    
    atoms.putXyzIntoList(al2,set2)

    return RMSD 

# atom building functions
def missingAtomXyz(atom, sp, bondLength):
    if sp == 'R':
        sp = '2'
    
    if sp == '1':
        if len(atom.bonds) == 1:
            atom2 = atom.bonds[0]
            lpairlist = [atom.xyz - bondLength * (atom2.xyz - atom.xyz) / numpy.linalg.norm(atom2.xyz - atom.xyz)]
        else:
            print >> sys.stderr, 'sp1 atom has to have 1 bond so that I can find new atom positions, but has %d' % len(atom.bonds)
            sys.exit(1)
    elif sp == '2':
        if len(atom.bonds) == 1:
            atom2 = atom.bonds[0]
            if len(atom2.bonds) == 1:
                 atom3xyz = numpy.array([10000.,0.,0.])
            else:
                for atom3 in atom2.bonds:
                    if atom3 != atom:
                        break
                atom3xyz = atom3.xyz
            lpairlist = sp2given1(atom.xyz,atom2.xyz,atom3xyz,bondLength)
        elif len(atom.bonds) == 2:
            atom2 = atom.bonds[0]
            atom3 = atom.bonds[1]
            lpairlist = [sp2given2(atom.xyz,atom2.xyz,atom3.xyz,bondLength)]
        else:
            print >> sys.stderr, 'sp2 atom has to have 1 or 2 bonds so that I can find new atom positions, but has %d' % len(atom.bonds)
            sys.exit(1)
    elif sp == '3':
        if len(atom.bonds) == 1:
            atom2 = atom.bonds[0]
            if len(atom2.bonds) == 1:
                 atom3xyz = numpy.array([10000.,0.,0.])
            else:
                for atom3 in atom2.bonds:
                    if atom3 != atom:
                        break
                atom3xyz = atom3.xyz
            lpairlist = sp3given1(atom.xyz,atom2.xyz,atom3xyz,bondLength)
        elif len(atom.bonds) == 2:
            atom2 = atom.bonds[0]
            atom3 = atom.bonds[1]
            lpairlist = sp3given2(atom.xyz,atom2.xyz,atom3.xyz,bondLength)
        elif len(atom.bonds) == 3:
            atom2 = atom.bonds[0]
            atom3 = atom.bonds[1]
            atom4 = atom.bonds[2]
            lpairlist = [sp3given3(atom.xyz,atom2.xyz,atom3.xyz,atom4.xyz,bondLength)]
        else:
            print >> sys.stderr, 'sp3 atom has to have 1,2 or 3 bonds so that I can find new atom positions, but has %d' % len(atom.bonds)
            sys.exit(1)
    else:
        print >> sys.stderr, 'sp has to be 1,2,3 or R, but it is %s' % (sp)
        sys.exit(1)
    return lpairlist 

def sp2given2(centerxyz,atom1xyz,atom2xyz,bondLength):
    '''returns one xyz'''
    a1 = atom1xyz-centerxyz
    a1 /= numpy.linalg.norm(a1)
    a2 = atom2xyz-centerxyz
    a2 /= numpy.linalg.norm(a2)
    a3 = -a1-a2
    a3 /= numpy.linalg.norm(a3)
    return centerxyz + a3 * bondLength

def sp2given1(centerxyz,bondedAtomxyz,atomInPlanexyz,bondLength):
    '''returns tuple of length 2
    first returned coordinate has dihedral 0 with atomInPlanexyz
    second returned coordinate has dihedral 180 with atomInPlanexyz
    '''
    a1 = bondedAtomxyz - centerxyz
    a1 = a1 / numpy.linalg.norm(a1) * bondLength
    a2 = atomInPlanexyz - bondedAtomxyz
    a3 = numpy.cross(a1,a2)
    a4 = rotate(a1,a3,120, radians=False)
    a5 = rotate(a1,a3,-120, radians=False)
    return a4+centerxyz,a5+centerxyz

def sp3given2(centerxyz,atom1xyz,atom2xyz,bondLength):
    '''returns tuple of length 2
    if called with (ca.xyz,n.xyz,c.xyz,bondLength)
    then the first returned direction is CB of L amino acid, second direction is direction of HA    
    '''
    a1 = atom1xyz-centerxyz
    a1 /= numpy.linalg.norm(a1)
    a2 = atom2xyz-centerxyz
    a2 /= numpy.linalg.norm(a2)
    a3 = -a1-a2
    a3 /= numpy.linalg.norm(a3)
    a4 = numpy.cross(a1,a2)
    a4 /= numpy.linalg.norm(a4)
    angle = 109.4712 / 2.0 / 180.0 * math.pi
    a5 = centerxyz + bondLength *(math.sin(angle) * a4 + math.cos(angle) *a3)
    a6 = centerxyz + bondLength *(-math.sin(angle) * a4 + math.cos(angle) *a3)
    return  a5,a6

def sp3given3(centerxyz,atom1xyz,atom2xyz,atom3xyz,bondLength):
    '''returns one xyz'''
    a1 = atom1xyz-centerxyz
    a1 /= numpy.linalg.norm(a1)
    a2 = atom2xyz-centerxyz
    a2 /= numpy.linalg.norm(a2)
    a3 = atom3xyz-centerxyz
    a3 /= numpy.linalg.norm(a3)
    a4 = a1+a2+a3
    a4 /= numpy.linalg.norm(a4)    
    return centerxyz - a4 * bondLength

def sp3given1(centerxyz,bondedAtomxyz,staggeredAtomxyz,bondLength):
    '''returns tuple of length 3'''
    axis = numpy.cross(centerxyz-bondedAtomxyz,staggeredAtomxyz-bondedAtomxyz)
    a1 = rotate(bondedAtomxyz,axis,109.4712,origin = centerxyz, radians=False)
    a2 = a1 - centerxyz
    a2 /= numpy.linalg.norm(a2)
    a3 = centerxyz + bondLength * a2
    a4 = rotate(a3,centerxyz-bondedAtomxyz,120,origin = centerxyz, radians=False)
    a5 = rotate(a3,centerxyz-bondedAtomxyz,-120,origin = centerxyz, radians=False)
    return a3,a4,a5

def angleXyz(atom1xyz,atom2xyz,atom3xyz):
    '''returns angle atom1-atom2-atom3 (atom2 is the central atom)'''
    v1 = atom1xyz-atom2xyz
    v2 = atom3xyz-atom2xyz
    return math.acos(numpy.dot(v1,v2) / numpy.linalg.norm(v1) / numpy.linalg.norm(v2) )

def angle(atom1,atom2,atom3,degrees=False):
    '''returns angle atom1-atom2-atom3 (atom2 is the central atom)'''
    if degrees:
        return angleXyz(atom1.xyz,atom2.xyz,atom3.xyz)*180.0/math.pi
    else:
        return angleXyz(atom1.xyz,atom2.xyz,atom3.xyz)

def torsionXyz(atom1xyz,atom2xyz,atom3xyz,atom4xyz):
    b3 = atom4xyz-atom3xyz
    b2 = atom3xyz-atom2xyz
    b1 = atom2xyz-atom1xyz
    b2length = numpy.linalg.norm(b2)
    return math.atan2( b2length*numpy.dot(b1,numpy.cross(b2,b3)), numpy.dot(numpy.cross(b1,b2),numpy.cross(b2,b3)) )

def findAngleOfRotation(axis,point1,point2):
    """finds angle of rotation to get point 2 from point 1 along the axis
    center of rotation is 0,0,0
    returned angle is in radians
    """  
    normaxis = numpy.linalg.norm(axis)
    if normaxis < 1.0e-10:
        print 'ERROR: the axis has zero length so I cannot find rotation angle'
        sys.exit(1) 
    naxis = axis/normaxis
    perp1 = point1 - axis * numpy.dot(axis, point1) 
    perp2 = point2 - axis * numpy.dot(axis, point2)
    cosarg = numpy.dot(perp1,perp2) / numpy.linalg.norm(perp1) / numpy.linalg.norm(perp2) 
    angle = math.acos( cosarg)
    if numpy.dot(axis,numpy.cross(perp1,perp2)) > 0:
        return angle
    else:
        return -angle

def findRotation(a1,b1,c1,a2,b2,c2):
    """for two rigid bodies 1 and 2 defined by 3 points a1,b1,c1 and a2,b2,c2 finds translations and rotation that maps body 1 to body 2
    body2 = rotation(body1+translation)
    translation is a2-a1
    output is translation, axis, angle in radians
    points are inputted as xyz (1x3) numpy array"""
    translation = a2-a1    
    d1 = b1-a1
    d2 = b2-a2
    e1 = c1-a1
    e2 = c2-a2
    f1 = numpy.cross(d1,e1)
    f2 = numpy.cross(d2,e2)
    nd1 = numpy.linalg.norm(d1)
    nd2 = numpy.linalg.norm(d2)
    ne1 = numpy.linalg.norm(e1)
    ne2 = numpy.linalg.norm(e2)
    nf1 = numpy.linalg.norm(f1)
    nf2 = numpy.linalg.norm(f2)
    assert abs(nd2-nd1)+abs(ne2-ne1)+abs(nf2-nf1) < 1.0e-10, 'ERROR: the inputted vectors do no represent srigid body rotation %f, %f, %f' % (abs(nd2-nd1),abs(ne2-ne1),abs(nf2-nf1))
    d1 /= nd1
    d2 /= nd2
    e1 /= ne1
    e2 /= ne2
    f1 /= nf1
    f2 /= nf2
    g1 = numpy.vstack([d1,e1,f1]).T    
    g2 = numpy.vstack([d2,e2,f2]).T    
    # want to solve for M:   g2 = M g1
    # M = g2*inv(g1)
    M = numpy.dot(g2,numpy.linalg.inv(g1))
    # we have a matrix of the rotation
    eigval,eigvec = numpy.linalg.eig(M) 
    eigvalreal = numpy.real(eigval)
    index = 0
    for i in range(len(eigval)):
        if abs(eigval[i]-1.0) < 1.0e-10:
            index = i
    if abs(eigval[index]-1.0) > 1.0e-10:
        print 'ERROR: i did not find eigenvalue 1 for this matrix, i.e. it is not rotation '
        sys.exit(1)
    axis = numpy.real(eigvec[:,index])
    axis /= numpy.linalg.norm(axis)    
    x = numpy.array([1.,0.,0.])
    y = numpy.array([0.,1.,0.])
    x = x - axis*numpy.dot(axis,x)
    y = y - axis*numpy.dot(axis,y)
    if numpy.linalg.norm(y) > numpy.linalg.norm(x):
        x = y 
    xrot = numpy.dot(M,x)
    angle = findAngleOfRotation(axis,x,xrot)
    return translation,axis,angle
   
def randomRotate(xyz):
    '''randomly rotates the matrix xzy around origin
    generating a random rotation in 3d is quite tricky'''
    # first get a random point on sphere A
    #   the loop would be unnecessary if i was certain how to use sin and cos 
    #   to preserve the uniform distribution in 3d 
    r = 0.0
    while not (0.1<r<1.):
        a = 2*numpy.random.random(3) - 1.
        r = numpy.linalg.norm(a)
    a = a / r
    # then get point B, such that angle AOB is 90deg
    #   for this need a vector which is at 90 deg to a
    x1 = numpy.array([1.,0.,0.])
    x2 = numpy.array([0.,1.,0.])
    nx1 = numpy.linalg.norm(numpy.cross(a,x1))
    nx2 = numpy.linalg.norm(numpy.cross(a,x2))
    if nx1 > nx2:
        b = numpy.cross(a,x1)
    else:
        b = numpy.cross(a,x2)
    b = b / numpy.linalg.norm(b)
    #   for random point on a spehre can use uniform distribution
    randomAngle = numpy.random.rand() * math.pi * 2.
    b = rotate(b,a,randomAngle,origin=None,radians=True)
    # then find rotation mappoing Ao=[1,0,0] to A and Bo=[0,1,0] to B
    origin = numpy.array([0.,0.,0.])
    Ao =  numpy.array([1.,0.,0.])
    Bo = numpy.array([0.,1.,0.])
    translation,axis,angle = findRotation(origin,Ao,Bo,origin,a,b)    

    return rotate(xyz, axis, angle, origin=origin, radians=True)

def randomVector(dims):
    '''returns normalized random direction'''
    vec = numpy.random.normal(size=(dims,))
    vec /= numpy.linalg.norm(vec)
    return vec

def rotationMatrixFromAxis(axis,angle,radians=True):
    vec = axis/numpy.linalg.norm(axis)
    vec = numpy.reshape(vec,(3,1))
    eye = numpy.eye(3)
    uu = numpy.dot(vec,vec.T)
    ux = numpy.zeros((3,3))
    ux[0,1] = -vec[2]
    ux[1,0] = vec[2]
    ux[0,2] = vec[1]
    ux[2,0] = -vec[1]
    ux[1,2] = -vec[0]
    ux[2,1] = vec[0]
    if not radians:
        angle = angle / 180. * math.pi
    return math.cos(angle)*eye + math.sin(angle)*ux + (1-math.cos(angle))*uu

def rotationAngleFromMatrix(rotmatrix,radians=True):
    '''this might come handy: angle of rotation of the 3d rotation
    rotationMatrix = U
    trace = numpy.trace(U)  # trace = 1+2*cos(angle)
    angleDeg = numpy.arccos( (trace-1.)/2. )*180./numpy.pi'''
    
    # we have a matrix of the rotation
    eigval,eigvec = numpy.linalg.eig(rotmatrix) 
    eigvalreal = numpy.real(eigval)
    index = 0
    for i in range(len(eigval)):
        if abs(eigval[i]-1.0) < 1.0e-10:
            index = i
    assert abs(eigval[index]-1.0) < 1.0e-10, 'ERROR: i did not find eigenvalue 1 for this matrix, i.e. it is not rotation '
    axis = numpy.real(eigvec[:,index])
    axis /= numpy.linalg.norm(axis)    

    trace = numpy.trace(rotmatrix)
    angle = numpy.arccos((trace-1.)/2.)

    # possibly switch the sign of the axis
    testvector = numpy.array([1.,0.,0.])
    if abs(numpy.dot(axis,testvector)) > 0.9:
        testvector = numpy.array([0.,1.,0.])      
    assert abs(numpy.dot(axis,testvector)) < 0.9, 'error: cannot find a good test vector'
    testvector1 = numpy.dot(rotmatrix, testvector)    
    testvector2 = rotate(testvector,axis,angle,radians=True)

    if numpy.linalg.norm(testvector1-testvector2) < 1e-10:
        pass
    else:
        axis *= -1.0
        testvector2 = rotate(testvector,axis,angle,radians=True)
        assert numpy.linalg.norm(testvector1-testvector2) < 1e-10, 'error: inverting axis does not help (mismatch is %e)'%numpy.linalg.norm(testvector1-testvector2)

    if not radians:
        angle = angle * 180. / math.pi
    return axis,angle
    


