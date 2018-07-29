''' solvate.py contains solvation tools and tools for periodic boxes
note that boxsize here is using the NAMD way, there is an example in defaults.py'''


import numpy, copy, random, sys, math
import bgf, defaults, atoms, periodic
# from packages import kdtree # is only used in removeOverlappingWater

def putCoordinatesIntoBoxButKeepHydrogensNearBase(al,boxSize):
    """NAMD needs to have the hydrogen atom close to it's base atom in absolute numbers, even though box boundary is crossed"""
    periodic.putCoordinatesIntoBox(al,boxSize)
    
    # go through possible relative positions of the hydrogen and it's base atom
    a = boxSize[0,:]
    b = boxSize[1,:]
    c = boxSize[2,:]
    z = numpy.array([0.0,0.0,0.0])
    pos = [a,b,c,-a,-b,-c,a+b,a-b,-a+b,-a-b,a+c,a-c,-a+c,-a-c,c+b,c-b,-c+b,-c-b,+a+b+c,+a+b-c,+a-b+c,+a-b-c,-a+b+c,-a+b-c,-a-b+c,-a-b-c]
    lpos = range(len(pos))
    
    totalH = 0
    improvedH = 0
    
    for a in al:
        if a.aName.strip()[0] == 'H':
            counted = False
            totalH = totalH + 1
            b = a.bonds[0]
            best = z
            bestval = numpy.linalg.norm(a.xyz-b.xyz)
            for p in pos:
               val = numpy.linalg.norm(a.xyz+p-b.xyz)
               if val < bestval:
                    print '    improving hydrogen: %s, from %f to %f' % (str(a),bestval,val)
                    best = p
                    bestval = val
                    if not counted:
                        counted = True
                        improvedH =  improvedH + 1
            a.xyz = a.xyz + best
    
    print '  total number of hydrogen atoms  %d' % totalH
    print '  total number of hydrogens which were corrected %d' % improvedH

    
def getWaterBox(boxSize):
    al = getSmallBox()
    sA = defaults.waterBoxSize[0,:]
    sB = defaults.waterBoxSize[1,:]
    sC = defaults.waterBoxSize[2,:]
    sO = defaults.waterBoxSize[3,:]
    lA = boxSize[0,:]
    lB = boxSize[1,:]
    lC = boxSize[2,:]
    lO = boxSize[3,:]    
    
    if lC[0]!=0 or lC[1]!=0 or lB[2]!=0 or lA[2]!=0:
        print 'only can lay basic water box when certain coordinates are 0'
        sys.exit(1)
    if lA[1] != 0:
        print 'only can lay basic water box when lA[1] is zero (simpler to implement)'
        sys.exit(1)
    if lB[0] < 0:
        print 'only can lay basic water box when lB[0] < 0 is zero (simpler to implement)'
        sys.exit(1)
    if lA[0] < 0 or lB[1] < 0:
        print 'only can lay basic water box when lA[0] > 0 and lB[1] > 0 (simpler to implement)'
        sys.exit(1)
    
    xyz = atoms.getXyzFromList(al)
    xyz = xyz - sO + (sA+sB+sC)/2.0 
    atoms.putXyzIntoList(al,xyz)
    
    nO = lO-(lA+lB+lC)/2.0
    nABC = nO + lA+lB+lC
    nal = []
    for i in numpy.arange(nO[0],nABC[0],sA[0]):
        for j in numpy.arange(nO[1],nABC[1],sB[1]):
            for k in numpy.arange(nO[2],nABC[2],sC[2]):
                vec = numpy.array([i,j,k])
                nal.extend( periodic.copyTo(al, vec) )
    # trim waters and renumber  their residue numbers
    print 'number of waters in the construction box: %d' % (len(nal)/3)
    sel = []
    rNo = 0
    for a in nal:
        if a.ffType.strip() == 'OW' and a.rName.strip() == 'WAT':
            if periodic.isInsideBox(a.xyz, boxSize):
                if len(a.bonds) != 2:
                    print 'this is not water!!!!! % s ' % str(a)
                    sys.exit(1)
                h1 = a.bonds[0]
                h2 = a.bonds[1]
                rNo = rNo + 1
                a.rNo = rNo
                h1.rNo = rNo
                h2.rNo = rNo
                sel.extend([a,h1,h2])

    atoms.renumberAtomsFrom(sel)
    print 'number of waters in the produced box: %d' % (len(sel)/3)
    return sel

def getSmallBox():
    al = bgf.readFile(defaults.waterBox)
    return al

def waterPositions(al):
    waterPos = []
    for i in range(len(al)):
        a = al[i]
        if a.ffType.strip() == 'OW' and a.rName.strip() == 'WAT':
            waterPos.append(i)
    return waterPos

def placeIons(al, concentration, extraCl = 0, extraNa = 0):
    waterPos = waterPositions(al)
    numWaters = len(waterPos)
    originalLength = len(al)
    #computenumber of places nacl for the concentration
    nIon = int( concentration * numWaters / (55.5 + 2.0 * concentration) )
    #
    print 'placing ions for molality %f mol / (kg of water)' % concentration 
    print '  there is %d waters' % numWaters
    print '  for this molality one needs %d NaCl groups' % nIon
    print '  total number of Cl requested :   %d + %d = %d' % (nIon,extraCl,extraCl+nIon)
    print '  total number of Na requested :   %d + %d = %d' % (nIon,extraNa,extraNa+nIon) 
        
    placeNaCl(al,waterPos,extraCl+nIon,extraNa+nIon)
    expected = originalLength - 2* (extraCl+nIon+extraNa+nIon)
    print 'checking the final length of the atom list (expected %d):' %expected
    if len(al) == expected:
        print '    OK'
    else:
        print '    wrong: got %d  !!!' %len(al)
        sys.exit(1)

def placeNaCl(al, waterPos, ncl, nna):
    waterPos = random.sample(waterPos,ncl+nna)
    print 'removing waters with the following indices'
    print waterPos
    for i in waterPos[0:ncl]:
        ox = al[i]
        ox.aName = 'Cl-'
        ox.rName = 'Cl-'
        ox.ffType = 'IM'
        ox.charge = -1.0        
    for i in waterPos[ncl:(ncl+nna)]:
        ox = al[i]
        ox.aName = 'Na+'
        ox.rName = 'Na+'
        ox.ffType = 'IP'
        ox.charge = +1.0
    waterPos.sort(reverse=True)
    for i in waterPos:
        ox = al[i]
        al.remove(ox.bonds[0])
        al.remove(ox.bonds[1])
        ox.deleteBond(ox.bonds[0])
        ox.deleteBond(ox.bonds[0])
    
def removeOverlappingWater(al,water):
    """remove waters which have oxygen closer than 2.5 to any atom in the list al"""
    from packages import kdtree
    radius = 2.5
    kd = kdtree.KDTree(3,10) # dimension = 3, bucket size = 10
    xyz = atoms.getXyzFromList(water)
    oxygens = []
    for i in range(len(water)):
        a = water[i]
        if a.ffType.strip() == 'OW' and a.rName.strip() == 'WAT':
            oxygens.append(i)    
    kd.set_coords(xyz[oxygens,:])
    watersToRemove = []
    for a in al:
        # test if there is any water near atom a 
        kd.search(a.xyz, radius)
        # collect the results
        result_indices = kd.get_indices()
        watersToRemove.extend( list(result_indices))
    watersToRemove = set(watersToRemove)
    watersToRemove = list(watersToRemove)
    keepwater = [True] * len(oxygens)
    for i in watersToRemove:
        keepwater[i] = False
    rNo = 0
    newwater = []
    for i in range(len(oxygens)):
        if keepwater[i]:
            ox = water[oxygens[i]]
            h1 = ox.bonds[0]
            h2 = ox.bonds[1]
            rNo = rNo + 1
            ox.rNo = rNo
            h1.rNo = rNo
            h2.rNo = rNo
            newwater.extend([ox,h1,h2])
    return newwater

