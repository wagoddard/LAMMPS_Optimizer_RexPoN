import numpy
import atoms, periodic
# from packages import kdtree # only used in findBondsBetween

def findBonds(al, cutoff = 2.0, maxBonds = 4):
    ''' finds bonds within atoms in list al'''
    xyz = atoms.getXyzFromList(al)
    findBondsBetween(al, al, xyz, cutoff = cutoff, maxBonds = maxBonds, watchForSelf = True)

def findBondsBetween(al, al2, xyz, cutoff = 2.0, maxBonds = 4, watchForSelf = True):
    from packages import kdtree
    if watchForSelf:  # atom itself will be found, so we need to increase the maxBonds counter
        maxBonds = maxBonds + 1
        
    kd = kdtree.KDTree(3,10) # dimension = 3, bucket size = 10
    kd.set_coords(xyz)
    
    counter = 0
    for a in al: 
        kd.search(a.xyz, cutoff)
        result_indices = kd.get_indices()
        
        distances = numpy.array( [a.distance(al2[b]) for b in result_indices] )
        best = distances.argsort()
        
        for i in best[0:maxBonds]:  
            b = al2[result_indices[i]]
            if not a.isBonded(b) and a!=b:
                a.makeBond(b)
                counter = counter + 1
    print 'Found %d new bonds.' % counter

def findBondsAcrossBoundary(al, boxSize, cutoff = 2.0, maxBonds = 4):
    
    periodic.putCoordinatesIntoBox(al,boxSize)
    xyz = atoms.getXyzFromList(al)
    
    # cell vectors are boxSize[i,:]   for i 0,1,2
    a = boxSize[0,:]
    b = boxSize[1,:]
    c = boxSize[2,:]
    
    n = len(al)

    points = numpy.zeros((0,3))
    
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                if i!=0 or j!=0 or k!=0:
                    v = i*a+j*b+k*c
                    findBondsBetween(al, al, xyz+v, cutoff = cutoff, maxBonds = maxBonds, watchForSelf = False)

