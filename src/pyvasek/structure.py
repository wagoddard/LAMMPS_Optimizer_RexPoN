import time, os, sys, time
import atoms, forceField

class Structure:
    """generates complete topology with paramater for atom list AL and force field FF"""
    def __init__(self,AL,FF,debug=False, exceptResidues=[]):
        if debug: print 'generating structure...'
        # first fix possible problems with teh atoms in the atom list        
        # 1. fix ffType so that they have proper number of characters for each force field

        forceField.fixAtomLabels(AL,FF)
        # 2. renumber atoms 
        atoms.renumberAtomsFrom(AL)
        # now generate the structure
        self.bonds = populateBonds(AL,FF)
        self.angles = populateAngles(AL,FF) # removed the exceptResidues parameter
        self.torsions = populateTorsions(AL,self.bonds,FF,exceptResidues) # each torsion is present, so there can be more torsions with the same 4 atoms 
        self.inversions = populateInversions(AL,FF)                    # each ivnersion is present 3 times for the symmetry
        # also get ready for the output: we want to know the list of used atom types
        self.usedAtomLabels = getUsedAtomLabels(AL)                   # dictionary mapping ffLabel to some index (from 0)
        self.usedBondTypes = getUsedTypes(self.bonds)                 # dictionary mapping atom BondType to some index (from 0)
        self.usedAngleTypes = getUsedTypes(self.angles)
        self.usedTorsionTypes = getUsedTypes(self.torsions)           
        self.usedInversionTypes = getUsedTypes(self.inversions)

        '''self.makeUsedLists(AL,FF)'''
        # also sort the dictionaries with the used types
        self.sortedAtomLabels = dictionaryToList(self.usedAtomLabels) # list containing the used atom types, index can be read from usedAtomTypes dictionary
        self.sortedBondTypes = dictionaryToList(self.usedBondTypes)
        self.sortedAngleTypes = dictionaryToList(self.usedAngleTypes)
        self.sortedTorsionTypes = dictionaryToList(self.usedTorsionTypes)
        self.sortedInversionTypes = dictionaryToList(self.usedInversionTypes)

        self.sortedHydrogenBondTypes = getUsedHydrogenBondTypes(FF.hydrogenBondTypes, self.usedAtomLabels) # list of hydrogen bond types

        # we are done
        if debug: 
            self.printInfo(AL)            
            print 'Structure generated.'

    def printInfo(self,AL):
            print 'Structure info:'
            print '    atoms      %6d      atom types      %6d' % (len(AL),len(self.sortedAtomLabels))
            print '    bonds      %6d      bond types      %6d' % (len(self.bonds),len(self.sortedBondTypes))
            print '    angles     %6d      angle types     %6d' % (len(self.angles),len(self.sortedAngleTypes))
            print '    torsions   %6d      torsion types   %6d' % (len(self.torsions),len(self.sortedTorsionTypes))
            print '    inversions %6d      inversion types %6d' % (len(self.inversions),len(self.sortedInversionTypes))
            print '                           hBond types     %6d' % (len(self.sortedHydrogenBondTypes))
        
class Bond:
    def __init__(self,i,j,FF):
        self.i = i
        self.j = j
        self.type = FF.getBondType(i.ffType,j.ffType)

class Angle:
    def __init__(self,i,j,k,FF):
        self.i = i
        self.j = j 
        self.k = k
        self.type = FF.getAngleType(i.ffType,j.ffType,k.ffType,skipMissing=FF.skipMissing)

class Torsion:
    def __init__(self,i,j,k,l,torsionType, avoid14 = False):
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.type = torsionType # this is a single torsion Type
        self.avoid14 = avoid14  # special value for amber, to cover the case of multiple torsions
    def hasSameAtoms(self, t2):
        if self.i == t2.i and self.j == t2.j and self.k == t2.k and self.l == t2.l:
            return True
        elif self.l == t2.i and self.k == t2.j and self.j == t2.k and self.i == t2.l:
            return True
        else:
            return False

class Inversion:
    def __init__(self,i,j,k,l,inversionType):
        self.i = i  
                    # amber official rule for atom ordering in inversion is not 
                    # kept so i include all three inversions even for amber 
                    #
                    # official amber rule:
                    #
                    #                                 J
                    #                                 |
                    #                                 K
                    #                                / \
                    #                               I L
                    #
                    #                          Improper I-J-K-L 
                    #
                    #where the central atom (K) is the third atom in the
                    #improper
                    #              and the order of the other three is determined
                    #alphabetically
                    #              by atom type and if types are the same by atom number.                  
                    # note: this order seems to be obeyed in most of the amebr library files (maybe not ARG)
                    #   but this order is not obeyed in new residues from gaff 
        self.j = j 
        self.k = k
        self.l = l
        self.type = inversionType

def populateBonds(AL,FF):
    bonds = []
    for a in AL:
        for b in a.bonds:
            if a.aNo < b.aNo: # lets count each bond only once
                bonds.append(Bond(a,b,FF))
    return bonds 

def populateAngles(AL,FF,exceptResidues=[]):
    angles = []
    for a in AL:
        if not a.rName.strip() in exceptResidues:
            for ib, b in enumerate(a.bonds):
                for c in a.bonds[ib+1:]:
                    thisangle = Angle(b,a,c,FF)
                    if thisangle.type != None:
                        angles.append(thisangle)
    return angles

def getCycles(al):
    '''helping function for identyfying exocyclic torsions in dreiding, only consider resonant atoms
    returns list of atoms which are in some cycle of length < 6'''
    haveCycle = []  
    reso = 'R'
    for a in al:
        if a.ffType[2] == reso:
            
            queue = []
            for b in a.bonds: 
                if b.ffType[2]==reso:
                    queue.append((a,b,1))
            while len(queue) > 0:
                (b,c,depth) = queue.pop()
                if c == a:
                    haveCycle.append(a)
                    break
                if depth <= 5:
                    for d in c.bonds:
                        if d!=b and c.ffType[2]==reso:
                            queue.append((c,d,depth+1))
    return haveCycle            
    

def isExocyclic(i,j,k,l,haveCycle):
    '''for dreiding only: want to find out if the torsion is exocyclic, but only care about resonant atoms'''
    reso = 'R'    
    # we are only checking if a torsion is exocyclic for resonant bonds
    if not (j.ffType[2] == reso and k.ffType[2] == reso):
        return False 
    # cases left: both j and k are resonant
    # first do simple heuristic base on number of neighbors
    # if both k and j have non sp2 neighbors then the bond is not exocyclic
    
    jNonRBond = False
    for a in j.bonds:
        if a.ffType[2] != reso:
            jNonRBond = True
            break
    kNonRBond = False
    for a in k.bonds:
        if a.ffType[2] != reso:
            kNonRBond = True
            break
    if jNonRBond and kNonRBond:
        return False

    jHasCycle = (j in haveCycle)
    kHasCycle = (k in haveCycle)
    
    if (not jHasCycle) and (not kHasCycle):
        return False
    if jHasCycle and (not kHasCycle):
        return True
    if kHasCycle and (not jHasCycle):
        return True
    # j and k are both in a cycle, need to check if one can get from j to k 
    queue = []
    for a in j.bonds: 
        if a!= k and a.ffType[2]==reso:
            queue.append((j,a,1))
    while len(queue) > 0:
        (a,b,depth) = queue.pop()
        if b == k:
            return False # we have found a loop from j to k
        if depth <= 4: # depth <=4 would be enought for look for loop from j to k, but <=5 is necessary for looking for loop around j
            for c in b.bonds:
                if c!=a and c.ffType[2]==reso and c!=j:
                    queue.append((b,c,depth+1))
    # we did not find a loop from j to k
    return True

def getTorsionsDreiding(i,j,k,l,FF,haveCycle):
    torsionTypeListRaw = FF.getTorsionType(i.ffType,j.ffType,k.ffType,l.ffType,skipMissing=FF.skipMissing) # this is a list
    if torsionTypeListRaw == None:
        return []
    # for dreiding we need to divide the torsion by the number of torsions
    # this is factor for torsion multiplicity
    fact1 = len(j.bonds) - 1
    fact2 = len(k.bonds) - 1
    factor = 1.0/float(fact1*fact2)
    # this is factor for exocyclic torsions
    if isExocyclic(i,j,k,l,haveCycle):
        factor = 0.4 * factor  # HARDCODED 0.4 FACTOR, hardcoded 0.4 factor !!!!!!!
    torsionTypeList = [t.dreiding(factor) for t in torsionTypeListRaw]

    torsionList = []
    for ii in range(len(torsionTypeList)):
        t = torsionTypeList[ii]
        torsionList.append( Torsion(i,j,k,l,t) )
    return torsionList

def haveCommonNeighbor(i,l):
    for a in i.bonds:
        if a in l.bonds:
            return True
    return False 

def getTorsionsAmber(i,j,k,l,FF,set14):
    torsionTypeList = FF.getTorsionType(i.ffType,j.ffType,k.ffType,l.ffType,skipMissing=FF.skipMissing) # this is a list 
    if torsionTypeList == None:
        return []
    torsionList = []
    if len(torsionTypeList) > 0:
        t = torsionTypeList[0]
        if (i,l) in set14 or (l,i) in set14:  avoid14 = True
        elif i in l.bonds:                    avoid14 = True # 4 cycle ring
        elif haveCommonNeighbor(i,l):         avoid14 = True # 5 cycle ring
        else:                                 avoid14 = False
        set14.add((i,l))
        torsionList.append( Torsion(i,j,k,l,t,avoid14=avoid14) )
    for ii in range(1,len(torsionTypeList)):
        t = torsionTypeList[ii]
        torsionList.append( Torsion(i,j,k,l,t,avoid14=True) )
    return torsionList

def populateTorsions(al,bonds,FF,exceptResidues):
    torsions = []

    if FF.forcefield == 'amber':
        set14 = set()  # keep track of which 1-4 interactions are already included, so that we get right the sign on list of atoms for amber torsion
        for b in bonds:
            j = b.i
            k = b.j
            if (not j.rName.strip() in exceptResidues) and (not k.rName.strip() in exceptResidues):
                for i in j.bonds:
                    for l in k.bonds:
                        if i!=k and l!=j and i!=l:
                            torsions.extend(getTorsionsAmber(i,j,k,l,FF,set14))
                            

    elif FF.forcefield == 'dreiding':
        haveCycle = getCycles(al) # list of atoms of X_R type which are in some cycle
        for b in bonds:
            j = b.i
            k = b.j
            if not j.rName.strip() in exceptResidues:
                for i in j.bonds:
                    for l in k.bonds:
                        if i!=k and l!=j and i!=l:
                            torsions.extend(getTorsionsDreiding(i,j,k,l,FF,haveCycle))        
    else:
        print >> sys.stderr, 'Unknown type of forcefield for torsion parsing: %s. This is not supported (only dreiding and amber). ' % FF.forcefield
        sys.exit(1)
    return torsions

def getInversions(i,j,k,l,possibleInversion):
    inv1 = Inversion(i,j,k,l,possibleInversion)
    inv2 = Inversion(i,k,l,j,possibleInversion)
    inv3 = Inversion(i,l,j,k,possibleInversion)
    return [inv1,inv2,inv3]

def populateInversions(AL,FF):
    inversions = []
    for a in AL:
        if len(a.bonds) == 3:
            i = a
            j = a.bonds[0]
            k = a.bonds[1]
            l = a.bonds[2]
            possibleInversion = FF.getInversionType(i.ffType,j.ffType,k.ffType,l.ffType)
            if not possibleInversion is None:
                inversions.extend(getInversions(i,j,k,l,possibleInversion))                
    return inversions

############# helping functions ##############

def dictionaryToList(myDict):
    myList = range(len(myDict))
    for i,j in myDict.iteritems():
        myList[j] = i
    return myList

def getUsedAtomLabels(AL):
    usedAtomLabels = {}
    for a in AL:
        if a.ffType not in usedAtomLabels:
            number =  len(usedAtomLabels) 
            usedAtomLabels[a.ffType] = number
    return usedAtomLabels

def getUsedTypes(dataWithType):
    usedTypes = {}
    for b in dataWithType:
        if b.type not in usedTypes:
            number =  len(usedTypes) 
            usedTypes[b.type] = number
    return usedTypes

def getUsedHydrogenBondTypes(hydrogenBondTypes, usedAtomLabels):
    usedHydrogenBondTypes = []
    for b in hydrogenBondTypes:
        if (b.donor in usedAtomLabels) and (b.hydrogen in usedAtomLabels) and (b.acceptor in usedAtomLabels):
            usedHydrogenBondTypes.append(b)
    return usedHydrogenBondTypes


