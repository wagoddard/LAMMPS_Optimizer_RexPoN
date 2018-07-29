"""atoms is the main module for molecular mechanics, it's main class is Atom
check out the fields of the Atom class in the __init__, they are use throughout this module
"""

import copy,numpy,math,sys

class Atom:
    """atom object stores information for one atom"""

    def __init__(self):
        """creates a blank atom
        
        Note:
        xyz position is stored in numpy array of length 3.
        fixed is optional
        """
        # atom attributes
        # initialized to default values:
        self.aTag   = 'HETATM'# atom tag, 'ATOM  ' or 'HETATM'
        self.aNo    = 0       # atom number
        self.aName  = ''      # atom name (cannot contain spaces)
        self.rName  = ''      # residue name
        self.chain  = ''      # chain name
        self.rNo    = 0       # residue number
        self.xyz    = numpy.zeros(3)# numpy array with x,y,z coordinates
        self.ffType = ''      # forcefield type (can contain spaces at the end)
        self.lpair  = 0       # number of lone pairs
        self.charge = 0.0     # atom charge
        self.fixed  = False   # fixed for dynamics
        self.occupancy = 1.0  # pdb occupancy
        self.beta   = 0.0     # pdb beta factor

        self.bonds  = []      # list of connected atoms
        self.order  = []      # bond order for connected atoms, amide bond is 1.3, aromatic bond is 1.5

    def __str__(self):
        """Short VMD style atom name"""
        return "%s%d:%s" % (self.rName, self.rNo, self.aName) 

    def copy(self):
        """returns a (deep)copy of the atom object """
        return copy.deepcopy(self)

    def copyxyzFrom(self, atom2):
        """(deep)copies coordinates from atom2 """
        self.xyz = atom2.xyz.copy()

    def delete(self):
        """deletes all bonds to this atom """
        for a in self.bonds[:]:
            self.deleteBond(a)        

    # Bond manipulation functions.
    def isBonded(self, atom2):
        """returns true if atom2 is in the bond list of self """
        if atom2 in self.bonds:       return True
        else:                         return False

    def isBondedSafe(self, atom2):
        """returns true if atom2 is in the bond list of self """
        t1 = self.isBonded(atom2)
        t2 = atom2.isBonded(self)
        if t1 and t2: return True
        elif (not t1) and (not t2): return False
        else: 
            print >> sys.stderr, 'Bond list is not consistent %s - %s' % (self, atom2)
            sys.exit(1)    

    def makeBond(self, atom2, bondOrder = 1):
        """makes a bond between self and atom2

        conect information of both atoms are updated;
        bondorder can be set as parameter, default is 1;
        does not check for existing bond
        """
        if atom2 == self:
            print >> sys.stderr, 'Can\'t make bond to self! %s - %s' % (self, atom2)
            sys.exit(1)

        self.bonds.append(atom2)
        self.order.append(bondOrder)
        atom2.bonds.append(self)
        atom2.order.append(bondOrder)

    def makeBondSafe(self, atom2, bondOrder = 1):
        if atom2 in self.bonds or self in atom2.bonds:
            print >> sys.stderr, 'Cannot make bond - atoms are already bonded: %s - %s' % (self, atom2)
            sys.exit(1)
        self.makeBond(atom2,bondOrder = bondOrder)

    def deleteBond(self, atom2):
        """deletes bond between self and atom2. if bond does not exist throws error."""
        if atom2 == self:
            print >> sys.stderr, 'Can\'t delete bond to self! %s - %s' % (self, atom2)
            sys.exit(1)
        if not self.isBonded(atom2):
            print >> sys.stderr, 'Atoms are not bonded -- cannot delete bond! %s - %s' % (self, atom2)
            sys.exit(1)
        self.order.pop( self.bonds.index(atom2) )
        self.bonds.remove(atom2)
        atom2.order.pop( atom2.bonds.index(self) )
        atom2.bonds.remove(self)

    def getBondOrder(self, atom2):
        if atom2 in self.bonds:
            return self.order[ self.bonds.index(atom2) ]
        else:
            sys.exit('Error: atoms %s and %s are not bonded so cannot get the bond order' % (self, atom2))
    
    # Functions on atom and residue names
    def atomName(self):
        """returns striped atom name"""
        return self.aName.strip()

    def resID(self):
        """returns the unique residue ID {chain}{res.no.} of the atom as string for sorting"""
        return "%s%5d" % (self.chain, self.rNo)

    # Geometry functions
    def distance(self, atom2):
        """returns distance to atom"""
        return numpy.linalg.norm(self.xyz - atom2.xyz)

    def within(self, atom2, epsilon):
        """returns True is atom2 is within epsilon of self"""
        if self.distance(atom2) <= epsilon:  return True
        else:                                return False

    def closestAtom(self, atom_list):
        """returns closest atom from a list to self, different from self"""
        return min(atom_list, key = self.distance)

    # organization functions
    def sortBondsList(self):
        """ sort the atoms in the atom.bonds list by the atom number """
        pairs = zip([ai.aNo for ai in self.bonds], self.bonds, self.order)
        pairs.sort()
        self.bonds = [ x[1] for x in pairs]
        self.order = [ x[2] for x in pairs]

    # browsing neighbor lists
    def lessThanThreeBondsFrom(self,b):
        """returns True if self is b, bonded to b, or at most two bonds away from b"""
        if self == b:
            return True
        for c in self.bonds:
            if c == b:
                return True
            for d in c.bonds:
                if d == b:
                    return True
        return False

    def neighborWithName(self,name):
        """returns the first atom with the requested name from hte bond list, if not found returns -1"""
        for a in self.bonds:
            if a.aName == name:
                return a
        return -1

################# end of class atom ##################

def renumberAtomsFrom(AL,fromNumber=1):
    for i in range(len(AL)):
        AL[i].aNo = i+fromNumber

def sortBondsLists(AL):
    for atom in AL:
        atom.sortBondsList()

def charge(AL):
    charge = 0.0
    for a in AL: charge = charge+a.charge
    return charge

def disconnectExternalBonds(al):
    '''deletes bonds to atoms not contained in the atom list'''
    for a in al:
        for b in a.bonds[:]:
            if b not in al:
                a.deleteBond(b) 
def delete(al):
    '''deletes all bonds for atoms in the atom list'''
    for a in al:
        a.delete()

################# manipulation with arrays of coordinates ##############

def putXyzIntoList(AL, xyz):
    """puts coordinates from nAtom x 3 array xyz into the atoms of the atom list AL; does not make copy"""
    if len(AL) != xyz.size/3: 
        print >> sys.stderr, 'atom list and coordinate array have different dimensions'
        sys.exit(1)
    for i in range(len(AL)):
        AL[i].xyz = xyz[i,:]

def getXyzFromList(AL):
    """makes an array nAtom x 3 with coordinates of atoms in the atom list AL; keeps the coordinates and the AL linked"""
    nAtom = len(AL)
    xyz = numpy.zeros((nAtom,3))
    for i in range(nAtom):
        xyz[i,:] = AL[i].xyz.copy()    
    putXyzIntoList(AL, xyz)
    return xyz 

################# MATH functions #####################

def normalized(v):
    """returns a unit vector in the same direction as v, or 0 if the length of vector is 0"""
    n = numpy.linalg.norm(v)
    if n == 0.0:
        return v.copy()
    else:
        return v/n    
    
def angle(v1, v2, degrees = False):
    """Calculates the angle between two vectors.  Always between 0 and pi or 180deg. """
    v1v2dot = numpy.dot(v1,v2)
    v1v2norm = math.sqrt( numpy.dot(v1,v1) * numpy.dot(v2,v2) )
    angle = math.acos(v1v2dot / v1v2norm)
    if degrees:  return angle * 180.0 / math.pi
    else:        return angle

def atomAngle(a1, a2, a3, degrees = False):
    """Calculates the angle a1-a2-a3.  Always between 0 and pi or 180deg."""
    v1 = a1.xyz-a2.xyz
    v2 = a1.xyz-a2.xyz
    return angle(v1, v2, degrees)

def dihedral(a1, a2, a3, a4, degrees = False):
    """Calculates the dihedral angle a1-a2-a3-a4. Between -180deg and 180deg"""
    c3 = normalized(a4.xyz-a3.xyz)
    c2 = normalized(a3.xyz-a2.xyz)
    c1 = normalized(a2.xyz-a1.xyz)
    angle = math.atan2( numpy.dot(c1, numpy.cross(c2,c3)) , numpy.dot(numpy.cross(c1,c2) , numpy.cross(c2,c3)) )
    if degrees:  return angle * 180.0 / math.pi
    else:        return angle

def rad2deg(radians):
    return 180.0*(radians/math.pi)

def deg2rad(degrees):
    return degrees*math.pi/180.0


