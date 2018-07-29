"""modul for reading and writing amber files SEE IMPROPER WARNING

note1: the exclusion list in amber prmtop is not complete is if there are bonds without torsions (as in a crystal) 


note2: that the order of atoms that come out from leap matters.
official rule:

                                 J
                                 |
                                 K
                                / \
                               I L

                          Improper I-J-K-L 

where the central atom (K) is the third atom in the
improper
              and the order of the other three is determined
alphabetically
              by atom type and if types are the same by atom number. 

but it is not always enforced
"""

import numpy, sys, time, math, tempfile, subprocess, os
import atoms, structure, pdb, residues

def readTopology(fileName, debug=False):
    """reading topology only from amber prmtop file
    bonds between H1 and H2 in WAT are removed """ 
    if debug: print 'Reading amber topology file %s...' % fileName
    with open(fileName,'r') as f:
        # parsing the file line by line
        # POINTERS
        line = f.readline()
        while not line.startswith('%FLAG POINTERS'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()

        NATOM    = int(f.read(8)) # total number of atoms 
        NTYPES   = int(f.read(8)) # total number of distinct atom types
        NBONH    = int(f.read(8)) # number of bonds containing hydrogen
        MBONA    = int(f.read(8)) # number of bonds not containing hydrogen
        NTHETH   = int(f.read(8)) # number of angles containing hydrogen
        MTHETA   = int(f.read(8)) # number of angles not containing hydrogen
        NPHIH    = int(f.read(8)) # number of dihedrals containing hydrogen
        MPHIA    = int(f.read(8)) # number of dihedrals not containing hydrogen
        NHPARM   = int(f.read(8)) # currently not used
        NPARM    = int(f.read(8)) # used to determine if addles created prmtop
        f.readline()
        NNB      = int(f.read(8)) # number of excluded atoms
        NRES     = int(f.read(8)) # number of residues
        NBONA    = int(f.read(8)) # MBONA + number of constraint bonds
        NTHETA   = int(f.read(8)) # MTHETA + number of constraint angles
        NPHIA    = int(f.read(8)) # MPHIA + number of constraint dihedrals
        NUMBND   = int(f.read(8)) # number of unique bond types
        NUMANG   = int(f.read(8)) # number of unique angle types
        NPTRA    = int(f.read(8)) # number of unique dihedral types
        NATYP    = int(f.read(8)) # number of atom types in parameter file, see SOLTY below
        NPHB     = int(f.read(8)) # number of distinct 10-12 hydrogen bond pair types
        f.readline()
        IFPERT   = int(f.read(8)) # set to 1 if perturbation info is to be read in
        NBPER    = int(f.read(8)) # number of bonds to be perturbed
        NGPER    = int(f.read(8)) # number of angles to be perturbed
        NDPER    = int(f.read(8)) # number of dihedrals to be perturbed
        MBPER    = int(f.read(8)) # number of bonds with atoms completely in perturbed group
        MGPER    = int(f.read(8)) # number of angles with atoms completely in perturbed group
        MDPER    = int(f.read(8)) # number of dihedrals with atoms completely in perturbed groups
        IFBOX    = int(f.read(8)) # set to 1 if standard periodic box, 2 when truncated octahedral
        NMXRS    = int(f.read(8)) # number of atoms in the largest residue
        IFCAP    = int(f.read(8)) # set to 1 if the CAP option from edit was specified
        f.readline()
        NUMEXTRA = int(f.read(8)) # number of extra points found in topology
        # NCOPY    = int(f.read(8)) # number of PIMD slices / number of beads
        
        # ATOM_NAME
        AL = [atoms.Atom() for i in range(NATOM)]
        for i in range(len(AL)):
            AL[i].aNo = i+1
        f.readline()
        f.readline()
        f.readline()
        for i in range(NATOM):
            AL[i].aName = f.read(4).strip()
            if (i+1) % 20 == 0: f.readline()
        # CHARGE
        line = f.readline()
        while not line.startswith('%FLAG CHARGE'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()
        for i in range(NATOM):
            AL[i].charge = float(f.read(16)) / 18.2223
            if (i+1) % 5 == 0: f.readline()
        # MASS 
        # ATOM_TYPE_INDEX 
        # NUMBER_EXCLUDED_ATOMS 
        # NONBONDED_PARM_INDEX  
        # RESIDUE_LABEL
        while not line.startswith('%FLAG RESIDUE_LABEL'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()
        resLabel = ['' for i in range(NRES)]
        for i in range(NRES):
            resLabel[i] = f.read(4)
            if (i+1) % 20 == 0: f.readline()
        # RESIDUE_POINTER
        while not line.startswith('%FLAG RESIDUE_POINTER'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()
        resPointer = [0 for i in range(NRES)]
        for i in range(NRES):
            resPointer[i] = int(f.read(8))
            if (i+1) % 10 == 0: f.readline()
        resPointer.append(1+NATOM)
        # now fill the residue info
        for i in range(NRES):
            for j in range(resPointer[i],resPointer[i+1]):
                AL[j-1].rName = resLabel[i].strip()
                AL[j-1].rNo   = i+1
        # BOND_FORCE_CONSTANT
        # BOND_EQUIL_VALUE
        # ANGLE_FORCE_CONSTANT
        # ANGLE_EQUIL_VALUE
        # DIHEDRAL_FORCE_CONSTANT
        # DIHEDRAL_PERIODICITY
        # DIHEDRAL_PHASE
        # SCEE_SCALE_FACTOR -- not even present in prmtop
        # SCNB_SCALE_FACTOR -- not even present in prmtop
        # SOLTY -- currently unused (reserved for future use)
        # LENNARD_JONES_ACOEF
        # LENNARD_JONES_BCOEF
        # BONDS_INC_HYDROGEN
        while not line.startswith('%FLAG BONDS_INC_HYDROGEN'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()
        for i in range(NBONH):
            atom1 = int(f.read(8))/3 
            if (3*i+1) % 10 == 0: f.readline()
            atom2 = int(f.read(8))/3 
            if (3*i+2) % 10 == 0: f.readline()
            f.read(8)
            if (3*i+3) % 10 == 0: f.readline()
            AL[atom1].makeBond(AL[atom2])
        # BONDS_WITHOUT_HYDROGEN
        while not line.startswith('%FLAG BONDS_WITHOUT_HYDROGEN'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()
        for i in range(NBONA):
            atom1 = int(f.read(8))/3 
            if (3*i+1) % 10 == 0: f.readline()
            atom2 = int(f.read(8))/3 
            if (3*i+2) % 10 == 0: f.readline()
            f.read(8)
            if (3*i+3) % 10 == 0: f.readline()
            AL[atom1].makeBond(AL[atom2])
        # ANGLES_INC_HYDROGEN
        # ANGLES_WITHOUT_HYDROGEN
        # DIHEDRALS_INC_HYDROGEN
        # DIHEDRALS_WITHOUT_HYDROGEN
        # EXCLUDED_ATOMS_LIST
        # HBOND_ACOEF
        # HBOND_BCOEF
        # HBCUT
        # AMBER_ATOM_TYPE
        while not line.startswith('%FLAG AMBER_ATOM_TYPE'):
            line = f.readline()
            if len(line) == 0: sys.exit(1)
        line = f.readline()
        for i in range(NATOM):
            AL[i].ffType = f.read(4)[0:2]
            if (i+1) % 20 == 0: f.readline()
        # TREE_CHAIN_CLASSIFICATION
        # JOIN_ARRAY
        # IROTAT
        # SOLVENT_POINTERS
        # ATOMS_PER_MOLECULE
        # BOX_DIMENSIONS
        # FLAG RADIUS_SET
        # RADII
        # SCREEN

        if IFCAP > 0:
            print >> sys.stderr, 'IFCAP > 0 is not implemented (IFCAP is %d here )' % IFCAP
            sys.exit(1)
        if IFPERT > 0:
            print >> sys.stderr, 'IFPERT > 0 is not implemented (IFPERT is %d here )' % IFPERT
            sys.exit(1)

    if debug: print '%d atoms read from %s.' %(len(AL), fileName)
    # removing bonds between H1 and H2 in WAT
    for a in AL:
        if a.rName == 'WAT' and a.aName == 'O':
            if len(a.bonds) != 2:
                print 'WARNING: atom %s is WATER oxygen, but is not bonded to two atoms' % a   
            else:
                h1 = a.bonds[0]
                h2 = a.bonds[1]
                if (h1.aName not in ['H1','H2'] ) or (h2.aName not in ['H1','H2'] ):
                    print 'WARNING: atom %s is WATER oxygen, the two atoms it is bonded to are not H1, H2' % a  
                else:
                    h1.deleteBond( h2 )
    return AL

def readCoordinatesXYZ(fileName, debug=False):
    """reading coordinates from amber inpcrd or restart file""" 
    if debug: print 'Reading amber coordinates file %s...' % fileName
    with open(fileName,'r') as f:
        f.readline()
        nAtoms  = int(f.readline())
        xyz = numpy.zeros((nAtoms,3))
        for i in range(nAtoms):
            xyz[i,0] = float(f.read(12))
            xyz[i,1] = float(f.read(12))
            xyz[i,2] = float(f.read(12))
            if (i+1)%2==0: f.readline()
    if debug: print 'Done reading amber coordinates file %s.' % fileName
    return xyz  

def readTopoCord(topologyFile,coordinatesFile,debug=False):
    """reads amber prmtop and inpcrd files at once, returns the produced atom list"""
    AL = readTopology(topologyFile,debug=debug)
    xyz = readCoordinatesXYZ(coordinatesFile, debug=debug)
    atoms.putXyzIntoList(AL,xyz)

    for a in AL:
        if a.rName in residues.name321:
            a.aTag = 'ATOM  '
    return AL

def makeResidueList(AL):
    """helping function for generating prmtop"""
    resList = [1]
    for a in AL[2:]:
        prevResidueAtom = AL[resList[-1]]
        if a.rNo != prevResidueAtom.rNo or a.chain != prevResidueAtom.chain :
            resList.append(a.aNo)
    return resList

def getExcludedAtoms(AL):
    """helping functions for generating prmtop"""
    atoms.sortBondsLists(AL)    
    numberExcludedAtoms = [0] * len(AL)
    excludedAtoms = []
    for i in range(len(AL)):
        a = AL[i]
        excluded = []
        for a12 in a.bonds:
            excluded.append(a12.aNo)
            for a13 in a12.bonds:
                excluded.append(a13.aNo)
                for a14 in a13.bonds:
                    excluded.append(a14.aNo)
        myset = set(excluded)  # this removes duplicates
        excluded = list(myset)
        mynumpy = numpy.array(excluded)
        excluded = list( mynumpy[mynumpy > a.aNo] )    # we only need atoms with larger index than this atom
        excluded.sort()
        # apparently amber stores 0 in the NATEX array if there is no excluded atoms 
        if len(excluded) == 0:
            excluded.append(0)
        excludedAtoms.extend(excluded)
        numberExcludedAtoms[i] = len(excluded)
    return numberExcludedAtoms,excludedAtoms

def isHydrogen(atom):
    """helping functions for generating prmtop"""
    if atom.ffType[0] == 'H' or atom.ffType[0] == 'h':
        return True
    else:
        return False
def hasHydrogenBond(b):
    """helping functions for generating prmtop"""    
    return isHydrogen(b.i) or isHydrogen(b.j)
def hasHydrogenAngle(b):
    """helping functions for generating prmtop"""    
    return isHydrogen(b.i) or isHydrogen(b.k)
def hasHydrogenTorsion(b):
    """helping functions for generating prmtop"""
    return isHydrogen(b.i) or isHydrogen(b.l)
def hasHydrogenInversion(b):
    """helping functions for generating prmtop"""    
    return isHydrogen(b.j) or isHydrogen(b.k) or isHydrogen(b.l)
def nothasHydrogenBond(b):
    """helping functions for generating prmtop"""    
    return not hasHydrogenBond(b)
def nothasHydrogenAngle(b):
    """helping functions for generating prmtop"""    
    return not hasHydrogenAngle(b)
def nothasHydrogenTorsion(b):
    """helping functions for generating prmtop"""    
    return not hasHydrogenTorsion(b)
def nothasHydrogenInversion(b):
    """helping functions for generating prmtop"""    
    return not hasHydrogenInversion(b)

def dToPhase(D):
    """helping functions for generating prmtop"""    
    if D ==-1:  return 0.0
    elif D==1: return 180.0
    else:      print >> sys.stderr, 'Invalid D: %f (should be +/-1)' % D; sys.exit(1)

def writeTopology(fileName,AL,FF,debug=False, exceptResidues=[]):
    """write topology to prmtop file""" 
    if debug: print 'Writing amber topology file %s...' % fileName
    
    # get ready for amber
    ST = structure.Structure(AL,FF,debug=debug, exceptResidues=exceptResidues)
    residueList = makeResidueList(AL)
    numberExcludedAtoms,excludedAtoms = getExcludedAtoms(AL)
    
    bondsWithH = filter(hasHydrogenBond,ST.bonds)
    anglesWithH = filter(hasHydrogenAngle,ST.angles)
    torsionsWithH = filter(hasHydrogenTorsion,ST.torsions)
    inversionsWithH = filter(hasHydrogenInversion,ST.inversions)
    bondsWithoutH = filter(nothasHydrogenBond,ST.bonds)
    anglesWithoutH = filter(nothasHydrogenAngle,ST.angles)
    torsionsWithoutH = filter(nothasHydrogenTorsion,ST.torsions)
    inversionsWithoutH = filter(nothasHydrogenInversion,ST.inversions)

    with open(fileName,'w') as f:
        f.write('%%VERSION  VERSION_STAMP = V0001.000  DATE = %s                  \n' % time.strftime('%x  %X') )
        f.write('%FLAG TITLE                                                                     \n')
        f.write('%FORMAT(20a4)                                                                   \n')
        f.write('title                                                                           \n')                                                       
        f.write('%FLAG POINTERS                                                                  \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        NATOM    = len(AL) # total number of atoms 
        NTYPES   = len(ST.usedAtomLabels) # total number of distinct atom types
        NBONH    = len(bondsWithH)       # number of bonds containing hydrogen
        MBONA    = len(bondsWithoutH)    # number of bonds not containing hydrogen
        NTHETH   = len(anglesWithH)      # number of angles containing hydrogen
        MTHETA   = len(anglesWithoutH)   # number of angles not containing hydrogen
        NPHIH    = len(torsionsWithH)+len(inversionsWithH) # number of dihedrals containing hydrogen
        MPHIA    = len(torsionsWithoutH)+len(inversionsWithoutH) # number of dihedrals not containing hydrogen
        NHPARM   = 0      # currently not used
        NPARM    = 0      # used to determine if addles created prmtop
        NNB      = len(excludedAtoms)  # number of excluded atoms 
        NRES     = len(residueList)  # number of residues
        NBONA    = MBONA  # MBONA + number of constraint bonds
        NTHETA   = MTHETA # MTHETA + number of constraint angles
        NPHIA    = MPHIA  # MPHIA + number of constraint dihedrals
        NUMBND   = len(ST.sortedBondTypes) # number of unique bond types
        NUMANG   = len(ST.sortedAngleTypes) # number of unique angle types
        NPTRA    = len(ST.sortedTorsionTypes)+len(ST.sortedInversionTypes) # number of unique dihedral types
        NATYP    = 1      # number of atom types in parameter file, see SOLTY below
        NPHB     = 1      # number of distinct 10-12 hydrogen bond pair types
        IFPERT   = 0      # set to 1 if perturbation info is to be read in
        NBPER    = 0      # number of bonds to be perturbed
        NGPER    = 0      # number of angles to be perturbed
        NDPER    = 0      # number of dihedrals to be perturbed
        MBPER    = 0      # number of bonds with atoms completely in perturbed group
        MGPER    = 0      # number of angles with atoms completely in perturbed group
        MDPER    = 0      # number of dihedrals with atoms completely in perturbed groups
        IFBOX    = 0      # set to 1 if standard periodic box, 2 when truncated octahedral, 0 when there is no box
        NMXRS    = len(AL) # number of atoms in the largest residue ???????
        IFCAP    = 0      # set to 1 if the CAP option from edit was specified
        NUMEXTRA = 0      # number of extra points found in topology
        # NCOPY    = int(f.read(8)) # number of PIMD slices / number of beads
        f.write('%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (NATOM,  NTYPES, NBONH,  MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM,))
        f.write('%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (NNB,    NRES,   NBONA,  NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,))
        f.write('%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (IFPERT, NBPER,  NGPER,  NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS,  IFCAP,))
        f.write('%8d\n' % (NUMEXTRA))
        f.write('%FLAG ATOM_NAME                                                                 \n')
        f.write('%FORMAT(20a4)                                                                   \n')
        for i in range(len(AL)):
            f.write('%-4.4s' % AL[i].aName)
            if (i+1)%20 == 0: f.write('\n')
        if len(AL) % 20 != 0: f.write('\n')
        f.write('%FLAG CHARGE                                                                    \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(AL)):
            f.write('%16.8E' % (AL[i].charge * 18.2223))
            if (i+1)%5 == 0: f.write('\n')
        if len(AL) % 5 != 0: f.write('\n')
        f.write('%FLAG MASS                                                                      \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(AL)):
            f.write('%16.8E' %  FF.atomTypes[AL[i].ffType].mass  )
            if (i+1)%5 == 0: f.write('\n')
        if len(AL) % 5 != 0: f.write('\n')
        f.write('%FLAG ATOM_TYPE_INDEX                                                           \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        for i in range(len(AL)):
            label = AL[i].ffType
            labelindex = ST.usedAtomLabels[label]
            f.write('%8d' %  (labelindex+1) )
            if (i+1)%10 == 0: f.write('\n')
        if len(AL) % 10 != 0: f.write('\n')
        f.write('%FLAG NUMBER_EXCLUDED_ATOMS                                                     \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        for i in range(len(AL)):
            f.write('%8d' %  numberExcludedAtoms[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(AL) % 10 != 0: f.write('\n')
        f.write('%FLAG NONBONDED_PARM_INDEX                                                      \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        array = [0]*(NTYPES*NTYPES)
        for i in range(1,NTYPES+1):
            for j in range(1,NTYPES+1):
                index = NTYPES * (i-1) + j
                a = max(i,j)
                b = min(i,j)
                array[index-1] = b + a*(a-1)/2
        data = array
        for i in range(len(data)):
            f.write('%8d' % data[i])
            if (i+1)%10 == 0: f.write('\n')
        if len(data) % 10 != 0: f.write('\n')
        f.write('%FLAG RESIDUE_LABEL                                                             \n')
        f.write('%FORMAT(20a4)                                                                   \n')
        for i in range(len(residueList)):
            f.write('%-4s' % AL[residueList[i]-1].rName)
            if (i+1)%20 == 0: f.write('\n')
        if len(residueList) % 20 != 0: f.write('\n')
        f.write('%FLAG RESIDUE_POINTER                                                           \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        for i in range(len(residueList)):
            f.write('%8d' % residueList[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(residueList) % 10 != 0: f.write('\n')
        f.write('%FLAG BOND_FORCE_CONSTANT                                                       \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(ST.sortedBondTypes)):
            f.write('%16.8E' %  (ST.sortedBondTypes[i].K / 2.0))
            if (i+1)%5 == 0: f.write('\n')
        if len(ST.sortedBondTypes) % 5 !=0 or len(ST.sortedBondTypes)==0: f.write('\n')
        f.write('%FLAG BOND_EQUIL_VALUE                                                          \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(ST.sortedBondTypes)):
            f.write('%16.8E' %  ST.sortedBondTypes[i].R )
            if (i+1)%5 == 0: f.write('\n')
        if len(ST.sortedBondTypes) % 5 != 0 or len(ST.sortedBondTypes)==0: f.write('\n')
        f.write('%FLAG ANGLE_FORCE_CONSTANT                                                      \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(ST.sortedAngleTypes)):
            f.write('%16.8E' %  (ST.sortedAngleTypes[i].K / 2.0))
            if (i+1)%5 == 0: f.write('\n')
        if len(ST.sortedAngleTypes) % 5 != 0 or len(ST.sortedAngleTypes)==0: f.write('\n')
        f.write('%FLAG ANGLE_EQUIL_VALUE                                                         \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(ST.sortedAngleTypes)):
            f.write('%16.8E' %  (ST.sortedAngleTypes[i].theta0 /180.0*math.pi ))
            if (i+1)%5 == 0: f.write('\n')
        if len(ST.sortedAngleTypes) % 5 != 0 or len(ST.sortedAngleTypes)==0: f.write('\n')
        f.write('%FLAG DIHEDRAL_FORCE_CONSTANT                                                   \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        numTor = len(ST.sortedTorsionTypes)
        numInv = len(ST.sortedInversionTypes)
        for i in range(numTor):
            f.write('%16.8E' %  (ST.sortedTorsionTypes[i].K/2.0 ))
            if (i+1)%5 == 0: f.write('\n')
        for i in range(numInv):
            f.write('%16.8E' %  (ST.sortedInversionTypes[i].amber()) )
            if (i+1+numTor)%5 == 0: f.write('\n')
        if (numTor+numInv) % 5 != 0 or (numTor+numInv)==0: f.write('\n')
        f.write('%FLAG DIHEDRAL_PERIODICITY                                                      \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        numTor = len(ST.sortedTorsionTypes)
        numInv = len(ST.sortedInversionTypes)
        for i in range(numTor):
            f.write('%16.8E' %  ST.sortedTorsionTypes[i].N )
            if (i+1)%5 == 0: f.write('\n')
        for i in range(numInv):
            f.write('%16.8E' %  2.0 )
            if (i+1+numTor)%5 == 0: f.write('\n')
        if (numTor+numInv) % 5 != 0 or (numTor+numInv)==0: f.write('\n')
        f.write('%FLAG DIHEDRAL_PHASE                                                            \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        numTor = len(ST.sortedTorsionTypes)
        numInv = len(ST.sortedInversionTypes)
        for i in range(numTor):
            dToPhase(ST.sortedTorsionTypes[i].D)
            f.write('%16.8E' %  (dToPhase(ST.sortedTorsionTypes[i].D)/180.0*math.pi ))
            if (i+1)%5 == 0: f.write('\n')
        for i in range(numInv):
            f.write('%16.8E' %  (180.0/180.0*math.pi ))
            if (i+1+numTor)%5 == 0: f.write('\n')
        if (numTor+numInv) % 5 != 0 or (numTor+numInv)==0: f.write('\n')
        f.write('%FLAG SOLTY                                                                     \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        f.write('%16.8E\n' %  0.0 )
        f.write('%FLAG LENNARD_JONES_ACOEF                                                       \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        array = [0]*(NTYPES*(NTYPES+1)/2)
        for a in range(1,NTYPES+1):
            for b in range(1,a+1):
                aT = FF.atomTypes[ST.sortedAtomLabels[a-1]]
                bT = FF.atomTypes[ST.sortedAtomLabels[b-1]]
                index = b + a*(a-1)/2
                array[index-1] = aT.mixWith(bT,FF.mean)
        data = array
        for i in range(len(data)):
            ac = data[i].vdwD * math.pow(data[i].vdwR,12)
            f.write('%16.8E' %  ac)
            if (i+1)%5 == 0: f.write('\n')
        if len(data) % 5 != 0: f.write('\n')
        f.write('%FLAG LENNARD_JONES_BCOEF                                                       \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        for i in range(len(data)):
            ac = 2.0 * data[i].vdwD * math.pow(data[i].vdwR,6)
            f.write('%16.8E' %  ac)
            if (i+1)%5 == 0: f.write('\n')
        if len(data) % 5 != 0: f.write('\n')
        f.write('%FLAG BONDS_INC_HYDROGEN                                                        \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        array = [0] * (3*len(bondsWithH))
        for i in range(len(bondsWithH)):
            b = bondsWithH[i]
            array[3*i] = (b.i.aNo-1) * 3
            array[3*i+1] = (b.j.aNo-1) * 3
            array[3*i+2] = ST.usedBondTypes[b.type]+1
        for i in range(len(array)):
            f.write('%8d' % array[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(array) % 10 != 0 or len(array) == 0: f.write('\n')
        f.write('%FLAG BONDS_WITHOUT_HYDROGEN                                                    \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        array = [0] * (3*len(bondsWithoutH))
        for i in range(len(bondsWithoutH)):
            b = bondsWithoutH[i]
            array[3*i] = (b.i.aNo-1) * 3
            array[3*i+1] = (b.j.aNo-1) * 3
            array[3*i+2] = ST.usedBondTypes[b.type]+1
        for i in range(len(array)):
            f.write('%8d' % array[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(array) % 10 != 0  or len(array) == 0: f.write('\n')
        f.write('%FLAG ANGLES_INC_HYDROGEN                                                       \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        array = [0] * (4*len(anglesWithH))
        for i in range(len(anglesWithH)):
            b = anglesWithH[i]
            array[4*i] = (b.i.aNo-1) * 3
            array[4*i+1] = (b.j.aNo-1) * 3
            array[4*i+2] = (b.k.aNo-1) * 3
            array[4*i+3] = ST.usedAngleTypes[b.type]+1
        for i in range(len(array)):
            f.write('%8d' % array[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(array) % 10 != 0  or len(array) == 0: f.write('\n')
        f.write('%FLAG ANGLES_WITHOUT_HYDROGEN                                                   \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        array = [0] * (4*len(anglesWithoutH))
        for i in range(len(anglesWithoutH)):
            b = anglesWithoutH[i]
            array[4*i] = (b.i.aNo-1) * 3
            array[4*i+1] = (b.j.aNo-1) * 3
            array[4*i+2] = (b.k.aNo-1) * 3
            array[4*i+3] = ST.usedAngleTypes[b.type]+1
        for i in range(len(array)):
            f.write('%8d' % array[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(array) % 10 != 0  or len(array) == 0: f.write('\n')
        f.write('%FLAG DIHEDRALS_INC_HYDROGEN                                                    \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        # there is minus on 3rd or 3rd and 4th atoms... 
        #    minus on the fourth atom means inversion
        #    minus on the third atom meands not to include this 1-4 interaction
        lenWithH = len(torsionsWithH)
        array = [0] * (5*lenWithH+5*len(inversionsWithH))
        for i in range(lenWithH):
            b = torsionsWithH[i]
            array[5*i] = (b.i.aNo-1) * 3
            array[5*i+1] = (b.j.aNo-1) * 3
            if b.avoid14:  avoid14 = -1
            else:          avoid14 = 1
            array[5*i+2] = (b.k.aNo-1) * 3 * avoid14
            array[5*i+3] = (b.l.aNo-1) * 3
            array[5*i+4] = ST.usedTorsionTypes[b.type]+1
        for i in range(lenWithH,lenWithH+len(inversionsWithH)):
            b = inversionsWithH[i-lenWithH]
            array[5*i] = (b.j.aNo-1) * 3
            array[5*i+1] = (b.k.aNo-1) * 3
            array[5*i+2] = -(b.i.aNo-1) * 3
            array[5*i+3] = -(b.l.aNo-1) * 3
            array[5*i+4] = ST.usedInversionTypes[b.type]+1+numTor
        for i in range(len(array)):
            f.write('%8d' % array[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(array) % 10 != 0  or len(array) == 0: f.write('\n')        
        f.write('%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        lenWithH = len(torsionsWithoutH)
        array = [0] * (5*lenWithH+5*len(inversionsWithoutH))
        for i in range(lenWithH):
            b = torsionsWithoutH[i]
            array[5*i] = (b.i.aNo-1) * 3
            array[5*i+1] = (b.j.aNo-1) * 3
            if b.avoid14:  avoid14 = -1
            else:          avoid14 = 1            
            array[5*i+2] = (b.k.aNo-1) * 3 * avoid14
            array[5*i+3] = (b.l.aNo-1) * 3
            array[5*i+4] = ST.usedTorsionTypes[b.type]+1
        for i in range(lenWithH,lenWithH+len(inversionsWithoutH)):
            b = inversionsWithoutH[i-lenWithH]
            array[5*i] = (b.j.aNo-1) * 3
            array[5*i+1] = (b.k.aNo-1) * 3
            array[5*i+2] = -(b.i.aNo-1) * 3
            array[5*i+3] = -(b.l.aNo-1) * 3
            array[5*i+4] = ST.usedInversionTypes[b.type]+1+numTor
        for i in range(len(array)):
            f.write('%8d' % array[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(array) % 10 != 0  or len(array) == 0: f.write('\n')        
        f.write('%FLAG EXCLUDED_ATOMS_LIST                                                       \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        data = excludedAtoms        
        for i in range(len(data)):
            f.write('%8d' % data[i] )
            if (i+1)%10 == 0: f.write('\n')
        if len(data) % 10 != 0: f.write('\n')
        f.write('%FLAG HBOND_ACOEF                                                               \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        f.write('  0.00000000E+00\n')
        f.write('%FLAG HBOND_BCOEF                                                               \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        f.write('  0.00000000E+00\n')
        f.write('%FLAG HBCUT                                                                     \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        f.write('  0.00000000E+00\n')
        f.write('%FLAG AMBER_ATOM_TYPE                                                           \n')
        f.write('%FORMAT(20a4)                                                                   \n')
        data = AL
        for i in range(len(data)):
            f.write('%-4s' % data[i].ffType )
            if (i+1)%20 == 0: f.write('\n')
        if len(data) % 20 != 0: f.write('\n')
        f.write('%FLAG TREE_CHAIN_CLASSIFICATION                                                 \n')
        f.write('%FORMAT(20a4)                                                                   \n')
        # just filling crap !!!!!!       
        data = AL
        for i in range(len(data)):
            f.write('BLA ' )
            if (i+1)%20 == 0: f.write('\n')
        if len(data) % 20 != 0: f.write('\n')
        f.write('%FLAG JOIN_ARRAY                                                                \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        data = AL
        for i in range(len(data)):
            f.write('       0' )
            if (i+1)%10 == 0: f.write('\n')
        if len(data) % 10 != 0: f.write('\n')
        f.write('%FLAG IROTAT                                                                    \n')
        f.write('%FORMAT(10I8)                                                                   \n')
        data = AL
        for i in range(len(data)):
            f.write('       0' )
            if (i+1)%10 == 0: f.write('\n')
        if len(data) % 10 != 0: f.write('\n')
        #f.write('%FLAG BOX_DIMENSIONS                                                            \n')
        #f.write('%FORMAT(5E16.8)                                                                 \n')
        #f.write('  9.00000000E+01  5.62095570E+01  5.66096690E+01  6.77711270E+01\n')
        f.write('%FLAG RADIUS_SET                                                                \n')
        f.write('%FORMAT(1a80)                                                                   \n')
        f.write('modified Bondi radii (mbondi)                                                   \n')
        f.write('%FLAG RADII                                                                     \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        # just filling crap !!!!!!
        data = AL
        for i in range(len(data)):
            f.write('%16.8E' % 0.0 )
            if (i+1)%5 == 0: f.write('\n')
        if len(data) % 5 != 0: f.write('\n')
        f.write('%FLAG SCREEN                                                                    \n')
        f.write('%FORMAT(5E16.8)                                                                 \n')
        # just filling crap !!!!!!
        data = AL
        for i in range(len(data)):
            f.write('%16.8E' % 0.0 )
            if (i+1)%5 == 0: f.write('\n')

    if debug: print '%d atoms written to %s.' %(len(AL), fileName)


def writeCoordinates(fileName,AL,debug=False):
    """write coorinates to inpcrd file""" 
    if debug: print 'Writing amber coordinates file %s...' % fileName

    with open(fileName,'w') as f:
        f.write('title\n')
        f.write('%6d\n' % len(AL))
        for i in range(len(AL)):
            xyz = AL[i].xyz
            f.write('%12.7F%12.7F%12.7F' % (xyz[0],xyz[1],xyz[2]) )
            if (i+1)%2 == 0: f.write('\n')
        if len(AL) % 2 != 0: f.write('\n')

        # do not write cell information
        #f.write('  56.2095570  56.6096690  67.7711270  90.0000000  90.0000000  90.0000000\n')

    if debug: print '%d atoms written to %s.' %(len(AL), fileName)

def writeTopoCord(fileNameRoot,AL,FF,debug=False, exceptResidues=[]):
    """write topology and coordinates to prmtop and inpcrd file""" 
    writeTopology(fileNameRoot+'.prmtop',AL,FF,debug=debug,exceptResidues=exceptResidues)
    writeCoordinates(fileNameRoot+'.inpcrd',AL,debug=debug)

def convertFromTleap(al):
    """retypes IP ions to Na, IM to Cl, removes third bond in water"""
    for a in al:
        t = a.ffType.strip()
        if t == 'IP':
            a.ffType = 'Na'
        elif t == 'IM':
            a.ffType = 'Cl'
        elif t == 'OW':
            # this is water removing third water bond
            if len(a.bonds) == 2:
                h1 = a.bonds[0]
                h2 = a.bonds[1]
                if len(h1.bonds) == 2 and len(h2.bonds) == 2:
                    h1.deleteBond(h2)

def addHydrogensAndBonds(proteinpdbal):
    """ takes atom list, writes pdb
    runs tleap to add hydrogens and bonds
    I order the atoms in the individual residues, because tleap uses oreder that does not work with scream
    """
    f = tempfile.NamedTemporaryFile(mode='w', prefix='tmpAmber', dir='.')
    f.close()
    filename = f.name
    al = proteinpdbal
    for a in al:
        if a.aTag != 'ATOM  ':
            print 'Warning: addHydrogensAndBonds: will run TLEAP with HETATM entries, so not sure what will happen'
            break   
    pdb.write(filename+'.pdb', al, withConect=False)
    with open(filename+'.tleap','w') as f:
        f.write('source leaprc.ff99SB \n')
        f.write('logFile %s.log \n' % filename)
        f.write('pdb = loadpdb %s.pdb \n' % filename)
        f.write('saveamberparm pdb %s.prmtop %s.inpcrd \n' % (filename,filename))
        f.write('quit \n')    
    subprocess.call(['tleap','-f',filename+'.tleap'])
    # check if wrong atoms have been added
    lines = open(filename+'.log').readlines()
    for line in lines:
        if line[0:5] == 'FATAL':
            sys.exit('Error: fatal error in tleap (probably there are some unknown or extra atoms in the pdb)')
        elif line[0:7] == '  Added':
            a = line.split()[-2].split('<')[-1]
            if a != 'OXT':
                print 'Warning: addHydrogensAndBonds: tleap added more heavy atoms than just OXT: %s in %s ' % (a,' '.join(line.split()[4:]))
    # read the amber files
    al2 = readTopoCord(filename+'.prmtop',filename+'.inpcrd')
    rl = residues.makeList(al)
    rl2 =  residues.makeList(al2)
    if len(rl) != len(rl2):
        sys.exit('Error: addHydrogensAndBonds: the file from tleap has different number of residues %d than the original file %d' % (len(rl2),len(rl)))
    for i,r in enumerate(rl):
        rl2[i].chain = r.chain
        rl2[i].rNo = r.rNo
        rl2[i].updateAtoms()
        rl2[i].sort()
    al2 = residues.makeAtomList(rl2)
    for a in al2:
        a.aTag = 'ATOM  '

    # clean up
    os.remove(filename+'.prmtop')
    os.remove(filename+'.inpcrd')
    os.remove(filename+'.tleap')
    os.remove(filename+'.pdb')
    os.remove(filename+'.log')
    os.remove('leap.log')
    return al2


