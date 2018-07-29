import sys
import forceField, defaults

def readAmberDat(FF, fileName = '', debug=False):
    """reading amber dat file""" 
    if fileName == '':
        fileName = defaults.amberDat

    if debug: print 'Reading amber dat file %s...' % fileName
    with open(fileName,'r') as f:
        FF.forcefield     = 'amber'  # dreiding
        FF.special14lj    = 0.5 # 1.0  # factor for 1-4 vdw
        FF.special14coul  = 1.0/1.2 # 1.0  # factor for 1-4 coulomb
        FF.mean           = 'arithmetic' # geometric for dreiding
        FF.dielectric     = 1.0 # 2.5
        FF.splineNonBond  = 180.0 # 10
        FF.cutoffNonBond  = 200.0 # 12
        FF.splineHBond    = 0.0 # 5
        FF.cutoffHBond    = 0.0 # 6
        FF.angleHBond     = 0.0 # 90

        FF.doHbonds = False

        # read atom types
        line = f.readline()
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberMass(line,FF)
            line = f.readline()

        # read bonds                
        line = f.readline()
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberBond(line,FF)
            line = f.readline()

        # read angles
        # no wild card X are recognized for amber angle term
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberAngle(line,FF)
            line = f.readline()

        # read torsions
        line = f.readline()
        prevList = []   # helping list for adding multiple torsions
        while line[0] != ' ' and line[0] != '\n':
            prevList = readAmberTorsion(line,prevList,FF)
            line = f.readline()

        # read inversions
        # first one or first two atoms can be wild cards
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberInversion(line,FF)
            line = f.readline()
            
        # read hydrogen bond terms
        # if encounter nonzero value report error
        line = f.readline()
        while line[3:4] != ' ' and line[3:4] != '' and line[3:4] != '\n':
            KT1 = line[2:4]
            KT2 = line[6:8]
            A   = float(line[10:20])
            B   = float(line[20:30])
            if A!=0.0 or B!=0.0:
                print >> sys.stderr, 'Read nonzero LJ 12-10 term, this is not supported. '
                sys.exit(1)
            line = f.readline()                
        # reading equivalent atom types
        equivalentAtomTypes = []
        line = f.readline()
        spl = line.split()
        while len(spl) > 0:
            newTypes = []
            for newType in spl:
                newTypes.append('%-2.2s' % newType)
            equivalentAtomTypes.append(newTypes)
            line = f.readline()
            spl = line.split()
        # reading VDW parameters
        line = f.readline()
        if not line.startswith('MOD4      RE'):
            print >> sys.stderr, 'Only RE potential is supported. Not: %s' % line
            sys.exit(1)
        line = f.readline()
        while  line[2:3] != ' ' and line[2:3] != '\n' and line[2:3] != '':
            readAmberVDW(line,FF)
            line = f.readline()
        # refill the info to the equivalent atoms
        for eqType in equivalentAtomTypes:
            eqHead = FF.atomTypes[eqType[0]]
            for i in range(1,len(eqType)):
                thisType = eqType[i]
                FF.atomTypes[thisType].vdwR = eqHead.vdwR
                FF.atomTypes[thisType].vdwD = eqHead.vdwD
        # check that all atoms have been assigned radius and well depth
        for aType in FF.atomTypes:
            if FF.atomTypes[aType].vdwR == -100.0 or FF.atomTypes[aType].vdwD == -100.0:
                print >> sys.stderr, 'This atoms does not have defined vdw radius or well depth: %s' % aType
                sys.exit(1)
        # we do not need to read the rest 
        
    if debug: 
        FF.printInfo()
        print 'Done reading amber dat file %s.' % fileName

def readAmberMass(line,FF):
    new = forceField.AtomType()            
    new.vdwType = 'lj12-6'
    new.mass = float(line[2:13])

    label = line[0:2]
    if label in FF.atomTypes:
        print >> sys.stderr, 'This atom is already in the forcefield: %s' % label; sys.exit(1)
    else:
        FF.atomTypes[label] = new                
        
def readAmberBond(line,FF):
    new = forceField.BondType()
    new.type = 'harmonic'            
    new.K = float(line[5:15]) * 2.0   # dat stores K from K*(r-r0)^2
    new.R = float(line[15:25])

    label = line[0:2]+line[3:5]
    label2 = line[3:5]+line[0:2]
    if label in FF.bondTypes or label2 in FF.bondTypes:
        print >> sys.stderr, 'This bond is already in the forcefield: %s' % label
        sys.exit(1)
    else:
        FF.bondTypes[label] = new 

def readAmberAngle(line,FF):
    new = forceField.AngleType()
    new.type = 'harmonic'
    new.K = float(line[8:18]) * 2.0   # dat stores K from K*(a-a0)^2
    new.theta0 = float(line[18:28])

    label = line[0:2]+line[3:5]+line[6:8]            
    label2 = line[6:8]+line[3:5]+line[0:2]
    if label in FF.angleType3 or label2 in FF.angleType3:
        print >> sys.stderr, 'This angle is already in the forcefield: %s' % label
        sys.exit(1)
    else:
        FF.angleType3[label] = new 

def readAmberTorsion(line,prevList,FF):
    IPT = line[0:2]
    JPT = line[3:5]
    KPT = line[6:8]
    LPT = line[9:11]
    IDIVF = int(line[11:15])
    PK = float(line[15:30])
    PHASE = float(line[30:45])
    PN = float(line[45:60])

    new = forceField.TorsionType()
    new.type = 'harmonic'
    new.K = 2.0*PK/IDIVF
    new.N = abs(PN)
    if PHASE == 0.0:
        new.D = -1
    elif PHASE == 180.0:
        new.D = 1
    else:
        print >> sys.stderr, 'Unknown phase in amber torsion %f (should be 0 or 180)' % PHASE
        sys.exit(1)
    
    prevList.append(new)  # just add the torsion into the temporary list
    if PN > 0: # now we are adding the list into the database
        if IPT == 'X ' and LPT =='X ':
            label = JPT+KPT
            label2 = KPT+JPT
            if label in FF.torsionType2 or label2 in FF.torsionType2:
                print >> sys.stderr, 'This torsion is already in the forcefield: X-%s-X' % label
                sys.exit(1)
            else:
                FF.torsionType2[label] = prevList
        elif IPT != 'X ' and LPT !='X ':
            label = IPT+JPT+KPT+LPT
            label2 = LPT+KPT+JPT+IPT
            if label in FF.torsionType4 or label2 in FF.torsionType4:
                print >> sys.stderr, 'This torsion is already in the forcefield: X-%s-X' % label
                sys.exit(1)
            else:
                FF.torsionType4[label] = prevList
        else:
            print >> sys.stderr, 'This torsion is not with 2 or no wildcards X: %s' % line; sys.exit(1)
        prevList = []                 
    return prevList

def readAmberInversion(line,FF):            
    IPT = line[0:2]
    JPT = line[3:5]
    KPT = line[6:8]
    LPT = line[9:11]
    PK = float(line[15:30])
    PHASE = float(line[30:45])
    PN = float(line[45:60])

    if PHASE != 180.0 or PN != 2:
        print >> sys.stderr, 'Amber inversion has to have phase 180 and PN=2, but has phase = %f and PN= %d ' % (PHASE,PN)
        sys.exit(1)

    new = forceField.InversionType()
    new.toAmber(PK / 3)

    # first one or first two atoms can be wild cards  in amber dat file         
    if IPT == 'X ' and JPT =='X ' and LPT =='X':
        print >> sys.stderr, 'This torsion has too many wildcards: %s' % line; sys.exit(1)
    elif IPT == 'X ' and JPT =='X ':
        label = KPT+LPT
        if label in FF.inversionType2:
            print >> sys.stderr, 'This inversion is already in the forcefield: %s-X-X' % label
            sys.exit(1)
        else:
            FF.inversionType2[label] = new
    elif IPT == 'X ':
        label = KPT+JPT+LPT
        label2 = KPT+LPT+JPT
        if label in FF.inversionType3 or label2 in FF.inversionType3:
            print >> sys.stderr, 'This inversion is already in the forcefield: %s-X' % label
            sys.exit(1)
        else:
            FF.inversionType3[label] = new
    else:
        label  = KPT+IPT+JPT+LPT
        label2 = KPT+IPT+LPT+JPT
        label3 = KPT+JPT+IPT+LPT
        label4 = KPT+JPT+LPT+IPT
        label5 = KPT+LPT+JPT+IPT
        label6 = KPT+LPT+IPT+JPT
        if label  in FF.inversionType4 or label2 in FF.inversionType4 or \
           label3 in FF.inversionType4 or label4 in FF.inversionType4 or \
           label5 in FF.inversionType4 or label6 in FF.inversionType4:
            print >> sys.stderr, 'This inversion is already in the forcefield: %s' % label
            sys.exit(1)
        else:
            FF.inversionType4[label] = new

def readAmberVDW(line,FF):
    label = line[2:4]
    R = float(line[10:22])
    D = float(line[22:34])
    FF.atomTypes[label].vdwR = R*2.0
    FF.atomTypes[label].vdwD = D

def readAmberFrcmod(fileName, FF, debug=False):
    """reading amber frcmod file""" 
    if debug: 
        print 'Reading amber frcmod file %s...' % fileName
        print 'Forcefield before reading frcmod file:'
        FF.printInfo()
        
    with open(fileName,'r') as f:
        # read atom types
        line = f.readline()
        line = f.readline()
        if not line.startswith('MASS'): 
            print >> sys.stderr, 'This should be MASS line, instead got %s ' % line; sys.exit(1)
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberMass(line,FF)
            line = f.readline()

        # read bonds                
        line = f.readline()
        if not line.startswith('BOND'): 
            print >> sys.stderr, 'This should be BOND line, instead got %s ' % line; sys.exit(1)        
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberBond(line,FF)
            line = f.readline()

        # read angles
        # no wild card X are recognized for amber angle term
        line = f.readline()
        if not line.startswith('ANGLE'): 
            print >> sys.stderr, 'This should be ANGLE line, instead got %s ' % line; sys.exit(1)                
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberAngle(line,FF)
            line = f.readline()

        # read torsions
        line = f.readline()
        if not line.startswith('DIHE'): 
            print >> sys.stderr, 'This should be DIHE line, instead got %s ' % line; sys.exit(1)                
        line = f.readline()
        prevList = []   # helping list for adding multiple torsions
        while line[0] != ' ' and line[0] != '\n':
            prevList = readAmberTorsion(line,prevList,FF)
            line = f.readline()

        # read inversions
        # first one or first two atoms can be wild cards
        line = f.readline()
        if not line.startswith('IMPROPER'): 
            print >> sys.stderr, 'This should be IMPROPER line, instead got %s ' % line; sys.exit(1)                
        line = f.readline()
        while line[0] != ' ' and line[0] != '\n':
            readAmberInversion(line,FF)
            line = f.readline()
            
        # reading VDW parameters
        line = f.readline()
        if not line.startswith('NONBON'): 
            print >> sys.stderr, 'This should be NONBON line, instead got %s ' % line; sys.exit(1)                
        line = f.readline()
        while  line[2:3] != ' ' and line[2:3] != '\n' and line[2:3] != '':
            readAmberVDW(line,FF)
            line = f.readline()
        # check that all atoms have been assigned radius and well depth
        for aType in FF.atomTypes:
            if FF.atomTypes[aType].vdwR == -100.0 or FF.atomTypes[aType].vdwD == -100.0:
                print >> sys.stderr, 'This atoms does not have defined vdw radius or well depth: %s' % aType
                sys.exit(1)
        # done
    if debug: 
        print 'Forcefield after reading frcmod file:'
        FF.printInfo()
        print 'Done reading amber dat file %s.' % fileName

