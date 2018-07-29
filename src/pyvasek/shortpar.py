
import sys, re
import forceField, defaults

def label(text):
    '''make dreiding label from a possibly shorter string'''
    return '%-5s' % text 

def nextline(f):
    line = f.readline()
    s = line.split()
    if len(s) == 0 and len(line)>0: s = ['#']
    return s 

def read(ff, fileName = '', debug=False):
    """reading dreiding shortpar file""" 
    if fileName == '':
        fileName = defaults.dreiding3

    if debug: print 'Reading amber dat file %s...' % fileName
    
    with open(fileName,'r') as f:
        s = ['']
        ff.doHbonds = True
        s = nextline(f)
        while s[0] != 'bond' and len(s) > 0:
            if s[0][0] == '#':            1 # dummy statement to do nothing
            elif s[0] == 'forcefield':    ff.forcefield = s[1]
            elif s[0] == 'special14lj':   ff.special14lj = float(s[1])
            elif s[0] == 'special14coul': ff.special14coul = float(s[1])
            elif s[0] == 'mean':          ff.mean = s[1]
            elif s[0] == 'dielectric':    ff.dielectric = float(s[1])
            elif s[0] == 'splineNonBond': ff.splineNonBond = float(s[1])
            elif s[0] == 'cutoffNonBond': ff.cutoffNonBond = float(s[1])
            elif s[0] == 'splineHBond':   ff.splineHBond = float(s[1])
            elif s[0] == 'cutoffHBond':   ff.cutoffHBond = float(s[1])
            elif s[0] == 'angleHBond':    ff.angleHBond = float(s[1])                
            else:
                print >> sys.stderr, 'Unknown header option in the shortpar file %s' % s[0]
                sys.exit(1)                    
            s = nextline(f)
        # read bonds
        # dreiding bonds are computed
        # 1. read the atom radii and atom types            
        s=nextline(f)
        bondRadii = {}                      
        while s[0] != 'angle' and len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               
                    atomLabel = label(s[0]) 
                    if atomLabel in ff.atomTypes:
                        print >> sys.stderr, 'Second bond entry for %s ' % atomLabel
                        sys.exit(1)
                    else:
                        bondRadii[atomLabel] = float(s[1])     
            s=nextline(f)           
        # 2. compute bonds
        atomlist = []
        for b in bondRadii: atomlist.append(b)
        expandBonds(atomlist,bondRadii, ff)
        
        # read angles
        angles = {}
        s=nextline(f)                      
        while s[0] != 'torsionI' and len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               readAngle(atomlist,s,ff)
            s=nextline(f)           
        
        # read torsions
        s=nextline(f)                      
        while s[0] != 'inversionI' and len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               readTorsion(atomlist,s,ff)
            s=nextline(f)            
        
        # read inversion
        s=nextline(f)                        
        while s[0] != 'atomType' and len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               readInversion(atomlist,s,ff)
            s=nextline(f)            
        
        # read atom types
        # they are are expanded bellow the bond section
        s=nextline(f)            
        while s[0] != 'offDiagonalVdw' and len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               readAtom(atomlist,s,ff)
            s=nextline(f)

        # read offDiagonalVDW            
        s=nextline(f)                        
        while s[0] != 'hBondDonor' and len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               readOffDiagVDW(atomlist,s,ff)
            s=nextline(f)            
        
        # read hydrogen bonds
        s=nextline(f)                        
        while len(s) > 0:
            if len(s) > 0:
                if s[0][0] == '#':  pass
                else:               readHBond(atomlist,s,ff)
            s=nextline(f)            
    
    if debug: 
        ff.printInfo()
        print 'Done reading amber dat file %s.' % fileName
        
def expandBonds(atomlist,bondRadii, ff):
    for i in range(len(atomlist)):
        for bj in atomlist[i:]:
            bi = atomlist[i]
            if (bi+bj in ff.bondTypes) or (bj+bi in ff.bondTypes):
                print >> sys.stderr, 'bond %s-%s is already in the forcefield ' % (bi,bj)
                sys.exit(1)
            else:
                ri = bondRadii[bi]
                rj = bondRadii[bj]
                # determine bonds order\
                bondorder = 1.0
                if bi[2] == bj[2] == '1':
                    bondorder = 3.0
                elif bi[2] == bj[2] == 'R':
                    bondorder = 1.5
                elif bi[2] == bj[2] == '2':
                    bondorder = 2.0
                elif (bi[2] == '2' and bj[2] == 'R') or (bj[2] == '2' and bi[2] == 'R'):
                    bondorder = 2.0

                new = forceField.BondType()
                new.type = 'harmonic'
                new.K = 700.0 * bondorder 
                new.R = ri+rj-0.01
                ff.bondTypes[bi+bj] = new

def readAngle(atomlist,s,ff):
    i = label(s[0])
    j = label(s[1])               
    k = label(s[2])
    
    new = forceField.AngleType()
    new.type = s[3]
    new.K = float(s[4]) 
    new.theta0 = float(s[5])

    if not new.type in ['harmonic','cosine']:
        print >> sys.stderr, 'unimplemented inversion type (only know harmonic and cosine): %s' % ' '.join(s)
        sys.exit(1)    
    if i != '.....' or k != '.....':
        print >> sys.stderr, 'only angles with two wildcards implemented for dreiding, not: %s' % ''.join(s)
        sys.exit(1)
    else:
        # matching to atom types
        compile = re.compile(j)
        for a in atomlist:
            if compile.match(a) != None:
                ff.angleType1[a] = new
                
def readTorsion(atomlist,s,ff):
    i = label(s[0])
    j = label(s[1])               
    k = label(s[2])
    l = label(s[3])

    new = forceField.TorsionType()
    new.type = 'harmonic'
    new.K = float(s[4]) 
    new.N = int(s[5])
    new.D = int(s[6])
    
    if l != '.....':
        print >> sys.stderr, 'only torsions with last atom being a wildcard implemented for dreiding not: %s' % ' '.join(s)
        sys.exit(1)
    elif i == '.....': # case of two wildcards
        cj = re.compile(j)
        ck = re.compile(k)
        for aj in atomlist:
            if cj.match(aj) != None:
                for ak in atomlist:
                    if ck.match(ak) != None:
                        ff.torsionType2[aj+ak] = [new]
    else: # case of one wildcard
        ci = re.compile(i)
        cj = re.compile(j)
        ck = re.compile(k)
        for aj in atomlist:
            if cj.match(aj) != None:
                for ak in atomlist:
                    if ck.match(ak) != None:
                        for ai in atomlist:
                            if ci.match(ai) != None:
                                ff.torsionType3[ai+aj+ak] = [new]
                           
def readInversion(atomlist,s,ff):
    i = label(s[0])
    j = label(s[1])               
    k = label(s[2])
    l = label(s[3])

    new = forceField.InversionType()
    new.type = s[4]
    new.K = float(s[5]) / 3.0
    theta0 = float(s[6])

    if theta0 != 0.0: 
        print >> sys.stderr, 'only inversion with theta0 = 0 is implemented, not  %f' % theta0
        sys.exit(1)
    if not new.type in ['dreiding','amber']:
        print >> sys.stderr, 'unimplemented inversion type (only know dreiding and amber): %s' % ' '.join(s)
        sys.exit(1)
    if j != '.....' or k != '.....' or l != '.....' :
        print >> sys.stderr, 'only inversions with three wild cards implemented for dreiding not: %s' % ' '.join(s)
        sys.exit(1)
    else: 
        ci = re.compile(i)
        for ai in atomlist:
            if ci.match(ai) != None:
                ff.inversionType1[ai] = new

def readAtom(atomlist,s,ff):
    i = label(s[0])
    
    new = forceField.AtomType()
    
    new.mass = float(s[1])
    new.vdwType = s[2]
    new.vdwR = float(s[3]) 
    new.vdwD = float(s[4])
    new.vdwAlpha = float(s[5])

    if new.vdwType not in ['lj12-6','morse','exp6']:
        sys.exit('Error: only lj12-6 and morse and exp6  are implemented, not: %s' % ' '.join(s))

    if len(ff.atomTypes) > 0:
        someType = ff.atomTypes.__iter__().next()
        if new.vdwType != ff.atomTypes[someType].vdwType:
            sys.exit('Error: this atom has different VDW type than previous atom (%s): %s' % (ff.atomTypes[someType].vdwType ,' '.join(s)))

    ci = re.compile(i)
    for ai in atomlist:
        if ci.match(ai) != None:
            ff.atomTypes[ai] = new

def readHBond(atomlist, s,ff):
    i = label(s[0])
    j = label(s[1])               
    k = label(s[2])
    
    ntype = s[3]
    nD = float(s[4]) 
    nR = float(s[5])
    nscale = float(s[6])
    npower = int(s[7])
    
    if not ntype in ['morse','lj12-10']:
        print >> sys.stderr, 'unimplemented hydrogen bond type (only know morse and lj12-10): %s' % ' '.join(s)
        sys.exit(1)    
    if j != 'H___A' and j != 'H___W':
        print >> sys.stderr, 'only H___A or H___W are hydrogen bond atoms, not: %s' % ' '.join(s)
        sys.exit(1)
    else: # case of two wildcards
        ci = re.compile(i)
        ck = re.compile(k)
        for ai in atomlist:
            if ci.match(ai) != None:
                for ak in atomlist:
                    if ck.match(ak) != None:
                        new = forceField.HydrogenBondType()
                        new.type = ntype
                        new.D = nD
                        new.R = nR
                        new.scale = nscale
                        new.power = npower
                        new.donor = i
                        new.hydrogen = j
                        new.acceptor = k                        
                        
                        ff.hydrogenBondTypes.append(new)

def readOffDiagVDW(atomlist,s,ff):
    i = label(s[0])
    j = label(s[1])               
    
    new = forceField.AtomType()
    
    new.mass = 0.0
    new.vdwType = s[2]
    new.vdwR = float(s[3]) 
    new.vdwD = float(s[4])
    new.vdwAlpha = float(s[5])

    if new.vdwType not in ['lj12-6','morse', 'exp6']:
        sys.exit('Error: only lj12-6, morse and exp6 are implemented for offdiag VDW, not: %s' % ' '.join(s))

    someType = ff.atomTypes.__iter__().next()
    if new.vdwType != ff.atomTypes[someType].vdwType:
        sys.exit('Error: this offdiagonal VDW rule has different VDW type than previous atom (%s): %s' %(ff.atomTypes[someType].vdwType, ' '.join(s)))

    ci = re.compile(i)
    cj = re.compile(j)
    for ai in atomlist:
        if ci.match(ai) != None:
            for aj in atomlist:
                if cj.match(aj) != None:
                    ff.offDiagVDWTypes[ai+aj] = new    
    
    
