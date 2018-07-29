import sys
import copy, math, numpy
import scipy.special.lambertw

class ForceField:
    """store here the force field parameters
    """
    def __init__(self):
        
        self.forcefield     = 'dreiding'  # dreiding or amber
        self.special14lj    = 1.0 # 1.0  # factor for 1-4 vdw
        self.special14coul  = 1.0 # 1.0  # factor for 1-4 coulomb
        self.mean           = 'geometric'  # geometric or arithmetic
        self.dielectric     = 1.0 # 2.5
        self.splineNonBond  = 180.0 # 10
        self.cutoffNonBond  = 200.0 # 12
        self.splineHBond    = 5.0 # 5
        self.cutoffHBond    = 6.0 # 6
        self.angleHBond     = 90.0 # 90

        self.doHbonds = True

        self.skipMissing = True   # by default missing angles and missing torsions are skipped
        
        self.atomTypes = {}       # dictionary 'ffType' -> AtomType
        self.bondTypes = {}       # dictionary 'ffType1ffType2' -> BondType
        self.angleType3 = {}
        self.angleType1 = {}      # bond and angle has to be unique so the dictionary returns the angleType directly
        self.torsionType4 = {}
        self.torsionType3 = {}
        self.torsionType2 = {}    # torsion does not have to be unique so the dictionary returns a list
        self.inversionType4 = {} 
        self.inversionType3 = {} 
        self.inversionType2 = {}  # central atom is first
        self.inversionType1 = {}  # central atom is first
        self.offDiagVDWTypes = {} # dictionary 'ffType1ffType2' -> AtomType
        self.hydrogenBondTypes = [] 

    def getBondType(self, i, j):
        """ find the actual bond type
        parameters are ffType strings
        """
        if i+j in self.bondTypes:
            return self.bondTypes[i+j]
        elif j+i in self.bondTypes:
            return self.bondTypes[j+i]
        else:
            print >> sys.stderr, 'Error: this bond is not in the force field: "%s"-"%s" ' % (i,j); sys.exit(1)

    def getAngleType(self,i,j,k,skipMissing=False):
        if i+j+k in self.angleType3:
            return self.angleType3[i+j+k]
        elif k+j+i in self.angleType3:
            return self.angleType3[k+j+i]
        elif j in self.angleType1:
            return self.angleType1[j]
        else:
            if skipMissing:
                print 'Warning: skipping angle, this angle is not in the force field: %s-%s-%s ' % (i,j,k)
                return None
            else:
                print >> sys.stderr, 'Error: this angle is not in the force field: %s-%s-%s ' % (i,j,k); sys.exit(1)

    def getTorsionType(self,i,j,k,l,skipMissing=False):
        ''' note that these are lists to cover case of multiple torsions'''
        if i+j+k+l in self.torsionType4:
            return self.torsionType4[i+j+k+l]
        elif l+k+j+i in self.torsionType4:
            return self.torsionType4[l+k+j+i]
        elif i+j+k in self.torsionType3:
            return self.torsionType3[i+j+k]
        elif l+k+j in self.torsionType3:
            return self.torsionType3[l+k+j]
        elif j+k in self.torsionType2:
            return self.torsionType2[j+k]
        elif k+j in self.torsionType2:
            return self.torsionType2[k+j]
        else:
            if skipMissing:
                print 'Warning: skipping torsion, this torsion is not in the force field: %s-%s-%s-%s ' % (i,j,k,l)
                return None
            else:
                print >> sys.stderr, 'Error: this torsion is not in the force field: %s-%s-%s-%s ' % (i,j,k,l); sys.exit(1)

    def getInversionType(self,i,j,k,l):
        """central atom is first
        do not require inversion to exist, returns -1 if there is no inversion
        """
        if i+j+k+l in self.inversionType4:
            return self.inversionType4[i+j+k+l]
        elif i+j+l+k in self.inversionType4:
            return self.inversionType4[i+j+l+k]
        elif i+l+j+k in self.inversionType4:
            return self.inversionType4[i+l+j+k]
        elif i+l+k+j in self.inversionType4:
            return self.inversionType4[i+l+k+j]
        elif i+k+j+l in self.inversionType4:
            return self.inversionType4[i+k+j+l]
        elif i+k+l+j in self.inversionType4:
            return self.inversionType4[i+k+l+j]
        elif i+j+k in self.inversionType3:
            return self.inversionType3[i+j+k]
        elif i+k+j in self.inversionType3:
            return self.inversionType3[i+k+j]
        elif i+j+l in self.inversionType3:
            return self.inversionType3[i+j+l]
        elif i+l+j in self.inversionType3:
            return self.inversionType3[i+l+j]
        elif i+k+l in self.inversionType3:
            return self.inversionType3[i+k+l]
        elif i+l+k in self.inversionType3:
            return self.inversionType3[i+l+k]
        elif i+j in self.inversionType2:
            return self.inversionType2[i+j]
        elif i+k in self.inversionType2:
            return self.inversionType2[i+k]
        elif i+l in self.inversionType2:
            return self.inversionType2[i+l]
        elif i in self.inversionType1:
            return self.inversionType1[i]
        else:
            return None

    def printInfo(self):
        print 'Forcefield info:'
        print '    forcefield:      %s' % self.forcefield
        print '    special14lj:     %f' % self.special14lj
        print '    special14coul:   %f' % self.special14coul
        print '    mean:            %s' % self.mean
        print '    dielectric:      %f' % self.dielectric
        print '    splineNonBond:   %f' % self.splineNonBond
        print '    cutoffNonBond:   %f' % self.cutoffNonBond
        print '    cutoffHBond:     %f' % self.cutoffHBond
        print '    angleHBond:      %f' % self.angleHBond
        print '    # of atom types:          %d' % len(self.atomTypes)
        print '    # of bond types:          %d' % len(self.bondTypes)
        print '    # of angle types:         %d' % (len(self.angleType1) + len(self.angleType3)) 
        print '    # of torsion types:       %d' % (len(self.torsionType2) + len(self.torsionType3)+ len(self.torsionType4)) 
        print '    # of inversion types:     %d' % (len(self.inversionType1) + len(self.inversionType2) + len(self.inversionType3)+ len(self.inversionType4)) 
        print '    # of hydrogen bond types: %d' % len(self.hydrogenBondTypes)
        
class AtomType:
    """ van der walls type
    lj12-6:        D * ( rho^-12 - 2 * rho^-6 )             where rho = r/R
    morse:         D * ( exp( -alpha/2*(r/R-1) ) - 1)^2 - D
    exp6:          D / (zeta-6) * ( 6 * exp(zeta*(1-rho)) - zeta * rho^-6 )   =  A*exp(-B*r) - C/r^6
        where rho = r/R
    """
    def __init__(self):
        self.mass = 0.0
        self.vdwType = ''   # lj12-6 or exp6 or morse
        self.vdwR = -100.0  # let -100 be the default,ie. uninitialized
        self.vdwD = -100.0  # let -100 be the default,ie. uninitialized; D is positive
        self.vdwAlpha = 0.0 # alpha for morse, zeta for exp6
    def __eq__(self, b2):
        return self.mass == b2.mass and self.vdwType == b2.vdwType and self.vdwR == b2.vdwR and self.vdwD == b2.vdwD and self.vdwAlpha == b2.vdwAlpha 
    def __hash__(self):
        return hash(self.mass)+hash(self.vdwType)+hash(self.vdwR)+hash(self.vdwD)+hash(self.vdwAlpha)
    def mixWith(self,type2,mixStyle):
        """returns new atom type which is mix of self and type2"""
        a = copy.deepcopy(self)
        if self.vdwType in ['lj12-6','morse']:
            a.vdwD = math.sqrt(self.vdwD * type2.vdwD)        
            a.vdwAlpha = (self.vdwAlpha + type2.vdwAlpha)/2.0
            if mixStyle == 'geometric':     
                a.vdwR = math.sqrt(self.vdwR * type2.vdwR)
            elif mixStyle == 'arithmetic':  
                a.vdwR = (self.vdwR + type2.vdwR)/2.0
            else: 
                sys.exit('Error: Unknown mixing style %s ' % mixStyle)
        elif self.vdwType == 'exp6':
            if mixStyle == 'geometric':   
                """
                #### THE MICING RULES FROM THE PAPER  
                # the mixing rule for EXP6 in the 1990 dreiding paper is defined for A,B,C
                # get A,B,C first
                A1 = self.vdwD / (self.vdwAlpha-6.0) * 6.0 * numpy.exp(self.vdwAlpha)
                C1 = self.vdwD / (self.vdwAlpha-6.0) * self.vdwAlpha * numpy.power(self.vdwR, 6)
                B1 = self.vdwAlpha / self.vdwR
                A2 = type2.vdwD / (type2.vdwAlpha-6.0) * 6.0 * numpy.exp(type2.vdwAlpha)
                C2 = type2.vdwD / (type2.vdwAlpha-6.0) * type2.vdwAlpha * numpy.power(type2.vdwR, 6)
                B2 = type2.vdwAlpha / type2.vdwR
                # the mixing is done here
                A = numpy.power(A1*A2, 0.5)
                C = numpy.power(C1*C2, 0.5)
                B = (B1+B2)/2.0
                # now inverse teh A,B,C back to R,D,zeta
                expression = -1./7.*numpy.power(6.*C/A*numpy.power(B,6),1./7.)
                # using lambert W function
                zeta = -7.0 * scipy.special.lambertw(expression, k=-1, tol=1e-15)
                if abs(float(numpy.imag(zeta))) > 1e-8:
                    sys.exit('Error: got imaginary root from lambert W function in exp6 geometric mixing of VDW parameters: %s' % zeta)
                a.vdwAlpha = float(numpy.real(zeta))
                a.vdwR = a.vdwAlpha / B
                a.vdwD = C * (a.vdwAlpha - 6.) / a.vdwAlpha / numpy.power(a.vdwR, 6)
                """
                ### Vaclav's new mixing rules
                # R = (Ri+Rj)/2
                # exp = (expi+Rj)/2
                # C6 = sqrt(C6i*C6j)
                # get A,B,C first
                C1 = self.vdwD / (self.vdwAlpha-6.0) * self.vdwAlpha * numpy.power(self.vdwR, 6)
                B1 = self.vdwAlpha / self.vdwR
                C2 = type2.vdwD / (type2.vdwAlpha-6.0) * type2.vdwAlpha * numpy.power(type2.vdwR, 6)
                B2 = type2.vdwAlpha / type2.vdwR
                # mixing
                C = numpy.power(C1*C2, 0.5)
                B = (B1+B2)/2.0                
                a.vdwR = (self.vdwR+type2.vdwR)/2.0
                #                 
                a.vdwAlpha = B*a.vdwR
                a.vdwD = C * (a.vdwAlpha - 6.) / a.vdwAlpha / numpy.power(a.vdwR, 6)
            else: 
                sys.exit('Error: Only geometric mixing implemented for exp6, not: %s ' % mixStyle)
        else:
            sys.exit('Error: Unknown vdw style encountered during mixing: %s ' % self.vdwType)
        return a
    def lammps(self):
        if self.vdwType == 'lj12-6':
            epsilon = self.vdwD
            sigma = self.vdwR * math.pow(2.0, -1.0/6.0)
            return '%.10f %.10f' %(epsilon, sigma)
        elif self.vdwType == 'morse':
            D0 = self.vdwD
            alpha = self.vdwAlpha / 2.0 / self.vdwR
            r0 = self.vdwR
            return '%.10f %.10f %.10f' %(D0,alpha,r0)
        elif self.vdwType == 'exp6':
            # lammps expression for buckingham potential is 
            #      A exp(-r/rho) - C / r^6
            A = self.vdwD / (self.vdwAlpha-6.0) * 6.0 * numpy.exp(self.vdwAlpha)
            rho = self.vdwR/self.vdwAlpha
            C = self.vdwD / (self.vdwAlpha-6.0) * self.vdwAlpha * numpy.power(self.vdwR, 6)
            return '%.10f %.10f %.10f' %(A,rho,C)
        else:
            print >> sys.stderr, 'Unknown atom type for lammps (probably unimplemented in forceField.py yet) "%s" ' % self.vdwType
            sys.exit(1)

class BondType:
    """bond type
    harmonic: 1/2 * k * (r-R)^2 
    """
    def __init__(self):
        self.type  = ''     # harmonic
        self.K = 0.0        # kcal/mol/A^2
        self.R = 0.0        # A
    def __eq__(self, b2):
        return self.type == b2.type and self.K == b2.K and self.R == b2.R
    def __hash__(self):
        return hash(self.type)+hash(self.K)+hash(self.R)
    def lammps(self):
        if self.type == 'harmonic':
            return '%f %f' % (self.K/2.0,self.R)
        else:
            print >> sys.stderr, 'Unknown bond type for lammps (probably unimplemented yet) %s ' % self.type
            sys.exit(1)
            
class AngleType:
    """angle type
    harmonic: 1/2 * K *(theta - theta_0)^2  
    cosine:
       when theta_0 == 180:   1/2 * K *(1 + cos(theta))
       when theta_0 != 180:   1/2 * C *(cos(theta) - cos(theta_0))^2
                                   where C = K / sin(theta_0)^2
    """
    def __init__(self):
        self.type = ''     # harmonic or cosine
        self.K = 0.0       # kcal/mol/rad^2  (for theta in radians)
        self.theta0 = 0.0  # in degrees
    def __eq__(self, b2):
        return self.type == b2.type and self.K == b2.K and self.theta0 == b2.theta0
    def __hash__(self):
        return hash(self.type)+hash(self.K)+hash(self.theta0)
    def lammps(self):
        if self.type == 'harmonic':
            return '%f %f' % (self.K/2.0,self.theta0)
        elif self.type == 'cosine':
            if self.theta0 == 180.0:
                C = self.K
                B = 1
                n = 1 
                return 'cosine/periodic %f %d %d' % (C,B,n)
            else:
                theta0 = self.theta0 
                theta0rad = theta0 * math.pi / 180.0 
                K = self.K / math.pow(math.sin(theta0rad),2) / 2.0
                return 'cosine/squared %f %f' % (K,theta0)
            return '%f %f' % (self.K/2.0,self.theta0) 
        else:
            print >> sys.stderr, 'Unknown angle type for lammps (probably unimplemented yet) %s ' % self.type
            sys.exit(1)

class TorsionType:
    """torsion type
    harmonic: 
     1/2*K*(1-D*cos(N*phi))
          note for D=-1 cis is maximum
               for D=1  cis is minimum
    note: in amber do not need to divide by any factor
    amber stores (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
    and PK is already divided by 2

    So: k = 2*pk/idivf
        n = pn
        d = -1 if phase = 0
        d = 1 if phase = 180
    """
    def __init__(self):
        self.type = '' # harmonic
        self.K = 0.0
        self.N = 0
        self.D = 0
    def __eq__(self, b2):
        return self.type == b2.type and self.K == b2.K and self.N == b2.N and self.D == b2.D
    def __hash__(self):
        return hash(self.type)+hash(self.K)+hash(self.N)+hash(self.D)
    def lammps(self):
        if self.type == 'harmonic':  # lammps formula: K*(1+d*cos(n*phi))
            return '%f %d %d' % (self.K/2.0,-1*self.D,self.N)
        else:
            print >> sys.stderr, 'Unknown torsion type for lammps (probably unimplemented yet) %s ' % self.type
            sys.exit(1)
    def dreiding(self, factor):
        new = TorsionType()
        new.type = self.type
        new.K = self.K * factor
        new.N = self.N
        new.D = self.D
        return new

class InversionType:
    """inversion type for atoms IJKL where I is the central atom
    add three inversions with 1/3 of the coefficient
    here I store the coeffient after dividing by 3
 
    dreiding:
        for psi0 != 0:
            E = Kidreiding / sin(psi0) * (cos(psi)-cos(psi0))^2
        for psi0 ==0    (the only implemented case):
            E = Kidreiding (1 - cos(psi))
        psi is angle between plane IJK and bond IL
            
        note: number of wild cards always 3
    amber:
        E = Pkamber * (1+cos(PN * phi - PHASE))   but PHASE =180deg and PN = 2 
          = Pkamber * (1-cos(2 * phi))
        phi is angle between plance ijk and ijl
        
        note: inputted as kijl with 0,1,2 wildcards: kijl, Xijl, XXjl 
    charmm:
        E = Kicharmm * (psi - psi0)^2  but psi0 = 0 always
          = Kicharmm * psi^2
        psi is angle between plane ijk and jkl
        
        note: has 0 or 4 wildcards ijkl or iXXl

    conversion factors:
    Kidreiding = 16 / 3 * Pkamber
    Kidreiding = 8/9 * Kpsicharmm
    Kpsicharmm = 6 * Pkamber
    """
    def __init__(self):
        self.type = ''  # dreiding or amber ('amber' is cvff in lammps)               
        self.K    = 0.0 # already divided by 3
    def __eq__(self, b2):
        return self.type == b2.type and self.K == b2.K
    def __hash__(self):
        return hash(self.type)+hash(self.K)
    def amber(self):
        return self.K/16.0*3.0
    def charmm(self):
        return self.K/9.0*8.0
    def toAmber(self,amberK):
        self.type = 'amber'
        self.K = 16.0/3.0*amberK
    def toCharmm(self,charmmK):
        self.type = 'charmm'
        self.K = 8.0/9.0*charmmK

    def lammps(self):
        if self.type == 'amber':  # lammps formula for cvff style: K*(1+d*cos(n*phi))
            n = 2
            d = -1
            k = self.amber()
            return '%f %d %d' % (k,d,n)
        elif self.type == 'dreiding':
            K = self.K
            omega0 = 0.0
            return '%f %f' % (K,omega0)
        else:
            print >> sys.stderr, 'Unknown inversion type for lammps (probably unimplemented yet): %s ' % self.type
            sys.exit(1)            
            

class HydrogenBondType:
    """ Hydrogen Bond
    hydrogen bond term replaces two-body VDW between donor and acceptor
    the VDW between hydrogen and acceptor, stayes because the hydrogen has reduced radius anyway
    lj12-10:  D * (5*(R/r)^12 - 6*(R/r)^6) * cos(theta)^anglePower 
    morse:    D * (chi^2 - 2 * chi) * cos(theta)^anglePower
      where chi = exp(-scale/2*(r/R-1)) 
      where r is distance: donor-acceptor
      where theta is angle: donor-hydrogen-acceptor
    """
    def __init__(self):
        self.type = '' # lj12-10 or morse
        self.D = 0.0   # kcal/mol
        self.R = 0.0   # A
        self.scale = 0.0 # no units
        self.power = 0
        self.donor = ''
        self.hydrogen = ''
        self.acceptor = ''
    def __eq__(self, b2):
        return self.type == b2.type and self.D == b2.D and self.R == b2.R and self.scale == b2.scale and self.power == b2.power and self.donor == b2.donor and self.hydrogen == b2.hydrogen and self.acceptor == b2.acceptor
    def __hash__(self):
        return hash(self.type)+hash(self.D)+hash(self.R)+hash(self.scale)+hash(self.power)+hash(self.donor)+hash(self.hydrogen)+hash(self.acceptor)
    def lammps(self):
        if self.type == 'morse':
            # D * (exp(-2a(r-r0)) - 2exp(-a(r-r0))) * cos(theta)^n
            if self.D != 0.0:
                D0 = self.D
                alpha = self.scale/2.0/self.R
                r0 = self.R
            else: # case of zero parameter for this hydrogen bond
                D0 = 0.0
                alpha = 1.0
                r0 = 1.0
            return '%f %f %f' % (D0, alpha, r0)            
        elif self.type == 'lj12-10':
            # 4 * epsilon * (5*(R/r)^12 - 6*(R/r)^6) * cos(theta)^n    FROM DOCUMENTATION BUT WRONG
            # epsilon * ((R/r)^12 - (R/r)^6) * cos(theta)^n    FROM DOCUMENTATION BUT WRONG
            epsilon = self.D   #*6.0*math.pow(6.0/5.0,5)
            r = self.R # * math.sqrt(5.0/6.0)
            return '%f %f' % (epsilon, r)
        else:
            print >> sys.stderr, 'Unknown hydrogen bonds type for lammps (probably unimplemented yet) %s ' % self.type
            sys.exit(1)  

############### functions #################
def fixAtomLabels(AL,FF):
    """fixes the length of ffType stings for teh atoms
    dreiding has 5 character atom types
    amber has 2 character atom types
    """
    if FF.forcefield == 'dreiding':
        format = '%-5.5s'
    elif FF.forcefield == 'amber':
        format  = '%-2.2s'
    else:
        print >> sys.stderr, 'Unknown force field type: %s' % FF.forcefield; sys.exit(1)
    for a in AL:
        a.ffType = format % a.ffType


