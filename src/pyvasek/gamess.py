""" writing GAMESS input .in and reading details from the output .out 

example: of GAMESS input file 
 
! gamess input file, comments have ! in the first column, comands have $ in the second column
 $contrl scftyp=rhf mplevl=2 ispher=1 runtyp=energy coord=unique units=angs $end
! scftyp = rhf       restricted hartree fock
! MP2: mplevl=2
! B3LYP: dfttyp=B3LYP1
! runtyp=energy or optimize
! icharg = integer charge on molecule, 0 is default
 $basis gbasis=n311 diffsp=.t. diffs=.t. ngauss=6 ndfunc=1 npfunc=1 $end
! 6311G**++     gbasis=n311 diffsp=.t. diffs=.t. ngauss=6 ndfunc=1 npfunc=1
! 631G**++      gbasis=n31 diffsp=.t. diffs=.t. ngauss=6 ndfunc=1 npfunc=1
! 631G**+       gbasis=n31 diffsp=.t. ngauss=6 ndfunc=1 npfunc=1
! 631G*+        gbasis=n31 diffsp=.t. ngauss=6 ndfunc=1
 $system mwords=200  $end
! maximum memory used per node, 1 word is 8MB
 $scf dirscf=.t. $end
! The default value of DIRSCF=.FALSE. causes integrals to be written to disk for later use.
! the line after $data is arbitrary comment line
! then there is Schoenflies symbol for symmetry group, C1 = no symmetry
! the atom label can have up to 10 characters, nuclear charge defines nucleus  
 $data
1.inp
C1
O 8  -4.006870  1.425010  0.085700
H 1  -4.761710  0.855490  -0.062800
H 1  -3.231940  0.873490  -0.022900
O 8  -2.969990  2.508940  0.098290
H 1  -3.724820  1.939420  -0.050220
H 1  -2.195060  1.957420  -0.010320
 $end
"""

import sys
import atoms

element2ztable={'H':1,
                'C':6,
                'O':8,
                'N':7,
                'P':15,
                'S':16,
                'B':5}

def element2z(element):
    if element in element2ztable:
        return element2ztable[element]
    else:
        print >> sys.stderr, 'i do not know what nuclear charge does this element have: %s (probably just have to add it in gamess.py)' % element
        sys.exit(1)

def write(fileName, AL, solvant = False, optimization = False, basis='6-311G**++', method='B3LYP', muliken=False):
    """write input gamess file for structure in atom list AL
    solvant - not implemented yet
    single point energy (default) or optimization
    basis 6-311G**++ or 6-31G**
    method B3LYP or MP2
    muliken - not implemented yet
    """ 
    with open(fileName,'w') as f:
        if solvant or muliken:
            print >> sys.stderr, 'solvant and muliken are not implemented yet for gamess input file'
            sys.exit(1)

        # working on control line
        if method == 'B3LYP':
            methodtext = 'dfttyp=B3LYP1'
        elif method == 'MP2':
            methodtext = 'mplevl=2'
        else:
            print >> sys.stderr, 'only have implemented B3LYP and MP2 but got: %s' % method
            sys.exit(1)
        if optimization:
            runtypetext = 'runtyp=optimize'
        else:
            runtypetext = 'runtyp=energy'
        if abs(atoms.charge(AL)) > 0.001:
            chargetext = 'icharg=%d' % round(atoms.charge(AL))
        else:
            chargetext = 'icharg=0'
        f.write(' $contrl scftyp=rhf %s \n ispher=1 %s %s coord=unique units=angs $end\n' % (methodtext, runtypetext, chargetext))
        
        # working on basis line
        if basis=='6-311G**++':
            f.write(' $basis gbasis=n311 diffsp=.t. diffs=.t. ngauss=6 ndfunc=1 npfunc=1 $end\n')
        elif basis=='6-31G**':
            f.write(' $basis gbasis=n31 ngauss=6 ndfunc=1 npfunc=1 $end\n')
        elif basis=='augCCPVTZ':
            f.write(' $basis gbasis=ACCT $end\n')
        elif basis=='augCCPVQZ':
            f.write(' $basis gbasis=ACCQ $end\n')
        elif basis=='augCCPV5Z':
            f.write(' $basis gbasis=ACC5 $end\n')
        elif basis=='augCCPV6Z':
            f.write(' $basis gbasis=ACC6 $end\n')
        else:
            print >> sys.stderr, 'unplemented basis: %s' % basis
            sys.exit(1)

        # more job control
        f.write(' $system mwords=200  $end\n')
        f.write(' $scf dirscf=.t. $end\n')
        f.write(' $data\n')
        f.write('%s\n' % fileName)
        f.write('C1\n')

        for a in AL:
            if a.fixed:
                print >> sys.stderr, 'atom %s is fixed, but i do not know how to fix atoms in gamess input file ' % str(a)
                sys.exit(1)
            else:
                format = '%-8s %d %20.15f %20.15f %20.15f\n'
            if a.aName.strip()[1:2].islower():
                element = a.aName.strip()[0:2]
            else:
                element = a.aName.strip()[0:1]
            z = element2z(element)
            f.write(format % (a.aName,z,a.xyz[0],a.xyz[1],a.xyz[2]))
        f.write(' $end\n')
