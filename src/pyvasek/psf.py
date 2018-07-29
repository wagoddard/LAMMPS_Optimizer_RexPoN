"""writing psf file to be read by namd
namd does not need column formatted input, but understands space separated input
(in both psf and par file) 
"""

import numpy, sys, time, math
import atoms, structure, prmtop, pdb

class Cont:
    def __init__(self):
        self.atomtypes = {}
        self.bonds = {}
        self.angles = {}
        self.torsions = {} # this contains torsions only
        self.inversions = {} # this contains inversions

def fixAtomNames(AL):
    for a in AL:
        if a.chain.strip() == '':
            a.chain = 'X'
        if a.rName.strip() == '':
            a.rName = 'XXX'
        if a.aName.strip() == '':
            a.aName = 'XX'

def writeFiles(fileBase,AL,FF,exceptResidues=[]):
    fixAtomNames(AL)
    ST = structure.Structure(AL,FF,exceptResidues=exceptResidues)
    ut = Cont() # keep track of the USED TYPES so that we can make a par file afterwards
    ut.atomtypes = ST.usedAtomLabels
    writePsf(fileBase+'.psf',AL,FF,ST,ut)
    writePrm(fileBase+'.prm',FF,ut)
    pdb.writeFile(fileBase+'.pdb',AL,withConect=False)

def writePsf(fileName,AL,FF, ST, ut):
    """write psf file file, force field is required to check if the structure can be parsed""" 
    
    with open(fileName,'w') as f:
        f.write('PSF\n')
        f.write('\n')
        f.write('1 !NTITLE\n')
        f.write("REMARKS generated by vaclav's psf.py on %s\n" % time.strftime('%x  %X'))
        f.write('\n')
        f.write('%d !NATOM\n' % len(AL))
        f.write('! atom ID, segment name, residue ID, residue name, atom name, atom type, charge, mass, and an unused 0\n')
        format = '%d %s %d %s %s %s %f %f 0\n'
        for a in AL:
            mass = FF.atomTypes[a.ffType].mass
            line = format % (a.aNo, a.chain, a.rNo, a.rName, a.aName, a.ffType, a.charge, mass)
            f.write(line)
            
        f.write('\n')
        f.write('%d !NBOND: bonds\n' % len(ST.bonds))
        i = 0
        for b in ST.bonds:
            f.write('%d %d ' % (b.i.aNo,b.j.aNo))
            i = i + 1
            if i % 4 == 0:
                f.write('\n')   
            # storing the bond type
            cd = ut.bonds 
            label1 = b.i.ffType+b.j.ffType
            label2 = b.j.ffType+b.i.ffType
            if (not label1 in cd) and (not label2 in cd):
                cd[label1] = b.type               
            
        f.write('\n')
        f.write('\n')
        f.write('%d !NTHETA: angles\n' % len(ST.angles))
        i = 0
        for b in ST.angles:
            f.write('%d %d %d ' % (b.i.aNo,b.j.aNo,b.k.aNo))
            i = i + 1
            if i % 3 == 0:
                f.write('\n')
            # storing the angle type
            cd = ut.angles 
            label1 = b.i.ffType+b.j.ffType+b.k.ffType
            label2 = b.k.ffType+b.j.ffType+b.i.ffType
            if (not label1 in cd) and (not label2 in cd):
                cd[label1] = b.type               
        
        f.write('\n')
        f.write('\n')
        # get unique dihedrals
        dih = []
        if len(ST.torsions) > 0:
            dih.append(ST.torsions[0])
        for i in range(1,len(ST.torsions)):
            if not ST.torsions[i-1].hasSameAtoms(ST.torsions[i]):
                dih.append(ST.torsions[i])
                
        # writing out torsions
        f.write('%d !NPHI: dihedrals\n' % len(dih))
        i = 0
        for b in dih:
            f.write('%d %d %d %d ' % (b.i.aNo,b.j.aNo,b.k.aNo, b.l.aNo))
            i = i + 1
            if i % 2 == 0:
                f.write('\n')
                
        # storing the torsion type
        for b in ST.torsions:
            cd = ut.torsions
            n = '%2d' %b.type.N
            label1 = b.i.ffType+b.j.ffType+b.k.ffType+b.l.ffType+n
            label2 = b.l.ffType+b.k.ffType+b.j.ffType+b.i.ffType+n
            if (not label1 in cd) and (not label2 in cd):
                cd[label1] = b.type  

        # now doing inversions
        f.write('\n')
        f.write('\n')
        f.write('%d !NIMPHI: impropers\n' % len(ST.inversions))
        i = 0
        for b in ST.inversions:
            # storing the inversion type
            cd = ut.inversions
            labeli = b.i.ffType+b.j.ffType+b.k.ffType+b.l.ffType
            if not labeli in cd:
                cd[labeli]=b.type
            # writing the 
            f.write('%d %d %d %d ' % (b.i.aNo,b.j.aNo,b.k.aNo,b.l.aNo))
            i = i + 1
            if i % 2 == 0:
                f.write('\n')
        
        f.write('\n')
        f.write('\n')
        f.write('0 !NDON: donors\n')
        f.write('\n')
        f.write('\n')
        f.write('0 !NACC: acceptors\n')
        f.write('\n')
        f.write('\n')
        f.write('0 !NNB: is number of explicit nonbond exclusions\n') 
        f.write('\n')
        i = 0
        for a in AL:
            f.write('0 ')
            i = i + 1
            if i % 8 == 0:
                f.write('\n')
                i = 0
        f.write('\n')
        f.write('1 0 !NGRP: is number of groups (not important for your NAMD simulation I think) \n')
        f.write('0 0 0\n')


def writePrm(fileName,FF,ut):
    """write this amber force field in charmm par file for NAMD """ 
    # n = number of characters in atom label
    if FF.forcefield == 'amber':
        n = 2
    elif FF.forcefield == 'dreiding':
        n = 5 
    else:
        print >> sys.stderr,  'ERROR: only can output psf for amber and for dreiding (without hbond), but have %s' % FF.forcefield
        sys.exit(1)
          
    with open(fileName,'w') as f:
        f.write('! par file readable by NAMD, contains amber data in charmm format\n')
        f.write('BONDS\n')
        f.write('!V(bond) = Kb(b - b0)**2\n')
        f.write('!Kb: kcal/mole/A**2\n')
        f.write('!b0: A\n')
        f.write('!atom type Kb b0\n')
        f.write('!\n')
        cd = ut.bonds
        for b in cd:
            i = b[0:n]
            j = b[n:(2*n)]
            t = cd[b]
            f.write('%s %s %f %f \n' % (i,j,t.K/2.0,t.R))
        f.write('\n')
        f.write('ANGLES\n')
        f.write('!V(angle) = Ktheta(Theta - Theta0)**2\n')
        f.write('!V(Urey-Bradley) = Kub(S - S0)**2\n')
        f.write('!Ktheta: kcal/mole/rad**2\n')
        f.write('!Theta0: degrees\n')
        f.write('!Kub: kcal/mole/A**2 (Urey-Bradley)\n')
        f.write('!S0: A\n')
        f.write('!atom types Ktheta Theta0 Kub S0  (last two are optional)\n')
        f.write('!\n')
        #if len(FF.angleType1) > 0:
        #    print >> sys.stderr, 'FF.angleType1 has entries, cannot convert this to the charmm format'
        #    sys.exit(1)
        cd = ut.angles
        for b in cd:
            i = b[0:n]
            j = b[n:(2*n)]
            k = b[(2*n):(3*n)]
            t = cd[b]
            f.write('%s %s %s %f %f \n' % (i,j,k,t.K/2.0,t.theta0))
        f.write('\n')
        f.write('DIHEDRALS\n')
        f.write('!V(dihedral) = Kchi(1 + cos(n(chi) - delta))\n')
        f.write('!Kchi: kcal/mole\n')
        f.write('!n: multiplicity\n')
        f.write('!delta: degrees\n')
        f.write('!atom types Kchi n delta\n')
        f.write('!\n')
        cds = []
        for b in ut.torsions: cds.append(b)
        cds.sort()
        cd = ut.torsions
        for b in cds:
            i = b[0:n]
            j = b[n:(2*n)]
            k = b[(2*n):(3*n)]
            l = b[(3*n):(4*n)]
            t = cd[b]        
            K = t.K / 2.0
            N = t.N
            D = prmtop.dToPhase(t.D)
            f.write('%s %s %s %s %f %d %f \n' % (i,j,k,l,K,N,D))

        f.write('\n')
        f.write('IMPROPER\n')
        f.write('!V(improper) = Kpsi(psi - psi0)**2\n')
        f.write('!Kpsi: kcal/mole/rad**2\n')
        f.write('!psi0: degrees\n')
        f.write('!note that the second column of numbers (0) is ignored\n')
        f.write('!atom types           Kpsi                   psi0\n')
        f.write('!\n')
        # the angle is probably defined the same say in charmm and ambrer
        # charmm: angle between the plane containing the first three atoms and last three atoms
        #            K*(psi-psi0)^2
        # amber: as dihedral where the third atom is the central atom
        #                   
        cd = ut.inversions
        for b in cd:
            i = b[0:n]
            j = b[n:(2*n)]
            k = b[(2*n):(3*n)]
            l = b[(3*n):(4*n)]            
            t = cd[b]
            coeff = t.charmm()
            f.write('%s %s %s %s %f 0 0.0 \n' % (i,j,k,l,coeff))
        
        f.write('\n')
        f.write('NONBONDED nbxmod 5 atom cdiel shift vatom vdistance vswitch -\n')
        f.write('! cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n') 
        f.write('! the above line used to be actual parameter skipped line\n')
        f.write('!\n')
        f.write('!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n')
        f.write('!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)\n')
        f.write('!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j\n')
        f.write('!atom ignored epsilon Rmin/2 ignored eps,1-4 Rmin/2,1-4\n')
        f.write('!\n')
        cd = ut.atomtypes        
        for b in cd:
            i = b
            t = FF.atomTypes[b]
            R = t.vdwR / 2.0
            E = - t.vdwD
            E14 = FF.special14lj*E
            f.write('%s 0 %f %f 0 %f %f \n' % (i,E,R,E14,R))
        f.write('\n')
        f.write('END\n')

def fixMartiniPsf(psffilein,psffileout,sequenceString):
    '''inputs the correct backbone beads into the martini coarse grain psf file

    the protein residue names must be VMD/charmm compatible, but this is the case
    if the psf file is generated in vmd 

    there should be no NTER or CTER patched -- not sure where they can come from in the
    coarse grain context

    sequence string encoding (not case sensitive):
        B bend
        C coil
        F free
        H helix
        T turn        
        E extended/beta
    '''
    # recognized protein residues
    protein = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HSD','HSE','HSP','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    # the psffile
    lines = open(psffilein).readlines()
    linesout = []
    # the header
    i = 0
    if lines[i].startswith('PSF'):
        linesout.append(lines[i])
    else:
        sys.exit('Error: fixMartiniPsf: first line is no PSF')
    i += 1
    if lines[i].strip() == '':
        linesout.append(lines[i])
    else:
        sys.exit('Error: fixMartiniPsf: second line is not blank: %s' % lines[i].strip())
    i += 1
    s = lines[i].split()
    nlines = int(s[0])
    line = '%8d !NTITLE\n' % (nlines+2)
    linesout.append(line)
    for j in range(nlines):
        i += 1
        linesout.append(lines[i])

    linesout.append(' REMARKS vcvicek added the secondary structure of the backbone beads\n')    
    linesout.append(' REMARKS vcvicek pyvasek.psf.fixMartiniPsf\n')    
    
    i += 1
    if lines[i].strip() == '':
        linesout.append(lines[i])
    else:
        sys.exit('Error: fixMartiniPsf: line is not blank: %s' % lines[i].strip())
    i += 1
    linesout.append(lines[i])
    s = lines[i].split()
    nlines = int(s[0])

    # the individual toms     
    for j in range(nlines):
        i += 1
        line = lines[i]
        rNo = int(line[14:19])
        rName = line[19:22]
        bType = line[24:27]
        bFF =   line[29:33]
        if (rName not in protein) or (bType != 'BAS'):
            linesout.append(line)
        else:
            if rNo > len(sequenceString):
                sys.exit('Error: fixMartiniPsf: sequenceString does not contain res %d' % (rNo))
            seq = sequenceString[rNo-1].upper()
            if seq in ['B','C','F']:
                if rName == 'PRO':
                    newbFF = 'NA'+seq
                elif rName == 'ALA':
                    newbFF = 'P4'+seq
                else:
                    newbFF = 'P5'+seq
            elif seq in ['H']:
                if (rNo == 1) or (sequenceString[rNo-2].upper() != 'H'):
                    # N end of helix
                    if rName == 'PRO':
                        newbFF = 'N0HP'
                    elif rName == 'ALA':
                        newbFF = 'N0H'
                    else:
                        newbFF = 'NDH'
                elif (rNo == len(sequenceString)) or (sequenceString[rNo].upper() != 'H'):
                    # C end of helix 
                    if rName == 'PRO':
                        newbFF = 'NAH'
                    elif rName == 'ALA':
                        newbFF = 'N0H'
                    else:
                        newbFF = 'NAH'
                else:
                    # just helix
                    if rName == 'PRO':
                        newbFF = 'C5HP'
                    elif rName == 'ALA':
                        newbFF = 'C5H'
                    else:
                        newbFF = 'N0H'
            elif seq in ['T','E']:
                if rName in ['PRO','ALA']:
                    newbFF = 'N0'+seq
                else:
                    newbFF = 'NDA'+seq
            else:
                sys.exit('Error: fixMartiniPsf: unknown character in sequenceString: %s' % seq)
            newline = line[:29]+'%-4s'%newbFF+line[33:]
            # finished lines
            linesout.append(newline)

    # finish the file
    for j in range(i+1,len(lines)):
        linesout.append(lines[j])

    with open(psffileout,'w') as f:
        f.writelines(linesout)


