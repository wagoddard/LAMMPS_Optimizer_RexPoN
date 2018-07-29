"""atom typing for Dreiding 3

note that S_4, P_4, O_3M, N_3P, etc, will never be produced by this typing algorithm (they also wouldn't be produced by Lingraf typing)
"""

import sys, os, tempfile, shutil, subprocess
#from packages import pybel # this is needed for sensing bond orders only
import pdb, atoms, bgf, defaults, mol2

def getElement(a):
    ''' this also creates field 'two letter' with the element name in the atom instance
    parameter a is one atoms.Atom instance
    Note:
      calcium has to be Ca, because CA is carbon
      sodium has to be Na, because NA is nitrogen
      mercury has to be Hg, because HG is hydrogen
      other two letter atom names can be capitalized
    '''
    aname = a.aName.strip()
    if aname == '':
        print >> sys.stderr, 'Error: This atom has no label, so I cannot determine the element %s ' % str(a) 
        sys.exit(1)        
    if aname[0:1].islower():
        newname = aname[0:1].upper() + aname[1:]
        print 'Warning: chaning atom name of atom %s from %s to %s' % (a, aname,newname)
        aname = newname
        a.aName = aname
    if aname[0:2] in ['CL',
                        'BR',
                        'FE',
                        'ZN',
                        'TI',
                        'TC',
                        'RU',
                        'MG',
                        'MN',
                        'GA',
                        'GE',
                        'AS',
                        'SE',
                        'IN',
                        'SN',
                        'SB',
                        'TE',
                        'AL',
                        'SI']:
        newname = aname[0:1] + aname[1:2].lower() + aname[2:]
        print 'Warning: chaning atom name of atom %s from %s to %s' % (a, aname,newname)
        aname = newname
        a.aName = aname
        a.twoletter = aname[0:2] 
    if aname[1:2].islower():
        a.twoletter = aname[0:2]
    else:
        a.twoletter = aname[0:1]+'_'

def getAtomHybridization(a):
    '''computation of hybridization is based on the original lingraf code gethyb.f
    note that here we allow bond order 1.3 (amide) and 1.5 (aromatic), so the code has to be a bit updated
    the hybridization is stored in the atom instance
    returns number>0 for C,O,S,N,P otherwise returns number 0
    '''
    numberOfBonds = len(a.bonds)                                # NBOND
    if numberOfBonds == 0:
        print 'Warning: atom %s has 0 bonds' % a    
        maximalOrder = 0                                        # MAXORD
        numberWithMaximalOrder = 0                              # NUMMAX        
    else:
        maximalOrder = max(a.order)                             # MAXORD
        numberWithMaximalOrder = a.order.count( maximalOrder )  # NUMMAX

    if numberOfBonds > 4:
        print >> sys.stderr, 'Error: Atom %s has more than 4 bonds: %d ' % (a,numberOfBonds) 
        sys.exit(1)  

    # default hybridization for elements not present here    
    hybridization = 0
    element = a.twoletter
    if element in ['C_']: # gethyb.f had the same rule for Si and Al
        # ketene
        if maximalOrder == 2 and numberWithMaximalOrder == 2 and numberOfBonds == 2:
            hybridization = 1
        # check just maximalOrder  
        elif maximalOrder == 3:
            hybridization = 1
        elif maximalOrder == 1:
            hybridization = 3
        else:
            hybridization = 2
    elif element in ['O_','S_']:
        if maximalOrder >= 3:
            hybridization = 1
        elif maximalOrder > 1:
            hybridization = 2
        else:
            hybridization = 3 
    elif element == 'N_':
        # the order of the rules here matters
        neighborOrders = [max(b.order) for b in a.bonds]
        maximalOrderOfNeighbors = max(neighborOrders)               # MAXNEI
        # sp3 nitrogens
        if numberOfBonds >= 4:
            hybridization = 3
        # sp nitrogens
        elif maximalOrder == 3:
            hybridization = 1
        elif maximalOrder == 2 and numberWithMaximalOrder == 2 and numberOfBonds == 3: # new case not in Lingraf, when -NO2 has two double bonds
            hybridization = 2
        elif maximalOrder == 2 and numberWithMaximalOrder == 2:
            hybridization = 1        
        # sp2 nitrogens
        # I am quite not sure about this rule
        elif maximalOrder == 1 and maximalOrderOfNeighbors > 1:
            hybridization = 2
        # check just maximalOrder
        elif maximalOrder == 1:
            hybridization = 3
        else:
            hybridization = 2
    elif element == 'P_':
        hybridization = 3            
    a.hybridization = hybridization 

def type(al): 
    '''gettign the atom type and number of lone pairs'''
    for a in al:
        getElement(a)
        getAtomHybridization(a) # returns number>0 for C,O,S,N,P otherwise returns number 0 
    for a in al:
        # the simplest case
        nbonds = len(a.bonds)
        if a.twoletter in ['F_',
                            'Cl',
                            'Br',
                            'I_',
                            'Na',
                            'Ca',
                            'Fe',
                            'Zn',
                            'Ti',
                            'Tc',
                            'Ru',
                            'Hg',
                            'Mg',
                            'Mn',
                            'K_']:
            a.ffType = a.twoletter
            a.lpair = 0
        # ffType has 3 attached
        elif a.twoletter in ['Ga',
                            'Ge',
                            'As',
                            'Se',
                            'In',
                            'Sn',
                            'Sb',
                            'Te',
                            'Al',
                            'Si']:
            a.ffType = a.twoletter+'3'
            a.lpair = 0
        # hydrogen 
        elif a.twoletter == 'H_':
            if nbonds == 2:
                a.ffType = 'H___b'     # hydrogen in diborane
                a.lpair = 0
            elif nbonds == 1:
                b = a.bonds[0]
                if b.twoletter in ['O_','N_','S_','F_','Cl','Br','I_']:
                    a.ffType = 'H___A' # polar hydrogen
                    a.lpair = 0
                else:
                    a.ffType = 'H_'   # normal hydrogen
                    a.lpair = 0                    
            else:
                print >> sys.stderr, 'Error: This hydrogen is not attached to anything, probably an error: %s' % (str(a)) 
                sys.exit(1)     
        # boron
        elif a.twoletter == 'B_':
            if nbonds == 4:
                a.ffType = 'B_3'
                a.lpair = 0
            elif nbonds in [2,3]:
                a.ffType = 'B_2'
                a.lpair = 3-nbonds
            else:
                print >> sys.stderr, 'Error: I only know B with 2, 3 or 4 bonds, but here there is %d bonds, %s' % (nbonds,str(a)) 
                sys.exit(1)     
        # phosphor
        elif a.twoletter == 'P_':
            if nbonds == 4:
                a.ffType = 'P_4'
                a.lpair = 0
            elif nbonds == 3:
                a.ffType = 'P_3'
                a.lpair = 1
            else:
                print >> sys.stderr, 'Error: I only know P with 4 or 3 bonds, but here there is %d bonds, %s' % (nbonds,str(a)) 
                sys.exit(1)                
        # sulphur
        elif a.twoletter == 'S_':
            if nbonds == 4:
                a.ffType = 'S_4'
                a.lpair = 0
            elif nbonds == 1 and a.hybridization == 2:
                a.ffType = 'S_2'
                a.lpair = 2
            elif nbonds == 2 and a.hybridization == 2:
                a.ffType = 'S_2'
                a.lpair = 1
            elif nbonds == 2 and a.hybridization == 3:
                a.ffType = 'S_3'
                a.lpair = 2
            elif nbonds == 3 and a.hybridization == 2 and sorted(a.order) == [1,1,2]: # 8/1/2013 added for sulfoxides on adam's request 
                a.ffType = 'S_4'
                a.lpair = 1
            else:
                print >> sys.stderr, 'Error: I do not know S with %d bonds and hybridization %d. %s' % (nbonds,a.hybridization,str(a)) 
                sys.exit(1)                
        # only C, O, N are left 
        # oxygen
        elif a.twoletter == 'O_':
            if a.hybridization == 1:
                if nbonds == 1:
                    a.ffType = 'O_1'
                    a.lpair = 1
                else:
                    print >> sys.stderr, 'Error: sp1 oxygen can only have one bond' 
                    sys.exit(1)                                                
            elif a.hybridization == 2:
                # count number of sp2 neighbors to decide if this oxygen is resonant
                sp2neighbors = 0
                for b in a.bonds:
                    if b.hybridization == 2:
                         sp2neighbors += 1
                if sp2neighbors >= 2: # O_R
                    if nbonds in [2,3]: 
                        a.ffType = 'O_R'
                        a.lpair = 3-nbonds
                    else:
                        print >> sys.stderr, 'Error: this resonant oxygen %s has too many bonds: %d ' % (str(a),nbonds) 
                        sys.exit(1)
                else:
                    if nbonds not in [1,2]:
                        print >> sys.stderr, 'Error: this sp2 oxygen %s should have 1 or 2 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                    a.ffType = 'O_2'
                    a.lpair = 3 - nbonds
            else: # oxygen hybridization 3
                # find out if this oxygen has any sp2 neighbor
                sp2neighbors = 0
                for b in a.bonds:
                    if b.hybridization == 2:
                         sp2neighbors += 1
                         break
                if sp2neighbors >= 1 and nbonds == 2: # pyran, or phenol
                    a.ffType = 'O_R'
                    a.lpair = 1
                elif nbonds == 1: # test for COOH
                    hasO_2neighbor = False
                    for b in a.bonds[0].bonds:
                        if b.twoletter == 'O_' and b.hybridization == 2:
                            hasO_2neighbor = True
                            break
                    if hasO_2neighbor:  
                        a.ffType = 'O_2'  # OH in COOH
                        a.lpair = 2
                    else:
                        a.ffType = 'O_3'
                        a.lpair = 3                    
                else:
                    if nbonds not in [1,2,3]:
                        print >> sys.stderr, 'Error: this sp3 oxygen %s should have 1,2 or 3 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                    a.ffType = 'O_3'
                    a.lpair = 4 - nbonds
        # carbon and nitrogen
        elif a.twoletter in ['C_','N_'] :
            if a.hybridization == 3:
                if a.twoletter == 'C_':
                    if nbonds == 4:
                        a.ffType = 'C_3'
                        a.lpair = 0
                    else:                                            
                        print >> sys.stderr, 'Error: sp3 carbon %s has to have 4 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                elif a.twoletter == 'N_':
                    if nbonds in [2,3,4]:
                        a.ffType = 'N_3'
                        a.lpair = 4-nbonds
                    else:
                        print >> sys.stderr, 'Error: sp3 nitrogen %s has to have 2,3 or 4 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                else:
                    print >> sys.stderr, 'Error: this part of code should never be reached 1234...' 
                    sys.exit(1)
            elif a.hybridization == 1:
                if a.twoletter == 'C_':
                    if nbonds == 2:
                        a.ffType = 'C_1'
                        a.lpair = 0
                    else:                                            
                        print >> sys.stderr, 'Error: sp1 carbon %s has to have 2 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                elif a.twoletter == 'N_':
                    if nbonds in [1,2]:
                        a.ffType = 'N_1'
                        a.lpair = 2-nbonds
                    else:
                        print >> sys.stderr, 'Error: sp1 nitrogen %s has to have 1 or 2 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                else:
                    print >> sys.stderr, 'Error: this part of code should never be reached 1235...' 
                    sys.exit(1)
            elif a.hybridization == 2:
                # sp2 can be _2 or _R (resonant)
                makeResonant = False
                # count sp2 neighbors and Osp3 oxygen neighbors
                sp2neighbors = 0
                Osp3neighbors = 0
                for b in a.bonds:
                    if b.hybridization == 2:
                         sp2neighbors += 1
                         sp2atom = b
                    if b.twoletter == 'O_' and b.hybridization == 3:
                         Osp3neighbors += 1
                if sp2neighbors>=2 or (sp2neighbors==1 and Osp3neighbors>=1):
                    # we found one patter in which we have to make sp2 resonant
                    makeResonant = True
                elif sp2neighbors==1:
                    # go through the neighbors and check their neighbors
                    for b in sp2atom.bonds:
                        if b!=a:
                            if b.hybridization==2 or (b.hybridization==3 and b.twoletter=='O_'):
                                makeResonant = True
                                break
                if makeResonant:
                    if a.twoletter == 'C_':
                        if nbonds == 3:
                            a.ffType = 'C_R'
                            a.lpair = 0
                        else:                                            
                            print >> sys.stderr, 'Error: sp2 carbon %s has to have 3 bonds, but has %d ' % (str(a),nbonds) 
                            sys.exit(1)
                    elif a.twoletter == 'N_':
                        if nbonds in [2,3]:
                            a.ffType = 'N_R'
                            a.lpair = 3-nbonds
                        else:
                            print >> sys.stderr, 'Error: sp2 (R_) nitrogen %s has to have 2 or 3 bonds, but has %d ' % (str(a),nbonds) 
                            sys.exit(1)
                    else:
                        print >> sys.stderr, 'Error: this part of code should never be reached 1237...' 
                        sys.exit(1)
                else:
                    if a.twoletter == 'C_':
                        if nbonds == 3:
                            a.ffType = 'C_2'
                            a.lpair = 0
                        else:                                            
                            print >> sys.stderr, 'Error: sp2 carbon %s has to have 3 bonds, but has %d ' % (str(a),nbonds) 
                            sys.exit(1)
                    elif a.twoletter == 'N_':
                        if nbonds in [1,2,3]:
                            a.ffType = 'N_2'
                            a.lpair = 3-nbonds
                        else:
                            print >> sys.stderr, 'Error: sp2 nitrogen %s has to have 1, 2 or 3 bonds, but has %d ' % (str(a),nbonds) 
                            sys.exit(1)
                    else:
                        print >> sys.stderr, 'Error: this part of code should never be reached 1238...' 
                        sys.exit(1)
            else:
                print >> sys.stderr, 'Error: this part of code should never be reached 1236...' 
                sys.exit(1)
            # end of carbon and nitrogen
        else:
            print >> sys.stderr, 'Error: There is no algorithm to type atom: "%s"' % a.twoletter 
            sys.exit(1)            

""" this is a less robust older version
def type(al): 
    '''gettign the atom type and number of lone pairs'''
    for a in al:
        getElement(a)
        getAtomHybridization(a)
    for a in al:
        # the simplest case
        if a.twoletter in ['F_',
                            'Cl',
                            'Br',
                            'I_',
                            'Na',
                            'Ca',
                            'Fe',
                            'Zn',
                            'Ti',
                            'Tc',
                            'Ru',
                            'Hg',
                            'Mg',
                            'Mn',
                            'K_']:
            a.ffType = a.twoletter
            a.lpair = 0
        # ffType has 3 attached
        elif a.twoletter in ['Ga',
                            'Ge',
                            'As',
                            'Se',
                            'In',
                            'Sn',
                            'Sb',
                            'Te',
                            'Al',
                            'Si']:
            a.ffType = a.twoletter+'3'
            a.lpair = 0
        # hydrogen 
        elif a.twoletter == 'H_':
            if len(a.bonds) == 2:
                a.ffType = 'H___b'
                a.lpair = 0
            elif len(a.bonds) == 1:
                b = a.bonds[0]
                if b.twoletter in ['O_','N_','S_','F_','Cl','Br','I_']:
                    a.ffType = 'H___A'
                    a.lpair = 0
                else:
                    a.ffType = 'H_'
                    a.lpair = 0                    
            else:
                a.ffType = 'H_'
                a.lpair = 0
        # boron
        elif a.twoletter == 'B_':
            if len(a.bonds) == 4:
                a.ffType = 'B_3'
                a.lpair = 0
            else:
                a.ffType = 'B_2'
                a.lpair = 0
        # phosphor
        elif a.twoletter == 'P_':
            a.ffType = 'P_3'
            a.lpair = 4-len(a.bonds) # this is just a guess since we are not analyzing all phosphors
        # sulphur
        elif a.twoletter == 'S_':
            a.ffType = 'S_3'
            a.lpair = 4-len(a.bonds) # this is just a guess since we are not analyzing all sulphurs
        elif a.twoletter in ['C_','O_','N_'] : # only C, O, N are left
            # oxygen
            if a.twoletter == 'O_':
                if a.hybridization == 1:
                    a.ffType = 'O_1'
                    a.lpair = 2 - len(a.bonds)
                elif a.hybridization == 2:
                    # count number of sp2 neighbors
                    sp2neighbors = 0
                    for b in a.bonds:
                        if b.hybridization == 2:
                             sp2neighbors += 1
                    if sp2neighbors >= 2:
                        a.ffType = 'O_R'
                        a.lpair = 3 - len(a.bonds)
                    else:
                        a.ffType = 'O_2'
                        a.lpair = 3 - len(a.bonds)
                else: # oxygen hybridization 3
                    # find out if this oxygen has some sp2 neighbor
                    sp2neighbors = 0
                    nbonds = len(a.bonds)
                    for b in a.bonds:
                        if b.hybridization == 2:
                             sp2neighbors += 1
                    if sp2neighbors >= 1 and nbonds >= 2:
                        a.ffType = 'O_R'
                        a.lpair = 3 - len(a.bonds)
                    elif nbonds == 1:
                        hasO_2neighbor = False
                        for b in a.bonds[0].bonds:
                            if b.twoletter == 'O_' and b.hybridization == 2:
                                hasO_2neighbor = True
                                break
                        if hasO_2neighbor:
                            a.ffType = 'O_2'
                            a.lpair = 3 - nbonds
                        else:
                            a.ffType = 'O_3'
                            a.lpair = 4 - nbonds                    
                    else:
                        a.ffType = 'O_3'
                        a.lpair = 4 - nbonds
            # carbon and nitrogen
            else:
                if a.hybridization == 3:
                    a.ffType = a.twoletter + '3'
                    a.lpair = 4 - len(a.bonds)
                elif a.hybridization == 1:
                    a.ffType = a.twoletter + '1'
                    a.lpair = 2 - len(a.bonds)
                else: # hybridization 2
                    # count sp2 neighbors and Osp3 oxygen neighbors
                    sp2neighbors = 0
                    Osp3neighbors = 0
                    nbonds = len(a.bonds)
                    for b in a.bonds:
                        if b.hybridization == 2:
                             sp2neighbors += 1
                             sp2atom = b
                        if b.twoletter == 'O_' and b.hybridization == 3:
                             Osp3neighbors += 1
                    if sp2neighbors>=2 or (sp2neighbors==1 and Osp3neighbors>=1):
                        a.ffType = a.twoletter + 'R'
                        a.lpair = 3-nbonds 
                    elif sp2neighbors==1:
                        # go through the neighbors and check their neighbors
                        makeResonant = False
                        for b in sp2atom.bonds:
                            if b!=a:
                                if b.hybridization==2 or (b.hybridization==3 and b.twoletter=='O_'):
                                    makeResonant = True
                                    break
                        if makeResonant:
                            a.ffType = a.twoletter + 'R'
                            a.lpair = 3-nbonds
                        else:
                            a.ffType = a.twoletter + '2'
                            a.lpair = 3-nbonds
                    else:
                        a.ffType = a.twoletter + '2'
                        a.lpair = 3-nbonds
        else:
            print >> sys.stderr, 'Error: There is no algorithm to type atom: "%d"' % a.twoletter 
            sys.exit(1)            
end of older version
"""

def typeBond(al, singleDouble = False):
    typeBondAntechamber(al, singleDouble = singleDouble)

def typeBondBabel(al):
    '''this assignes bond order to atom list using openbabel
    first pdb is created
    then it is read by python openbabel (pybel)
    then it is written into bgf
    then only the bond orders from the new bgf are used
    '''
    from packages import pybel

    mylines = pdb.writeLines(al)
    mylines = ''.join(mylines)
    pmol = pybel.readstring('pdb',mylines)
    mylines = pmol.write(format='bgf').split('\n')
    al2 = bgf.readLines(mylines)
    '''
    # alternative version
    pdb.write('tempMol.pdb',al)
    babel ---errorlevel 2 -ipdb tempMol.pdb -obgf tempMol.bgf
    al2 = bgf.read('tempMol.bgf')
    '''
    # 
    atoms.sortBondsLists(al)
    atoms.sortBondsLists(al2)
    for i in range(len(al)):
        a1 = al[i]
        a2 = al2[i]
        a1.order = a2.order[:]

def typeBondShell(al):
    '''runs on shell
    this assignes bond order to atom list using openbabel
    first pdb is created
    then it is read by python openbabel (pybel)
    then it is written into bgf
    then only the bond orders from the new bgf are used
    '''
    import os
    from pyvasek import mol2

    pdb.write('tempMol.pdb',al)
    os.system('babel ---errorlevel 2 -ipdb tempMol.pdb -omol2 tempMol.mol2')
    al2 = mol2.read('tempMol.mol2')
    # 
    atoms.sortBondsLists(al)
    atoms.sortBondsLists(al2)
    for i in range(len(al)):
        a1 = al[i]
        a2 = al2[i]
        a1.order = a2.order[:]

def typeBondAntechamber(al, singleDouble = False):
    """this uses bondtype program from ambertools package
    bonds have to be present in the input structure, but the bond orders may be incorrect
    the bondtype program fixes the bond orders"""
    if singleDouble:
        bondDict = {1:1, # single bond
            2:2, # double bond 
            3:3, # tripple bond 
            7:1, # aromatic single 
            8:2, # aromatic double
            9:1.5, # delocalized bond (C-O bonds in carboxyl anions and nitro groups)
            6:1} # conjugated
    else:
        bondDict = {1:1, # single bond
            2:2, # double bond 
            3:3, # tripple bond 
            7:1.5, # aromatic single 
            8:1.5, # aromatic double
            9:1.5, # delocalized bond (C-O bonds in carboxyl anions and nitro groups)
            6:1} # conjugated

    amberhome = os.getenv("AMBERHOME")
    if amberhome == None:
        sys.exit('Error: enrironment variable $AMBERHOME is not defined, but the bond typing uses bondtype program from ambertools')
    atoms.renumberAtomsFrom(al)

    # will run the bondtype job in a temporary directory
    tmpdir = tempfile.mkdtemp(prefix='tmpBondtype',dir='.')
    os.chdir(tmpdir)
    # 
    mol2.write('tempMol.mol2',al)
    # os.system('%s/bin/bondtype -i tempMol.mol2 -o tempOut.ac -f mol2 -j full' % amberhome)
    command ='%s/bin/bondtype -i tempMol.mol2 -o tempOut.ac -f mol2 -j full' % amberhome
    subprocess.call(command.split())
    lines = open('tempOut.ac').readlines()
    for line in lines:
        if line.startswith('BOND'):
            s = line.split()
            atom1 = al[int(s[2])-1]
            atom2 = al[int(s[3])-1]
            orderKey = int(s[4])
            atom1.deleteBond(atom2)
            if orderKey not in bondDict:
                sys.exit('Error in typeBondAntechamber: bondtype returned bond type %d, which is unknown '%orderKey)
            atom1.makeBond(atom2, bondOrder = bondDict[orderKey])
    os.chdir('..')
    shutil.rmtree(tmpdir)

def readCnv(cnvfile=''):
    '''returns dictionary  [RESatom.strip()] -> list of : atom type, number of bonds, number of lone pair, charge (but sometiems charge is not present)'''
    if cnvfile == '':
        cnvfile = defaults.dreiding3cnv
    lines = open(cnvfile).readlines()
    name2charge = {}
    for line in lines:
        s = line.split()
        if line[0:1] not in ['*','C','!'] and len(s)>0: # comment characters
            if len(s) < 5:
                print >> sys.stderr, 'Error: this line of cnv file does not have the right number of items: %s' % line
                sys.exit(1)                  
            rName = s[0]
            aName = s[1]
            ffType = '%-5s'% s[2]
            nbonds = int(s[3])
            lpair = int(s[4])
            if len(s) >= 6:
                charge = float(s[5])
                value = [ffType,nbonds,lpair,charge]
            else:
                value = [ffType,nbonds,lpair]
            if rName+aName in name2charge:
                usedvalue = name2charge[rName+aName]
                if len(value)!=len(usedvalue):
                    print >> sys.stderr, 'Error: res %s and atoms %s are twice in the cnv file ' % (rName, aName)
                    sys.exit(1)                
                for i in range(len(usedvalue)):
                    if usedvalue[i]!=value[i]:
                        print >> sys.stderr, 'Error: res %s and atoms %s are twice in the cnv file ' % (rName, aName)
                        sys.exit(1)                
            else:
                name2charge[rName+aName] = value                 
    return name2charge
            
def assignCharge(al,assignFfType=False,cnvfile=''):
    '''assign cahrge and optionally ffType from the cnv file for all atoms in the atom list al'''
    name2charge = readCnv(cnvfile=cnvfile)
    for a in al:
        aName = a.aName.strip()
        rName = a.rName.strip()
        key = rName + aName
        if key in name2charge:
            value = name2charge[key]
            if assignFfType:
                a.ffType = value[0]
            if len(value) >= 4:
                a.charge = value[3]
            else:
                print >> sys.stderr, 'Warning: Atom %s is in the cnv file without charge, so its charge is left %f' % (str(a),a.charge)
            if len(a.bonds) != value[1] or a.lpair != value[2]:
                print >> sys.stderr, 'Warning: atom %s had %d bonds and %d lone pairs, but cnv has %d bonds and %d lone pairs' % (str(a),len(a.bonds),a.lpair,value[1],value[2])                   
        elif '***'+aName in name2charge:
            value = name2charge['***'+aName]
            if assignFfType:
                a.ffType = value[0]
            if len(value) >= 4:
                a.charge = value[3]
            else:
                print >> sys.stderr, 'Warning: Atom %s is in the cnv file without charge, so its charge is left %f' % (str(a),a.charge)
            if len(a.bonds) != value[1] or a.lpair != value[2]:
                print >> sys.stderr, 'Warning: atom %s had %d bonds and %d lone pairs, but cnv has %d bonds and %d lone pairs' % (str(a),len(a.bonds),a.lpair,value[1],value[2])                   
        else:
            print >> sys.stderr, 'Warning: Atom %s is not in the cnv file, so its charge is left %f' % (str(a),a.charge)               
