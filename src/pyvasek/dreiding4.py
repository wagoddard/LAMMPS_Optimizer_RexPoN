"""atom typing for Dreiding 4

The atom typing is done quite well for organic molecules, and there are even more checks than
in Dreiding 3 typing. However, the parameters were never optimized for Dreiding 4, so currently
the only use for Dreiding 4 typing is to check integrity of molecules (eg. when they are coming
from VLS with possible broken bonds or incorrect atoms).
"""

import sys
import dreiding3 

def valenceTest(atom,atomName,desiredValence, upperBound = 0):
    """utility for the type function"""
    realValence = sum(atom.order)
    if upperBound == 0:
        if realValence != desiredValence:
            print >> sys.stderr, 'Warning: Valence problem: This %s atom does not have sum of bond orders %d, but %f: %s' % (atomName, desiredValence, realValence, str(atom)) 
    else:
        if not (desiredValence <= realValence <= upperBound):
            print >> sys.stderr, 'Warning: Valence problem: This %s atom does not have sum of bond orders between %f and %f, but %f: %s' % (atomName, desiredValence, upperBound,realValence, str(atom)) 


def type(al): 
    '''gettign the atom type and number of lone pairs'''
    for a in al:
        dreiding3.getElement(a)
        dreiding3.getAtomHybridization(a) # returns number>0 for C,O,S,N,P otherwise returns number 0 
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
                valenceTest(a,'H',2)
            elif nbonds == 1:
                b = a.bonds[0]
                if b.twoletter == 'O_': # perhaps water hydrogen
                    if len(b.bonds) == 2:
                        if b.bonds[0].twoletter == 'H_' and b.bonds[1].twoletter == 'H_':
                            a.ffType = 'H___W' # water hydrogen
                            a.lpair = 0                           
                        else:
                            a.ffType = 'H___A' # polar hydrogen
                            a.lpair = 0
                    else:
                        a.ffType = 'H___A' # polar hydrogen
                        a.lpair = 0
                elif b.twoletter in ['N_','S_','F_','Cl','Br','I_']:
                    a.ffType = 'H___A' # polar hydrogen
                    a.lpair = 0
                else:
                    a.ffType = 'H_'   # normal hydrogen
                    a.lpair = 0  
                valenceTest(a,'H',1)
            else:
                print >> sys.stderr, 'Error: This hydrogen is not attached to anything, probably an error: %s' % (str(a)) 
                sys.exit(1)     
        # boron
        elif a.twoletter == 'B_':
            if nbonds == 4:
                a.ffType = 'B_3'
                a.lpair = 0
            elif nbonds == 3:
                a.ffType = 'B_2'
                a.lpair = 0
            # minimize the the number of atom types for now            
            #elif nbonds == 2:
            #    a.ffType = 'B_2'
            #    a.lpair = 1
            else:
                print >> sys.stderr, 'Error: I only know B with 2, 3 or 4 bonds, but here there is %d bonds, %s' % (nbonds,str(a)) 
                sys.exit(1)     
        # phosphor
        elif a.twoletter == 'P_':
            if nbonds == 4:
                a.ffType = 'P_4'
                a.lpair = 0
                valenceTest(a,'P',5)
            elif nbonds == 3:
                a.ffType = 'P_3'
                a.lpair = 1
                valenceTest(a,'P',3)
            else:
                print >> sys.stderr, 'Error: I only know P with 4 or 3 bonds, but here there is %d bonds, %s' % (nbonds,str(a)) 
                sys.exit(1)                
        # sulphur
        elif a.twoletter == 'S_':
            if nbonds == 4:
                a.ffType = 'S_4'
                a.lpair = 0
                valenceTest(a,'S',6)
            elif nbonds == 1 and a.hybridization == 2:
                a.ffType = 'S_2'
                a.lpair = 2
                valenceTest(a,'S',2)
            elif nbonds == 2 and a.hybridization == 2: # atom type requested by Adam 8/2011
                a.ffType = 'S_2'
                a.lpair = 1
                print >> sys.stderr, 'Warning: This is a rare case for sulphur: %s, typed as %s' % (str(a), a.ffType)  
            elif nbonds == 2 and a.hybridization == 3:
                a.ffType = 'S_3'
                a.lpair = 2
                valenceTest(a,'S',2)
            elif nbonds == 3 and a.hybridization == 2 and sorted(a.order) == [1,1,2]: # 8/1/2013 added for sulfoxides on adam's request 
                a.ffType = 'S_4'
                a.lpair = 1
                print >> sys.stderr, 'Warning: This is a rare case for sulphur: %s, typed as %s' % (str(a), a.ffType)  
                valenceTest(a,'S',4)
            else:
                print >> sys.stderr, 'Error: I do not know S with %d bonds and hybridization %d. %s' % (nbonds,a.hybridization,str(a)) 
                sys.exit(1)                
        # only C, O, N are left 
        # oxygen
        elif a.twoletter == 'O_':
            if a.hybridization == 1:
                if nbonds == 1:
                    a.ffType = 'O_1P'
                    a.lpair = 1
                    valenceTest(a,'O',3)
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
                    #if nbonds==3: 
                    #    a.ffType = 'O_RP'
                    #    a.lpair = 0
                    #el
                    if nbonds==2: 
                        a.ffType = 'O_RP'
                        a.lpair = 1
                        valenceTest(a,'O',3)
                    else:
                        print >> sys.stderr, 'Error: this resonant oxygen %s has too many bonds: %d ' % (str(a),nbonds) 
                        sys.exit(1)
                else: # O_2
                    if nbonds==2: 
                        a.ffType = 'O_2P'
                        a.lpair = 1
                        valenceTest(a,'O',3)
                    elif nbonds==1: 
                        # check if this is carboxylate or nitro or not
                        numberOfOxygensWithOneNeighbor = 0
                        for b in a.bonds[0].bonds:
                            if b.twoletter == 'O_' and len(b.bonds) == 1:
                                numberOfOxygensWithOneNeighbor += 1
                        if numberOfOxygensWithOneNeighbor == 2:
                            a.ffType = 'O_2M'
                            a.lpair = 2 
                            valenceTest(a,'O',1,upperBound=2)
                        elif numberOfOxygensWithOneNeighbor == 1: 
                            a.ffType = 'O_2'
                            a.lpair = 2
                            valenceTest(a,'O',2)
                        else:
                            print >> sys.stderr, 'Error: this sp2 oxygen %s is attached to atom with more than two oxygens ??? ' % (str(a)) 
                            sys.exit(1)
                    else:
                        print >> sys.stderr, 'Error: this sp2 oxygen %s should have 1 or 2 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
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
                    valenceTest(a,'O',2)
                elif nbonds == 1: # test for COOH
                    hasO_2neighbor = False
                    for b in a.bonds[0].bonds:
                        if b.twoletter == 'O_' and b.hybridization == 2:
                            hasO_2neighbor = True
                            break
                    if hasO_2neighbor:  
                        a.ffType = 'O_2M'  # O in COO
                        a.lpair = 2
                        valenceTest(a,'O',1,upperBound=2)
                    elif a.bonds[0].hybridization == 2:
                        a.ffType = 'O_RM'
                        a.lpair = 2
                        valenceTest(a,'O',1)
                    else:
                        a.ffType = 'O_3M'
                        a.lpair = 3      
                        valenceTest(a,'O',1)              
                elif nbonds == 2:
                    if a.bonds[0].twoletter == 'H_' and a.bonds[1].twoletter == 'H_':
                        a.ffType = 'O_3W'
                        a.lpair = 2
                        valenceTest(a,'O',2)
                    else:
                        a.ffType = 'O_3'
                        a.lpair = 2
                        valenceTest(a,'O',2)
                elif nbonds == 3:
                    a.ffType = 'O_3P'
                    a.lpair = 1
                    valenceTest(a,'O',3)
                else:
                    print >> sys.stderr, 'Error: this sp3 oxygen %s should have 1,2 or 3 bonds, but has %d ' % (str(a),nbonds) 
                    sys.exit(1)
        # carbon and nitrogen
        elif a.twoletter in ['C_','N_'] :
            if a.hybridization == 3:
                if a.twoletter == 'C_':
                    if nbonds == 4:
                        a.ffType = 'C_3'
                        a.lpair = 0
                        valenceTest(a,'C',4)
                    else:                                            
                        print >> sys.stderr, 'Error: sp3 carbon %s has to have 4 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                elif a.twoletter == 'N_':
                    if nbonds == 2:
                        a.ffType = 'N_3M'
                        a.lpair = 2
                        valenceTest(a,'N',2)
                    elif nbonds == 3:
                        a.ffType = 'N_3'
                        a.lpair = 1
                        valenceTest(a,'N',3)
                    elif nbonds == 4:
                        a.ffType = 'N_3P'
                        a.lpair = 0
                        valenceTest(a,'N',4)
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
                        valenceTest(a,'C',4)
                    elif nbonds == 1:
                        a.ffType = 'C_1'
                        a.lpair = 1
                        print 'Warning: sp1 carbon %s has only 1 bond, so I will assume it has one negative formal charge and one lone pair' % (str(a)) 
                    else:                                            
                        print >> sys.stderr, 'Error: sp1 carbon %s has to have 2 bonds, but has %d ' % (str(a),nbonds) 
                        sys.exit(1)
                elif a.twoletter == 'N_':
                    if nbonds == 1:
                        a.ffType = 'N_1'
                        a.lpair = 1
                        valenceTest(a,'N',3)
                    elif nbonds == 2:
                        a.ffType = 'N_1P'
                        a.lpair = 0
                        valenceTest(a,'N',4)
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
                            valenceTest(a,'C',4,upperBound=4.5)
                        else:                                            
                            print >> sys.stderr, 'Error: sp2 carbon %s has to have 3 bonds, but has %d ' % (str(a),nbonds) 
                            sys.exit(1)
                    elif a.twoletter == 'N_':
                        if nbonds == 2:
                            if a.order[0]<1.4 and a.order[1]<1.4:
                                a.ffType = 'N_RM'
                                a.lpair = 1
                                valenceTest(a,'N',2,upperBound=2.4)
                            else:
                                a.ffType = 'N_RL'
                                a.lpair = 1
                                valenceTest(a,'N',3)
                        elif nbonds == 3:
                            if a.order[0]<1.4 and a.order[1]<1.4 and a.order[2]<1.4:
                                a.ffType = 'N_R'
                                a.lpair = 0
                                valenceTest(a,'N',3,upperBound=3.4)
                            else:
                                a.ffType = 'N_RP'
                                a.lpair = 0
                                valenceTest(a,'N',4)
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
                            valenceTest(a,'C',4)
                        else:                                            
                            print >> sys.stderr, 'Error: sp2 carbon %s has to have 3 bonds, but has %d ' % (str(a),nbonds) 
                            sys.exit(1)
                    elif a.twoletter == 'N_':
                        if nbonds == 1:
                            a.ffType = 'N_2M'
                            a.lpair = 2
                            valenceTest(a,'N',2)
                        elif nbonds == 2:
                            a.ffType = 'N_2'
                            a.lpair = 1
                            valenceTest(a,'N',3)
                        elif nbonds == 3:
                            a.ffType = 'N_2P'
                            a.lpair = 0
                            valenceTest(a,'N',4)
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
            print >> sys.stderr, 'Error: There is no algorithm to type atom: "%d"' % a.twoletter 
            sys.exit(1)            

