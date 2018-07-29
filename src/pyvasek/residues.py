'''This modul includes tools to manipulate residues, mostly for amino acids

allowed N termini:
1. not completed     N, H
2. neutral           N, H, H2
3. charged           N, H1, H2, H3
(PRO does not have H and H1)

allowed C termini:
1. not completed     C, O
3. charged           C, O, OXT
2. neutral           C, O, HC      (ffType O2 used for HC and O for amber)
'''


#import re, copy, math, sys
import math, sys, numpy, copy, os, time
import atoms, geometry, defaults, bgf, pdb

#####################
# RESIDUE NAMES
#####################
name123 = {'A':'ALA',
    'R':'ARG',
    'N':'ASN',
    'D':'ASP', 
    'C':'CYS', 
    'E':'GLU',
    'Q':'GLN',
    'G':'GLY',
    'H':'HIS', 
    'I':'ILE', 
    'L':'LEU',
    'K':'LYS',
    'M':'MET',
    'F':'PHE',
    'P':'PRO',
    'S':'SER',
    'T':'THR',
    'W':'TRP',
    'Y':'TYR',
    'V':'VAL'}

name321 = {'ALA':'A', 
    'ARG':'R', 'ARN':'R',
    'ASN':'N',
    'ASP':'D', 'ASH':'D',
    'CYS':'C', 'CYX':'C',
    'GLU':'E', 'GLH':'E',
    'GLN':'Q',
    'GLY':'G',
    'HIS':'H', 'HID':'H', 'HIE':'H',  'HIP':'H',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K', 'LYN':'K',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V'}

name321withOld = {'ALA':'A', 
    'ARG':'R', 'ARN':'R',
    'ASN':'N',
    'ASP':'D', 'ASH':'D', 'APP':'D',
    'CYS':'C', 'CYX':'C',
    'GLU':'E', 'GLH':'E', 'GLP':'E',
    'GLN':'Q',
    'GLY':'G',
    'HIS':'H', 'HID':'H', 'HIE':'H', 'HIP':'H', 'HSE':'H', 'HSP':'H','HSD':'H',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K', 'LYN':'K',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V'}

name321Scream = name321.copy()
name321Scream['HID'] = 'H'
name321Scream['HIE'] = 'J'
name321Scream['HIP'] = 'B'
name321Scream['CYX'] = 'X'

oldBiogroupNames = {'GLP':'GLH',
    'APP':'ASH',
    'HIS':'HID',
    'HSE':'HIE',
    'HSP':'HIP'}

water = [ 'HOH', 'WAT', 'H3O', 'OH-' ]

knownResidues = name321.keys()
knownResidues.extend( water )

#####################
# ATOM NAMES
#####################
backbone = ['N','CA','C','O']
aaAtoms = {'ALA': ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O'],
    'ARG': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O'],
    'ARN': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'C', 'O'],
    'ASN': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O'],
    'ASP': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2', 'C', 'O'],
    'ASH': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O'],
    'CYS': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'HG', 'C', 'O'],
    'GLU': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'C', 'O'],
    'GLH': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O'],
    'GLN': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O'],
    'GLY': ['N', 'H', 'CA', 'HA2', 'HA3', 'C', 'O'],
    'HID': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'C', 'O'],
    'HIE': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'HE2', 'C', 'O'],
    'HIP': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'HE2', 'C', 'O'],
    'ILE': ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'CD1', 'HD11', 'HD12', 'HD13', 'C', 'O'],
    'LEU': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O'],
    'LYS': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O'],
    'LYN': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ2', 'HZ3', 'C', 'O'],
    'MET': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O'],
    'PHE': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'HZ', 'C', 'O'],
    'PRO': ['N', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'C', 'O'],
    'SER': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OG', 'HG', 'C', 'O'],
    'THR': ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'OG1', 'HG1', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O'],
    'TRP': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'NE1', 'HE1', 'CE2', 'CE3', 'HE3', 'CZ2', 'HZ2', 'CZ3', 'HZ3', 'CH2', 'HH2', 'C', 'O'],
    'TYR': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'OH', 'HH', 'C', 'O'],
    'VAL': ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O'],
    'CYX': ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'C', 'O'],
    'WAT': ['O','H1','H2'],
    'HOH': ['O','H1','H2'],
    'H3O': ['O','H1','H2','H3'],
    'OH-': ['O','H1']}


orderGreek = {'A':1,'B':2,'G':3,'D':4,'E':5,'Z':6,'H':7}
orderSpecial = {'N': -10,
    'H': -9,
    'H1': -8,
    'H2': -7,
    'H3': -6,
    'HC': 102,
    'OXT': 103}   


class Residue:
    """residue object stores information for one residues"""
    def __init__(self):
        """creates a blank residue"""
        self.rName  = ''      # residue name
        self.chain  = ''      # chain name
        self.rNo    = 0       # residue number
        self.atoms  = []      # list of atoms

    def __str__(self):
        """Short residue name with chain"""
        return "%s%d:%s" % (self.rName, self.rNo, self.chain) 

    def __getitem__(self, index):
        return self.atoms[index]

    def __len__(self):
        return len(self.atoms)

    def __cmp__(self, other):
        if other.__class__.__name__ == 'Residue':
            if self.chain < other.chain:
                return -1
            elif self.chain > other.chain:
                return 1
            elif self.rNo < other.rNo:
                return -1
            elif self.rNo > other.rNo:
                return 1
            else:
                return 0
        else:
            return -1

    def initialize(self):
        '''Initilizes attributes from atoms stored in atoms list'''
        if len(self.atoms) != 0:
            self.rName = self.atoms[0].rName.strip()
            self.chain = self.atoms[0].chain
            self.rNo = self.atoms[0].rNo
        else:
            print 'Warning: atoms list is empty in residue initialize!'

    def updateAtoms(self):
        '''updates rName, rNo and chain for all atoms in this residue'''
        for atom in self.atoms:
            atom.rNo = self.rNo
            atom.rName = self.rName
            atom.chain = self.chain

    def addAtom(self, atom):
        self.atoms.append(atom)
        atom.rNo = self.rNo
        atom.rName = self.rName
        atom.chain = self.chain

    def disconnect(self):
        atoms.disconnectExternalBonds(self.atoms)

    # charge
    def charge(self):
        charge = 0
        for atom in self.atoms:
            charge += atom.charge
        if math.fabs(charge) <= 1e-10:
            return 0.0
        else:
            return charge

    # functions based on residue name
    def screamLetter(self):
        if self.rName in name321Scream:
            return name321Scream[self.rName]
        else:
            return '?'

    def screamString(self):
        aa = self.screamLetter()
        if aa == '?':
            return '?'
        else:
            return aa + str(self.rNo) + '_' + self.chain    

    def oneLetter(self):
        if self.rName in name321:
            return name321[self.rName]
        else:
            return '?'

    def isAminoAcid(self):
        if self.rName in name321:
            return True
        else:
            return False

    def isPolar(self):
        aa = self.screamLetter()
        polar_res = ['R','H','B','J','K','D','E','S','T','N','Q','Y']
        if aa in polar_res:
            return True
        else:
            return False

    def canonicalName(self):
        try:
            canon_name = name123[self.oneLetter()]
        except KeyError, e:
            print 'Warning: in canonicalName(): residue %s is not amino acid, so cannot assign canonical name'
            canon_name = self.rName
        return canon_name

    # atom look up
    def atomWithName(self, atom_label):
        '''Finds first instance of atom with atom label name. 
        Does not guard against residues with atoms with multiple atom labels. 
        Returns -1 if there is no such atom'''
        atom_label = atom_label.strip()        
        for atom in self.atoms:
            if atom.aName.strip() == atom_label:
                return atom
        return -1

    def listAtomsWithName(self, atom_label):
        '''Finds all atoms with atom label name.  Output in a list.'''
        atom_label = atom_label.strip()         
        mylist = []
        for atom in self.atoms:
            if atom.aName.strip() == atom_label:
                mylist.append(atom)
        return mylist
    # finding new chain
    def isConnectedAtN(self, res0):
        '''checks if atom C of residue res0 makes bond to atom N of this residue'''
        n = self.atomWithName('N')
        c = res0.atomWithName('C')
        if n !=-1 and c!=-1:
            answer  = n.isBonded(c)
            return answer
        return False

    def NisCloseToC(self, res0):
        '''checks if atom C of residue res0 is close enough to atom N of this residue to make a bond
        i.e. if they are closer than 1.6 A'''
        n = self.atomWithName('N')
        c = res0.atomWithName('C')
        if n !=-1 and c!=-1:
            answer = n.distance(c) < 1.6
            return answer
        return False

    # advanced amino acid functions
    def sort(self):
        '''sorts atoms in this residue to be in the standard order (amino acids and water only)'''
        if self.rName in knownResidues:
            self.atoms.sort( key = atomOrder )

    def checkAtomPresence(self,position='body',onlyHeavy=False):
        '''returns two tuplet of sets and one integer: (extraAtomNames, missingAtomNames, dupliciteAtomNumber)
        prints if there are missign or etra atoms
        position is 'body', 'Nterminus', 'Cterminus', 'Nterminus+','Cterminus-' 
        only works for known amino acids and water'''
        r = self
        if r.rName not in knownResidues:
            sys.exit('Error: Residue %s is not amino acid or water so I cannot perform atom presence test' % (r.rName))
        listthis = [a.aName for a in r.atoms]
        setthis = set(listthis)
        duplicite = len(listthis) - len(setthis)
        if duplicite > 0:
            print 'Error: residue %s has duplicite atoms' % r
        setref = set(aaAtoms[r.rName])

        if position == 'body':
            pass
        elif position == 'Cterminus':
            setref |= set(['HC'])
        elif position == 'Cterminus-':
            setref |= set(['OXT'])
        elif position == 'Nterminus':
            setref |= set(['H2'])
        elif position == 'Nterminus+':
            if r.rName == 'PRO':
                setref |= set(['H2','H3'])
            else:
                setref -= set(['H'])
                setref |= set(['H1','H2','H3'])
        else:
            sys.exit('Error: Unknown value %s of option position in function checkAtomPresence' % (position))

        missing = setref - setthis
        extra = setthis - setref
        
        if onlyHeavy:
            for atomname in list(extra):
                if atomname.startswith('H'):
                    extra.remove(atomname)
            for atomname in list(missing):
                if atomname.startswith('H'):
                    missing.remove(atomname)

        if extra != set([]) or missing != set([]):
            print 'Warning: residue: %s has   extra atoms: %s   missing atoms:   %s  (position: %s)' % (r,extra,missing, position)
        return extra, missing, duplicite
    
    def checkBonds(self,store, position='body'):
        '''checks if there any extra or any missigns bonds 
        returns True if there is no error
        position is 'body', 'Nterminus', 'Cterminus', 'Nterminus+','Cterminus-' 
        position is irrelevant for water'''
        isOK = True
        r = self
        if r.rName not in knownResidues:
            print >> sys.stderr, 'Error: Residue %s is not amino acid or water ' % (r.rName)
            sys.exit(1)
        else: # r is water or amino acid
            for a in r.atoms:
                bondListIn = [b.aName.strip() for b in a.bonds if b in r.atoms]
                bondListOut = [b.aName.strip() for b in a.bonds if b not in r.atoms]
                bondListIn.sort()
                bondListOut.sort()
                # resolve possible cysteine bridge
                if a.rName == 'CYX' and a.aName.strip() == 'SG':
                    if 'SG' not in bondListOut:
                        isOK = False
                        print 'Error: this residues %s is CYX, but cysteine bridge is not formed from atom %s' % (r, a)
                    else:
                        bondListOut.pop(bondListOut.index('SG'))
                # resolve N terminal
                elif a.aName.strip() == 'N':
                    if position in ['body','Cterminus','Cterminus-']:
                        if bondListOut != ['C']:
                            isOK = False
                            print 'Error: %s is connected to these atoms outside of the residue %s: %s' % (a, r, bondListOut)
                        else:
                            bondListOut.pop(bondListOut.index('C'))
                    elif position == 'Nterminus':
                        if 'H2' not in bondListIn:
                            isOK = False
                            print 'Error: %s is not connected to H2' % (a)
                        else:
                            bondListIn.pop(bondListIn.index('H2'))
                    elif position == 'Nterminus+':
                        if ('H2' not in bondListIn) or ('H3' not in bondListIn):
                            isOK = False
                            print 'Error: %s is not connected to H2 or H3' % (a)
                        else:
                            bondListIn.pop(bondListIn.index('H2'))
                            bondListIn.pop(bondListIn.index('H3'))
                        if a.rName != 'PRO':
                            if 'H1' not in bondListIn:
                                isOK = False
                                print 'Error: %s is not connected to H1' % (a)
                            else:
                                bondListIn.pop(bondListIn.index('H1'))
                                bondListIn.append('H')
                                bondListIn.sort()
                    else:
                        sys.exit('Error: unexpected position keyword in residues.Residue.checkBonds: %s' % position)
                # resolve C terminal
                elif a.aName.strip() == 'C':
                    if position in ['body','Nterminus','Nterminus+']:
                        if bondListOut != ['N']:
                            isOK = False
                            print 'Error: %s is connected to these atoms outside of the residue %s: %s' % (a, r, bondListOut)
                        else:
                            bondListOut.pop(bondListOut.index('N'))
                    elif position == 'Cterminus':
                        if 'HC' not in bondListIn:
                            isOK = False
                            print 'Error: %s is not connected to HC' % (a)
                        else:
                            bondListIn.pop(bondListIn.index('HC'))
                    elif position == 'Cterminus-':
                        if 'OXT' not in bondListIn:
                            isOK = False
                            print 'Error: %s is not connected to OXT' % (a)
                        else:
                            bondListIn.pop(bondListIn.index('OXT'))
                    else:
                        sys.exit('Error: unexpected position keyword in residues.Residue.checkBonds: %s' % position)
                else:
                    lookUpName = a.aName.strip()
                    if a.rName in name321:
                        if lookUpName in ['H1','H2','H3']:
                            lookUpName = 'H'
                        if lookUpName in ['OXT','HC']:
                            lookUpName = 'O'
                    if a.rName == 'PRO' and lookUpName == 'H':
                        if bondListIn != ['N']:
                            isOK = False
                            print 'Error: atom %s is not only connected to N but to: %s' %(a,bondListOut)
                    elif store[r.rName+lookUpName] != bondListIn:
                        isOK = False
                        setRef = set(store[r.rName+lookUpName])
                        setReal = set(bondListIn)
                        missing = setRef - setReal
                        extra = setReal - setRef
                        if len(missing) > 0:
                            print 'Error: atom %s is not connected to atoms: %s' % (a, missing)
                        if len(extra) > 0:
                            print 'Error: atom %s is connected to extra atoms: %s' % (a, missing)
                if len(bondListOut) > 0:
                    isOK = False
                    print 'Error: atom %s is connected to extra atoms outside of the residue %s: %s' %(a,r,bondListOut)
        return isOK

    def makeBonds(self,store, position='body'):
        '''make bonds within this residue based on the database 'store'
        position is 'body', 'Nterminus', 'Cterminus', 'Nterminus+','Cterminus-' 
        position is irrelevant for water

        allowed N termini:
        1. not completed     N, H
        2. neutral           N, H, H2
        3. charged           N, H1, H2, H3
        (PRO does not have H and H1)

        allowed C termini:
        1. not completed     C, O
        3. charged           C, O, OXT
        2. neutral           C, O, HC      (ffType O2 used for HC and O for amber)        
        '''
        specialN = 'body', 'Nterminus', 'Cterminus', 'Nterminus+'
        isOK = True
        r = self
        if r.rName not in knownResidues:
            print >> sys.stderr, 'Error: Residue %s is not amino acid or water ' % (r.rName)
            sys.exit(1)
        else: # r is water or amino acid
            for a in r.atoms:
                if a.aName[0:1] != 'H': # hydrogens will get always attached since there is no H-H bond
                    # correction for terminals
                    if a.aName == 'N':
                        atomstore = store[r.rName][a.aName][:]
                        if position == 'Nterminus':
                            atomstore.append('H2')
                        elif position == 'Nterminus+':
                            if r.rName == 'PRO':
                                atomstore.extend(['H1','H2'])
                            else:
                                atomstore.remove('H')
                                atomstore.extend(['H1','H2','H3'])
                    elif a.aName == 'C':
                        atomstore = store[r.rName][a.aName][:]
                        if position == 'Cterminus':
                            atomstore.append('HC')
                        elif position == 'Cterminus-':
                            atomstore.append('OXT')
                    else:
                        atomstore = store[r.rName][a.aName]

                    # make the bonds
                    for bname in atomstore:
                        b = r.atomWithName(bname)
                        if not a.isBonded(b):
                            a.makeBond(b)

    def assignFFCharge(self,store,position='body'):  
        for a in self.atoms:
            if a.aName not in ['HC','OXT','H1','H2','H3']:
                value = store[self.rName+a.aName]
                a.ffType = value[0]
                a.charge = value[1]
            else:
                a.ffType = ''
                a.charge = 0.0
        if position == 'body':
            pass
        elif position == 'Cterminus':
            a = self.atomWithName('HC')
            value = store['ALAN']
            if value[0] == 'N': # amber
                a.ffType = 'O2'
                self.atomWithName('O').ffType = 'O2'
            else:               # dreiding
                a.ffType = 'H_'
            a.charge = 0.0     
        elif position == 'Cterminus-':
            h1 = self.atomWithName('O')
            h2 = self.atomWithName('OXT')
            # ffType
            value = store['ALAO']
            h1.ffType = value[0]
            h2.ffType = value[0]
            # charge                
            ch = self.charge()-h1.charge-h2.charge
            dch = -2.0 - ch + round(ch)  # I had -1.0 -ch+rouch(ch) here, 
                # but charge of O is between -0.5 and -1.0 so I made it -2.0 - ch + round(ch) 
            partch = round(dch*10000.0/2.0)/10000.0
            h1.charge = dch-partch
            h2.charge = partch
        elif position == 'Nterminus':
            if self.rName == 'PRO':
                h2 = self.atomWithName('H2')
                value = store['ALAH']
                h2.ffType = value[0]
                h2.charge = 0.0
            else:
                h1 = self.atomWithName('H')
                h2 = self.atomWithName('H2')
                # ffType
                value = store['ALAH']
                h1.ffType = value[0]
                h2.ffType = value[0]
                # charge                
                ch = self.charge()-h1.charge-h2.charge
                dch = - ch + round(ch)
                partch = round(dch*10000.0/2.0)/10000.0
                h1.charge = dch-partch
                h2.charge = partch
        elif position == 'Nterminus+':
            if self.rName == 'PRO':
                h2 = self.atomWithName('H2')
                h3 = self.atomWithName('H3')
                value = store['ALAH']
                h2.ffType = value[0]
                h3.ffType = value[0]
                h2.charge = +0.5
                h3.charge = +0.5
            else:
                h1 = self.atomWithName('H1')
                h2 = self.atomWithName('H2')
                h3 = self.atomWithName('H3')
                # ffType
                value = store['ALAH']
                h1.ffType = value[0]
                h2.ffType = value[0]
                h3.ffType = value[0]
                # charge                
                ch = self.charge()-h1.charge-h2.charge-h3.charge
                dch = 1.0 - ch + round(ch)
                partch = round(dch*10000.0/3.0)/10000.0
                h1.charge = dch-2*partch
                h2.charge = partch
                h3.charge = partch

    # termination 
    def terminateN(self,mode='neutral'):
        '''the residue has to end with a dangling bond from N, this does not attach force field type or charge''' 
        if mode == 'neutral':
            a = atoms.Atom()
            a.aTag = 'ATOM  '
            a.aName = 'H2'
            bond = 0.97            
            nitrogen = self.atomWithName('N')
            a.makeBond(nitrogen)
            self.addAtom(a)
            # put the atom to the right position
            if self.rName == 'PRO':
                otheratom = self.atomWithName('CD')
            else:
                otheratom = self.atomWithName('H')
                if otheratom == -1:
                    # there is no H in this residue, let's just make it
                    h = atoms.Atom()
                    h.aTag = 'ATOM  '
                    h.aName = 'H'
                    ca = self.atomWithName('CA')
                    c = self.atomWithName('C')
                    h.xyz = geometry.sp3given1(nitrogen.xyz,ca.xyz,c.xyz,bond)[0]
                    nitrogen.makeBond(h)
                    self.addAtom(h)
                    otheratom = h
            self.sort()
            a.xyz = geometry.sp3given2(nitrogen.xyz,self.atomWithName('CA').xyz,otheratom.xyz,bond)[0]
        elif mode == 'charged':
            a = atoms.Atom()
            a.aTag = 'ATOM  '
            a.aName = 'H2'
            b = atoms.Atom()
            b.aTag = 'ATOM  '
            b.aName = 'H3'
            bond = 0.97
            
            nitrogen = self.atomWithName('N')
            a.makeBond(nitrogen)
            b.makeBond(nitrogen)
            self.addAtom(a)
            self.addAtom(b)
            self.sort()
            # put the atom to the right position
            if self.rName == 'PRO':
                otheratom = self.atomWithName('CD')
            else:
                otheratom = self.atomWithName('H')
                otheratom.aName = 'H1'
            pos = geometry.sp3given2(nitrogen.xyz,self.atomWithName('CA').xyz,otheratom.xyz,bond)
            a.xyz = pos[0]
            b.xyz = pos[1]
        else:
            print 'Error: unknown mode %s' % mode
            sys.exit(1)

    def terminateC(self,mode='neutral'):
        '''the residue has to end with a dangling bond from C'''
        carbon = self.atomWithName('C')
        if mode == 'neutral': # creates HC with charge 0, type H_, and C-HC bond of length ....
            a = atoms.Atom()
            a.aTag = 'ATOM  '
            a.aName = 'HC'
            bond = 1.22       # H-C bond should be 0.99, but this will be most likely used for amber,
                              # where I use oxygen type for HC, so the bond length is set to 1.22
            a.makeBond(carbon)
        elif mode == 'charged': # creates OXT with charge -1, type O..., and C-OXT bond of length ....
            a = atoms.Atom()
            a.aTag = 'ATOM  '
            a.aName = 'OXT'
            bond = 1.22
            a.makeBond(carbon,bondOrder = 1.5)
        else:
            print 'Error: unknown mode %s' % mode
            sys.exit(1)
        self.addAtom(a)
        self.sort()
        # put the atom to the right position
        a.xyz = geometry.sp2given2(carbon.xyz,self.atomWithName('O').xyz,self.atomWithName('CA').xyz,bond)

    def unterminate(self):
        # c terminal first
        a = self.atomWithName('HC')
        if a != -1:
            a.delete()
            self.atoms.remove(a)
        a = self.atomWithName('OXT')
        if a != -1:
            a.delete()
            self.atoms.remove(a)
        a = self.atomWithName('O')
        if a != -1:
            ch = self.charge()-a.charge
            dch = round(ch)-ch
            a.charge = dch
        # n terminal
        a = self.atomWithName('H3')
        if a != -1:
            a.delete()
            self.atoms.remove(a)
        a = self.atomWithName('H2')
        if a != -1:
            a.delete()
            self.atoms.remove(a)
        a = self.atomWithName('H1')
        if a != -1:
            a.aName = 'H'
        a = self.atomWithName('H')
        if a != -1:
            ch = self.charge()-a.charge
            dch = round(ch) - ch
            a.charge = dch
    
    # geometry of amino acids    
    def isCis(self,r2):
        ''' CA, C takend from self and N, CA taken from r2
        99.9% of amino acids are trans, 60-90% of PRO is trans
        trans... dihedral 180
        cis ... dihedral 0'''
        CA = self.atomWithName('CA')
        C  = self.atomWithName('C')
        CA2= r2.atomWithName('CA')
        N  = r2.atomWithName('N')
        if CA==-1 or C==-1:
            print 'Error: could not find backbone atoms to determine Cis/Trans on %s' % self
            sys.exit(1)
        if CA2==-1 or N==-1:
            print 'Error: could not find backbone atoms to determine Cis/Trans on %s' % r2
            sys.exit(1)        
        a1 = CA.xyz-C.xyz                 
        a2 = CA2.xyz-N.xyz                 
        a3 = N.xyz-C.xyz                 
        a3 = a3 / numpy.linalg.norm(a3)
        a1 = a1 - numpy.dot(a1,a3) * a3
        a2 = a2 - numpy.dot(a2,a3) * a3
        cos = numpy.dot(a1,a2) / numpy.linalg.norm(a1) / numpy.linalg.norm(a2)
        if cos < 0:
            return False
        else:
            return True

    def LorD(self):
        ''' returns L,D,G (for glycine) or ? 
        amino acids in proteins are L'''       
        if self.rName == 'GLY':
            return 'G'
        CA = self.atomWithName('CA')
        C  = self.atomWithName('C')
        CB = self.atomWithName('CB')
        N  = self.atomWithName('N')
        HCA= self.atomWithName('HA')
        if CA==-1 or C==-1 or CB==-1 or N==-1 or HCA==-1:
            print 'Error: could not find atoms to determine D/L on %s' % self
            return '?'
        b1=  C.xyz-CA.xyz
        b2= CB.xyz-CA.xyz
        b3=  N.xyz-CA.xyz
        b4=HCA.xyz-CA.xyz
        c = numpy.cross(b1,b2)
        if numpy.linalg.norm(c)<1e-9:
            return '?'
        c3= numpy.dot(c,b3)
        c4= numpy.dot(c,b4) 
        if c3 > 0.0 and c4 < 0.0:
          return 'L'
        elif c4 > 0.0 and c3 < 0.0:
          return 'D'
        else:
          return '?'

################# end of class Residue ##################

def makeList(al):
    '''make residue list from the list of atoms
    assume only atoms that are immediatelly following each other in the list belong to the same residue '''
    
    if len(al) == 0:
        return []
    rl = [Residue()] # the new residue list

    r = rl[0]
    a = al[0]
    r.rName = a.rName
    r.chain = a.chain
    r.rNo = a.rNo
    r.atoms.append(a)
    
    for a in al[1:]:
        if a.chain == r.chain and a.rNo == r.rNo and a.rName == r.rName:
            r.atoms.append(a)
        else:
            r = Residue()
            r.rName = a.rName
            r.chain = a.chain
            r.rNo = a.rNo
            r.atoms.append(a)
            rl.append(r)    
    return rl

def makeAtomList(rl):
    al = []
    for r in rl:
        al.extend(r.atoms)
    return al

def getResidue(rl, resno, chain):
    '''looks for residue with given number and chain from the residus list
    returns -1 if there is no such residue'''
    for r in rl:
        if r.rNo == resno and r.chain == chain:
            return r
    return -1

# atom ordering 

def atomOrder(a):
    if a.rName in aaAtoms:
        if a.aName in orderSpecial and a.rName not in water:
            return orderSpecial[a.aName]
        elif a.aName in aaAtoms[a.rName]:
            return aaAtoms[a.rName].index(a.aName)
        else:
            print >> sys.stderr, 'Error: Residue %s %s %s does not contain atom name: %s' % (a.rName, a.chain, a.rNo, a.aName)
            sys.exit(1)  
    else:
        print >> sys.stderr, 'Error: Can only sort amino acids or water, but this is: %s' % a.rName
        sys.exit(1)    

# atom naming

def namesOldToNew(al):
    '''translates from old biogroup atom names'''
    rl = makeList(al)
    for r in rl:
        # fix residues 
        if r.rName == 'GLP':
            h = r.atomWithName('HOE1')
            if h != -1:
                o1 = r.atomWithName('OE1')
                o2 = r.atomWithName('OE2')
                if o1 == -1 or o2 == -1:
                    sys.exit('Error: residue %s does not have both oxygens OE1 and OE2' % r)
                h.aName = 'HOE2'
                o1.aName = 'OE2'
                o2.aName = 'OE1'
        elif r.rName == 'APP':
            h = r.atomWithName('HOD1')
            if h != -1:
                o1 = r.atomWithName('OD1')
                o2 = r.atomWithName('OD2')
                if o1 == -1 or o2 == -1:
                    sys.exit('Error: residue %s does not have both oxygens OD1 and OD2' % r)
                h.aName = 'HOD2'
                o1.aName = 'OD2'
                o2.aName = 'OD1'
        elif r.rName == 'ARN':
            h22 = r.atomWithName('HH22') 
            if h22 != -1:
                h21 = r.atomWithName('HH21')
                h1 = r.atomWithName('HH11')
                if h1 == -1:
                    h1 = r.atomWithName('HNH1')
                n1 = r.atomWithName('NH1')
                n2 = r.atomWithName('NH2')
                if h21 == -1 or h1 == -1 or n1 == -1 or n2 == -1:
                    sys.exit('Error: residue %s is messed up (atoms: %s)' % (r,[a.aName for a in r.atoms]))
                n1.aName = 'NH2'
                n2.aName = 'NH1'
                h1.aName = 'HH21'
                h21.aName = 'HH11'
                h22.aName = 'HH12'        

        if r.rName in oldBiogroupNames:
            # change the residue name to the new one
            r.rName = oldBiogroupNames[r.rName]
            r.updateAtoms()

        # fix atoms 
        if r.rName in name321:
            counts = {}
            for a in r.atoms:
                if a.aName[0] == 'H':
                    if a.aName not in counts: counts[a.aName] = 0
                    counts[a.aName] += 1
            used = counts.copy()
            for i in counts:
                j = counts[i]
                if j == 2:
                    used[i] = 2
                elif j == 3:
                    used[i] = 1
            for a in r.atoms:
                if a.aName[0] == 'H':
                    # special cases first
                    # terminals
                    if a.aName == 'HN' and counts[a.aName] == 2:
                        if used[a.aName] == 2:
                            used[a.aName] += 1
                            a.aName = 'H'
                        else:
                            a.aName = 'H2'
                    elif a.aName == 'HC':
                        pass
                    # ARN
                    elif r.rName == 'ARN' and a.aName[0:2] == 'HH':
                        pass
                    elif r.rName == 'ARN' and a.aName[0:3] == 'HNH' and counts[a.aName] == 1:
                        a.aName = 'HH'+a.aName[3:4]+'1'
                    # PRO
                    elif r.rName == 'PRO' and a.aName == 'HN':
                        a.aName = 'H2'
                    # general case
                    elif counts[a.aName] == 1:
                        a.aName = a.aName[0]+a.aName[2:] 
                    else:
                        index = used[a.aName]
                        used[a.aName] += 1
                        if a.aName[1] == 'N' and len(a.aName[2:]) == 2:
                            index -= 1
                        a.aName = a.aName[0]+a.aName[2:]+'%d'%index

def namesNewToOld(al):
    '''translates to old biogroup atom names'''
    rl = makeList(al)
    for r in rl:
        # fix atoms
        if r.rName in name321:
            heavy = {}
            for a in r.atoms: 
                if a.aName[0] != 'H':
                    heavy[a.aName[1:]] = a.aName[0]
            for a in r.atoms:
                if a.aName[0] == 'H':
                    if len(a.aName) == 1:
                        a.aName = 'HN'
                    elif a.aName in ['H1','H2','H3']:
                        a.aName = 'HN'
                    elif a.aName in ['HC']:
                        a.aName = 'HC'                
                    elif len(a.aName) == 2:
                        a.aName = 'H'+heavy[a.aName[1:2]]+a.aName[1:2]
                    elif len(a.aName) == 3 and a.aName[1:] in heavy:
                        a.aName = 'H'+heavy[a.aName[1:3]]+a.aName[1:3]
                    elif len(a.aName) == 3:
                        a.aName = 'H'+heavy[a.aName[1:2]]+a.aName[1:2]
                    elif len(a.aName) == 4:
                        a.aName = 'H'+heavy[a.aName[1:3]]+a.aName[1:3]
                    else:
                        print 'Error: I cannot convert new atom name %s to old name in residue %s' % (a.aName,r) 
            # fix residue name
            oldValues = oldBiogroupNames.values()
            if r.rName in oldValues:
                oldKeys = oldBiogroupNames.keys()
                r.rName = oldKeys[oldValues.index(r.rName) ]
                r.updateAtoms()    

# checking integrity of proteins
def checkAtomPresence(al,termini='neutral',distanceTestC=False):
    '''returns True if all the atom names are fine,
    otherwise returns False and prints the type of the problem
    if distanceTestC = True the chain breaks are not determined from connectivity from from the distance to the next C atom'''
    if termini == 'neutral':
        modeN = 'Nterminus'
        modeC = 'Cterminus'
    elif termini == 'charged':
        modeN = 'Nterminus+'
        modeC = 'Cterminus-'
    elif termini == 'none':
        modeN = 'body'
        modeC = 'body'
    else:
        print >> sys.stderr, "Error: termini = '%s' but only: neutral, charged and none allowed" % termini
        sys.exit(1)    
    isOK = True
    rl = makeList(al)
    for ir in range(len(rl)):
        r = rl[ir]
        if r.rName in name321:
            a = r.atomWithName('N')        
            if a == -1:
                print "Error: residue %s misses atom N" % r
                isOK = False
            else:
                mode = modeN              
                if distanceTestC:
                    if ir >= 1:
                        prevr = rl[ir-1]
                        if r.NisCloseToC(prevr):
                            mode = 'body'
                else:
                    for b in a.bonds:
                        if b.aName == 'C':
                            mode = 'body'
                            break
                if mode == modeN:
                    extra, missing, duplicite = r.checkAtomPresence(position=mode)
                    if extra!=set([]) or missing!=set([]) or duplicite>0:
                        isOK = False
                else:
                    a = r.atomWithName('C')        
                    if a == -1:
                        print "Error: residue %s misses atom C" % r
                        isOK = False
                    else:
                        mode = modeC
                        if distanceTestC:
                            if ir < len(rl)-1:
                                nextr = rl[ir+1]
                                if nextr.NisCloseToC(r):
                                    mode = 'body'
                        else:
                            for b in a.bonds:
                                if b.aName == 'N':
                                    mode = 'body'
                                    break
                        extra, missing, duplicite = r.checkAtomPresence(position=mode)
                        if extra!=set([]) or missing!=set([]) or duplicite>0:
                            isOK = False
        elif r.rName in knownResidues:
            extra, missing, duplicite = r.checkAtomPresence()
            if extra!=set([]) or missing!=set([]) or duplicite>0:
                isOK = False
        else:
            print 'Warning: skipping atom presence test for residue %s' % r
    return isOK

def checkAtomPresenceHeavy(al):
    '''returns True if all the heavy atoms are present,
    otherwise returns False and prints the type of the problem;
    this test does not care about termination'''
    isOK = True
    rl = makeList(al)
    stringScwrl = ''
    for r in rl:
        if r.rName in name321:
            if r.rName == 'HIS':
                r.rName = 'HID'
            extra, missing, duplicite = r.checkAtomPresence(position='body',onlyHeavy=True)
            if extra!=set([]) or missing!=set([]) or duplicite>0:
                isOK = False
                stringScwrl += name321[r.rName].lower()
            else:
                stringScwrl += name321[r.rName]
        else:
            print 'Warning: skipping atom presence test for residue %s' % r

    return isOK,stringScwrl

def checkBonds(al,termini='neutral'):
    '''this assumes that atom presence test was passed
    returns True if all bonds are connected ok
    for ligans or unknown residues check if there are bonds only within that residue'''
    # create reference table    
    refBgf = defaults.refBgf['dreiding3-ff99SB']
    refal = bgf.read(refBgf)
    refrl = makeList(refal)
    store = {}    
    for r in refrl:
        if r.rName in knownResidues:
            for a in r.atoms:
                store[a.rName+a.aName.strip()] = [b.aName.strip() for b in a.bonds if b in r.atoms]
                store[a.rName+a.aName.strip()].sort()
    # treatment of the terminals
    if termini == 'neutral':
        modeN = 'Nterminus'
        modeC = 'Cterminus'
    elif termini == 'charged':
        modeN = 'Nterminus+'
        modeC = 'Cterminus-'
    elif termini == 'none':
        modeN = 'body'
        modeC = 'body'
    else:
        print >> sys.stderr, "Error: termini = '%s' but only: neutral, charged and none allowed" % termini
        sys.exit(1)   
    rl = makeList(al)
    isOK = True
    for r in rl:
        if r.rName in name321:
            a = r.atomWithName('N')        
            if a == -1:
                print "Error: residue %s misses atom N" % r
            else:
                mode = modeN      
                for b in a.bonds:
                    if b.aName == 'C':
                        mode = 'body'
                        break    
                if mode == modeN:
                    if not r.checkBonds(store,position=mode):
                        isOK = False
                else:
                    a = r.atomWithName('C')        
                    if a == -1:
                        print "Error: residue %s misses atom C" % r
                    else:
                        mode = modeC              
                        for b in a.bonds:
                            if b.aName == 'N':
                                mode = 'body'
                                break
                        if not r.checkBonds(store,position=mode):
                            isOK = False
        elif r.rName in knownResidues:
            if not r.checkBonds(store):
                isOK = False
        else: # for ligands only checking if there are bonds only within this residue
            for a in r.atoms:
                for b in a.bonds:
                    if b not in r.atoms:
                        print 'Warning: atoms %s and %s are connected by a bond, but they are not in the same residue' % (a,b)
                        isOK = False
    return isOK            

def makeBonds(al,termini='neutral'):
    '''this assumes that atom presence test was passed
    returns True if all bonds are connected ok
    for ligans or unknown residues check if there are bonds only within that residue'''
    # create reference table    
    refBgf = defaults.refBgf['dreiding3-ff99SB']
    refal = bgf.read(refBgf)
    refrl = makeList(refal)
    store = {}    
    for r in refrl:
        if r.rName in knownResidues and r.rName not in store:
            resstore = {}
            for a in r.atoms:
                resstore[a.aName] = [b.aName for b in a.bonds if b in r.atoms]
            store[r.rName] = resstore
    # treatment of the terminals
    if termini == 'neutral':
        modeN = 'Nterminus'
        modeC = 'Cterminus'
    elif termini == 'charged':
        modeN = 'Nterminus+'
        modeC = 'Cterminus-'
    elif termini == 'none':
        modeN = 'body'
        modeC = 'body'
    else:
        sys.exit("Error: termini = '%s' but only: neutral, charged and none allowed" % termini)
    rl = makeList(al)
    for ir in range(len(rl)):
        r = rl[ir]
        if r.rName in name321:
            a = r.atomWithName('N')        
            if a == -1:
                sys.exit("Error: residue %s misses atom N" % r)
            mode = modeN     
            if ir >= 1:
                prevr = rl[ir-1]
                if r.NisCloseToC(prevr):
                    mode = 'body'
                    n = a
                    c = prevr.atomWithName('C')
                    if not n.isBonded(c):
                        n.makeBond(c)
            if mode == modeN:
                r.makeBonds(store,position=mode)
            else:
                a = r.atomWithName('C')        
                if a == -1:
                    sys.exit("Error: residue %s misses atom C" % r)
                mode = modeC        
                if ir < len(rl)-1:
                    nextr = rl[ir+1]
                    if nextr.NisCloseToC(r):
                        mode = 'body'
                        n = nextr.atomWithName('N')
                        c = a
                        if not n.isBonded(c):
                            n.makeBond(c)
                r.makeBonds(store,position=mode)
        elif r.rName in knownResidues:
            r.makeBonds(store)

def assignFFCharge(al, ff = 'dreiding3', charges = 'ff99SB', termini='neutral', skipErrors=False):
    '''note: if residue is both on terminal C and N (chain of length 1) treat it as terminal N'''
    # create reference table    
    refBgf = defaults.refBgf['%s-%s'%(ff,charges)]
    refal = bgf.read(refBgf)
    store = {}    
    for a in refal:
        if a.rName in knownResidues:
            store[a.rName+a.aName] = (a.ffType.strip(),a.charge)
    # check integrity of the file
    if not checkAtomPresence(al,termini=termini):
        if not skipErrors:
          sys.exit('Error: the atom list did not pass atom presence test, so I quit (in assignFFCharge)')
        else:
          print 'Warning: assignFFCharge: skipping errors'
    # assign the charges
    if termini == 'neutral':
        modeN = 'Nterminus'
        modeC = 'Cterminus'
    elif termini == 'charged':
        modeN = 'Nterminus+'
        modeC = 'Cterminus-'
    elif termini == 'none':
        modeN = 'body'
        modeC = 'body'
    else:
        print >> sys.stderr, "Error: termini = '%s' but only: neutral, charged and none allowed" % termini
        sys.exit(1)   
    rl = makeList(al)
    for r in rl:
        if r.rName in name321:
            # find out if the residue is on terminal N (first) or terminal C (second)
            # if it is both on terminal C and N, treat it as terminal N
            a = r.atomWithName('N')        
            if a == -1:
                sys.exit("Error: residue %s misses atom N" % r)
            isTerminalN = True
            for b in a.bonds:
                if b.aName == 'C':
                    isTerminalN = False
                    break
            if isTerminalN:
                r.assignFFCharge(store,position=modeN)
            else:
                a = r.atomWithName('C')        
                if a == -1:
                    sys.exit("Error: residue %s misses atom C" % r)
                isTerminalC = True
                for b in a.bonds:
                    if b.aName == 'N':
                        isTerminalC = False
                        break
                if isTerminalC:
                    r.assignFFCharge(store,position=modeC)
                else:
                    r.assignFFCharge(store,position='body')
        elif r.rName in knownResidues:
            r.assignFFCharge(store)   
        else:
            print 'Warning: skipping charge and ffType assignment for residue %s' % r

def findCis(rl, name = 'residue list'):
    '''prints a warning if residues rl[i] and rl[i+1] have a cis bond
    only checked for amino acids
    goes throught the whole residue list'''
    for i in range(len(rl)-1):
        r1 = rl[i]
        r2 = rl[i+1]
        if r1.rName in name321 and r2.rName in name321:
            if r1.atomWithName('C').isBonded(r2.atomWithName('N')):
                if r1.isCis(r2):
                    print 'Warning: %s: Residues %s and %s make cis bond.' % (name,r1,r2)

def findCisMany(al,xyzlist = [], name = 'atom list'):
    '''checks if the protein in the atom list al contains cis bonds
    if so, prints a warning
    optionally: runs such query for all sets of xyz coordinates in the xyzlist
      in this case, the name should be list with the individual frame names
      NOTE: in this case the coordinates left in the al list are from the last frame'''    
    rl = makeList(al)
    if len(xyzlist) == 0:
        findCis(rl,name = name)
    else:
        if len(name) != len(xyzlist):
            sys.exit('Error: the variable name in findCisMany() should provide the list of names for xyz frames')
        else:
            for i in range(len(xyzlist)):
                atoms.putXyzIntoList(al,xyzlist[i])
                findCis(rl,name='frame '+name[i])

def sequence(rl):
    seq = ''
    for r in rl:
        seq += r.oneLetter()
    return seq

def alignResidue(rbase,rtoalign,fraction = 1.0):
    rbaselong = copy.deepcopy(rtoalign)
    backbone = ['N','CA','C','O']
    for atomName in backbone:
        rbaselong.atomWithName( atomName ).xyz = rbase.atomWithName( atomName ).xyz.copy()
    selection = []
    for i,a in enumerate(rbaselong.atoms):
        if a.aName in backbone:
            selection.append(i)
    xyz = atoms.getXyzFromList(rtoalign.atoms)
    xyzold = atoms.getXyzFromList(rbaselong.atoms)
    xyz2, RMSD, RMSDnonsel, RMSDall  = geometry.align(xyzold, xyz, selection = selection)
    xyz3 = fraction * xyz2 + (1.0 - fraction) * xyz
    atoms.putXyzIntoList(rtoalign.atoms, xyz3)

    rtoalign.rNo = rbase.rNo
    rtoalign.chain = rbase.chain
    rtoalign.updateAtoms()
    
    return rtoalign

def extendHelixC(al, endrNo, dropResidues = 1, residuesToAlign = 4, warningLimit = 0.6):
    """see doc of extendHelix()"""

    rla = makeList(al)
    if dropResidues > 0:
        rNomax = rla[-1].rNo
        rla = rla[:(-dropResidues)]
        print '    -> the chain ends with residue %d, but to keep nice alpha helix I am dropping the last %d residues' % (rNomax, dropResidues)
    rNomax = rla[-1].rNo
    print '    -> will extend the chain at C terminus from residue %d to residue %d' % (rNomax, endrNo)

    canonicalHelixFile = os.getenv('VASEK') + '/amino/canonicalHelix.pdb'
    cl = pdb.read( canonicalHelixFile )
    rlc = makeList(cl)
    xyzd = atoms.getXyzFromList(cl).copy()  
    backbone = ['N','CA','C','O']
    for i in range(residuesToAlign):
        irl = -residuesToAlign+i
        for atomName in backbone:
            rlc[i].atomWithName( atomName ).xyz= rla[irl].atomWithName( atomName ).xyz.copy()
    xyzc = atoms.getXyzFromList(cl)
    selection = [i for i in range(len(cl)) if (cl[i].rNo - cl[0].rNo < residuesToAlign) and (cl[i].aName in backbone) ]
    xyzd2, RMSD, RMSDnonsel, RMSDall  = geometry.align( xyzc, xyzd, selection = selection)

    # compute alignment error as the mismatch in position of the last C atom and the aligned one
    atoms.putXyzIntoList(cl, xyzd2)     
    fitError = numpy.linalg.norm( rla[-1].atomWithName('C').xyz - rlc[residuesToAlign-1].atomWithName('C').xyz)
    if fitError > warningLimit:
        print '    WARNING the alignment of the helix extension is not good: the last C atom is shifted by %f (fit RMSD %f)'% (fitError,RMSD)
    else: 
        print '    -> note: the mismatch of the alignmnet (on last C)  is %f (fit RMSD %f)'% (fitError,RMSD)

    for i in range(endrNo - rNomax):
        r =rlc[residuesToAlign+i] 
        rla.append( r )
        r.chain = rla[-2].chain
        r.rNo = rNomax + 1 + i
        r.updateAtoms()

    return makeAtomList(rla), fitError > warningLimit

def extendHelixN(al, endrNo, dropResidues = 1, residuesToAlign = 4, warningLimit = 0.6):
    """see doc of extendHelix()"""

    rla = makeList(al)
    if dropResidues > 0:
        rNomax = rla[0].rNo
        rla = rla[dropResidues:]
        print '    -> the chain ends with residue %d, but to keep nice alpha helix I am dropping the first %d residues' % (rNomax, dropResidues)
    rNomax = rla[0].rNo
    print '    -> will extend the chain at N terminus from residue %d to residue %d' % (rNomax, endrNo)

    canonicalHelixFile = os.getenv('VASEK') + '/amino/canonicalHelix.pdb'
    cl = pdb.read( canonicalHelixFile )
    rlc = makeList(cl)
    xyzd = atoms.getXyzFromList(cl).copy()  
    backbone = ['N','CA','C','O']
    for i in range(residuesToAlign):
        irl = len(rlc)-residuesToAlign+i
        for atomName in backbone:
            rlc[irl].atomWithName( atomName ).xyz = rla[i].atomWithName( atomName ).xyz.copy()
    xyzc = atoms.getXyzFromList(cl)
    selection = [i for i in range(len(cl)) if (cl[-1].rNo - cl[i].rNo < residuesToAlign) and (cl[i].aName in backbone) ]
    xyzd2, RMSD, RMSDnonsel, RMSDall  = geometry.align( xyzc, xyzd, selection = selection)

    # compute alignment error as the mismatch in position of the last C atom and the aligned one
    atoms.putXyzIntoList(cl, xyzd2)     
    fitError = numpy.linalg.norm( rla[0].atomWithName('N').xyz - rlc[-(residuesToAlign)].atomWithName('N').xyz)
    if fitError > warningLimit:
        print '    WARNING the alignment of the helix extension is not good: the last N atom is shifted by %f (fit RMSD %f)'% (fitError,RMSD)
    else: 
        print '    -> note: the mismatch of the alignmnet (on last N)  is %f (fit RMSD %f)'% (fitError,RMSD)

    for i in range(-endrNo + rNomax):
        r =rlc[-(residuesToAlign+1+i)] 
        rla.insert(0, r)
        r.chain = rla[-1].chain
        r.rNo = rNomax - 1 - i
        r.updateAtoms()

    return makeAtomList(rla), fitError > warningLimit

def extendHelix(al, startrNo, endrNo, NdropResidues = 1, NresiduesToAlign = 4, CdropResidues = 1, CresiduesToAlign = 4, warningLimit = 0.6):
    """extends helices by aligning canonical alpha helix to the last two turns (8residues)
    - by default removes the last residue from the source alpha helix
    - works with pdb files -- ie. there is no connectivity
    - counting of the residues relies on the canonicalHelix.pdb to have residues numbered from 1
    - second returned arguement announces if there was a bad fit"""
    chains = set([])
    for a in al:
        chains.add(a.chain)
    chains = list(chains)
    if len(chains) != 1:
        sys.exit('Error in extendHelix: the provided atom list should have exactly one chain, but now has %s' % str(chains))
    print 'CHAIN %s: Modyfing helix to the length: %d-%d (current length %d-%d)' % (chains[0],startrNo,endrNo,al[0].rNo,al[-1].rNo)

    # first throw away extra residues
    al = [a for a in al if startrNo <= a.rNo <= endrNo]
    # now extend residues
    warn1 = warn2 = False
    if al[-1].rNo < endrNo:
        for i in range(0,4):
            print '-> C term try #%d ' % (i+1)
            bl,warn1=extendHelixC(al,endrNo,dropResidues=CdropResidues+i, residuesToAlign=CresiduesToAlign, warningLimit=warningLimit)
            if not warn1:
                break
        al = bl
    if al[0].rNo > startrNo:
        for i in range(0,4):
            print '-> N term try #%d ' % (i+1)
            bl,warn2=extendHelixN(al,startrNo,dropResidues=NdropResidues+i,residuesToAlign=NresiduesToAlign,warningLimit=warningLimit)
            if not warn2:
                break
        al = bl
    return al, (warn1 or warn2)
 
def writePdbWithTER(filename,AL):
    '''writes pdb file, so that protein mainchain atoms have TER when they are not continuous (measured by distance)
    SSBOND lines are printed
    CONECT lines are not printed'''
    lines = []
    getlogin = os.getenv('LOGNAME')
    getenv = os.getenv('HOSTNAME')
    gettime = time.strftime('%X %x %Z')
    lines.append('REMARK saved by %s@%s on %s\n' %(getlogin,getenv,gettime))
    atoms.renumberAtomsFrom(AL)
    rl = makeList(AL)

    print 'Hello world'
    rNo2index = {'%s%d'%(rl[i].chain,rl[i].rNo):(i+1) for i in range(len(rl))} # amber changes residue numbers

    # look at SSBONDS
    ssbonds = []
    for r in rl:
        if r.rName in ['CYS','CYX']:
            sg = r.atomWithName('SG')
            if sg != -1:
                for a in sg.bonds:
                    if a.aName == 'SG':
                        sg2 = a
                        if [sg,sg2] not in ssbonds and [sg2,sg] not in ssbonds:
                            ssbonds.append([sg,sg2])
                        break
    i = 0
    for ssbond in ssbonds:
        i += 1
        s1 = ssbond[0]
        s2 = ssbond[1]
        line = 'SSBOND %3d %s %s %4d    %s %s %4d                                           \n' % (i,s1.rName,s1.chain,s1.rNo,s2.rName,s2.chain,s2.rNo)
        lines.append(line)

    # first residue
    for a in rl[0].atoms:
        lines.append(pdb.returnLine(a))
    # the rest of residues
    for i in range(1,len(rl)):
        r1 = rl[i-1] 
        r2 = rl[i]
        if r1.rName in name321 and r2.rName in name321:
            if not r2.NisCloseToC(r1):
                lines.append("TER   \n")    
        for a in r2.atoms:
            lines.append(pdb.returnLine(a))
    lines.append("TER   \n")    
    lines.append("END   \n")    

    fout = open(filename,'w')
    fout.writelines(lines)
    fout.close()


