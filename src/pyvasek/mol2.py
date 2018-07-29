import numpy, os, time, sys
import atoms

def read(fileName):
    return readFile(fileName)

def write(fileName, al):
    writeFile(fileName, al)

def readFile(fileName):
    with open(fileName,'r') as f:
        line=f.readline()
        while not (line.startswith('@<TRIPOS>MOLECULE') or line == ''):
            line=f.readline()
        line=f.readline()
        line=f.readline()
        s = line.split()
        numAtoms = int(s[0])
        numBonds = int(s[1])
        numRes = int(s[2])
        AL = [atoms.Atom() for i in range(numAtoms)]
        while not (line.startswith('@<TRIPOS>ATOM') or line == ''):
            line=f.readline()
        for i in range(numAtoms):
            line=f.readline()
            readAtomLine(line,AL[i])
        while not (line.startswith('@<TRIPOS>BOND') or line == ''):
            line=f.readline()
        for i in range(numBonds):
            line=f.readline()
            readBondLine(line, AL)
        while not (line.startswith('@<TRIPOS>SUBSTRUCTURE') or line == ''):
            line=f.readline()
        resRoot = []
        resChain = []
        for i in range(numRes):
            line=f.readline()
            s = line.split()
            root = int(s[2])
            chain = s[5]
            resRoot.append(root)
            resChain.append(chain)
        chainIndex = 0
        if len(resChain) > 0:
            for a in AL:
                if len(resRoot) > chainIndex+1:
                    if a.aNo > resRoot[chainIndex+1]:
                        chainIndex = chainIndex + 1 
                chain = resChain[chainIndex]
                a.chain = chain    
    return AL

def readAtomLine(line, a):
    """parameter a is a blank atom"""
    s = line.split()
    a.aNo   = int(s[0])
    a.aName = s[1]
    a.xyz   = numpy.array( [float(s[2]), float(s[3]), float(s[4]) ])
    a.ffType = s[5]
    ls = len(s)
    if ls > 6:
        a.rNo = int(s[6])
    if ls > 7:
        a.rName = s[7]
    if ls > 8:
        a.charge = float(s[8])

def readBondLine(line, AL):
    s = line.split() 
    origin_atom_id = int(s[1]) 
    target_atom_id = int(s[2])
    bond_type = s[3]
    if bond_type == 'nc':
        return
    elif bond_type == 'un' or bond_type == 'du' or bond_type == '1':
        bondorder = 1
    elif bond_type == '2':
        bondorder = 2
    elif bond_type == '3':
        bondorder = 3
    elif bond_type == 'am':
        bondorder = 1.3
    elif bond_type == 'ar':
        bondorder = 1.5
    else:
        sys.exit('Error in mol2.readBondLine: Unknown bond type: %s' % bond_type)
    AL[origin_atom_id-1].makeBond( AL[target_atom_id-1], bondorder)

def writeAtomLine(a):
    # changed to 5 decimal places (as bgf on 1/30/2012)
    # return '%7d %-8s %9.4f %9.4f %9.4f %-5s %5d %3s %13.6f\n'%(a.aNo, a.aName, a.xyz[0], a.xyz[1], a.xyz[2], a.ffType, a.rNo,a.rName,a.charge)
    return '%7d %-8s %15.10f %15.10f %15.10f %-5s %5d %3s %13.6f\n'%(a.aNo, a.aName, a.xyz[0], a.xyz[1], a.xyz[2], a.ffType, a.rNo,a.rName,a.charge)

def writeBondLine(a,i):
    line = ''
    for ib in range(len(a.bonds)):
        b = a.bonds[ib]
        order = a.order[ib]
        if a.aNo < b.aNo: 
            if order == 1.5:
                code = 'ar'
            elif order == 1.3:
                code = 'am'                
            else:
                code = '%d' % order 
            line = line + '%6d %4d %4d %-4s\n' % (i,a.aNo,b.aNo,code) 
            i = i+1
    return i,line

def writeFile(fileName, AL):
    atoms.renumberAtomsFrom(AL)
    with open(fileName,'w') as f:
        getlogin = os.getenv('LOGNAME')
        getenv = os.getenv('HOSTNAME')
        gettime = time.strftime('%X %x %Z')
        f.write('# saved by %s@%s on %s\n' %(getlogin,getenv,gettime))
        f.write('@<TRIPOS>MOLECULE\n')
        f.write('title\n')
        
        atoms.renumberAtomsFrom(AL)
        
        numAtoms = len(AL)
        numBonds = 0
        numRes = 1
        prevRes = 0
        for a in AL:
            numBonds = numBonds + len(a.bonds)
            b = AL[prevRes]
            if b.rNo!=b.rNo or b.chain!=b.chain:
                numRes = numRes + 1
                prevRes = a.aNo - 1
        numBonds = numBonds / 2  

        f.write('%5d %5d %5d     0     0\n' % (numAtoms, numBonds, numRes))
        f.write('SMALL\n')
        f.write('USER_CHARGES\n\n\n')
        f.write('@<TRIPOS>ATOM\n')
        for a in AL:
            f.write(writeAtomLine(a))
        f.write('@<TRIPOS>BOND\n')
        i = 1
        for a in AL:
            i,line = writeBondLine(a,i) 
            f.write(line)
        f.write('@<TRIPOS>SUBSTRUCTURE\n')
        # the last thing is to write residues
        i = 1
        format = '%6d %3s %9d ****              0 %4s  **** \n'
        line = format % (i, AL[0].rName, 1, AL[0].chain)
        f.write(line)
        prevRes = 0
        for a in range(1,numAtoms):
            b = AL[prevRes]
            if b.rNo!=b.rNo or b.chain!=b.chain:
                i = i+1
                line = format % (i, a.rName, ia+1, a.chain)                
                prevRes = a.aNo - 1

                
