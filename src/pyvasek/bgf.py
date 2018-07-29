""" reading and outputting bgf file

BGF formatting:

BIOGRF  332
DESCRP thing
FORCEFIELD DREIDING
REMARK
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,f10.5)
ATOM      21  C    ILE B    5  -18.90354  13.45419   0.97836 C_R    3 0  0.51000 0   0
ATOM      10  N    ARG E  226   94.46000  41.33500  27.46800        0 0 -0.47000
ATOM       7 HA2   GLY E  225   93.74400  41.75700  24.97700        0 0  0.09000
HETATM     6 H6    RES X  999    0.63971  -1.43177   4.88943 H_     1 0 -0.10729
01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
          1         2         3         4         5         6         7         8
 ATOM or HETATM lines 
     atom number -- has to be sorted
            5 letter atom name, has to be unique in a ligand, same symetric hydrogens are ok in protein  
                   3 letter residue name, RES or XXX typically used for ligand
                       1 letter chain, can be number
                          up to 5 digit residue number, has to be sorted
                                coordinates
                                                             5 letter DREIDING atom type
                                                                    number of bond - double bond count as one
                                                                      number of lone pairs
                                                                         cahrge
                                                                                 optional: 1 fixed, 0 movable
                                                                                     not used, probably for solvable area 
FORMAT CONECT (a6,12i6)
CONECT     3     2     4     5
ORDER      3     1     2     1                    note: connect lines should be sorted 
CONECT     4     3                                      numbers within line should be sorted too
ORDER      4     2       
END


PDB formatting:
     
        1          2          3            4      5        6          7
01234567890123456789012345678901234567890123456789012345678901234567890123
ATOM      1  N   GLY   225      92.103  42.564  25.823  1.00  0.00      E1
"""

import os,sys,time,copy,numpy,math
import atoms, residues

# population functions for one atom
def readLine(line):
    """store bgf line info into atom and return the atom"""
    a = atoms.Atom()

    if line[0:6] != 'ATOM  ' and line[0:6] != 'HETATM':
        print >> sys.stderr, 'This line is not ATOM/HETATM line: %s' % line
        sys.exit(1)
    a.aTag  = line[0:6]
    a.aNo   = int(line[6:12])
    a.aName = line[13:18].strip()
    a.rName = line[19:23].strip()
    a.chain = line[23:24]
    try:
        a.rNo   = int(line[24:30]) # vmd saves the residue number until characted 30
    except ValueError: 
        if line[24:29].strip() == '':
            a.rNo = 0
        else:
            a.rNo   = int(line[24:29]) # but normally character 29 can be a letter
    a.xyz   = numpy.array( [float(line[30:40]), float(line[40:50]), float(line[50:60]) ])
    a.ffType= line[61:66] 
    a.lpair = int(line[69:71])
    a.charge= float(line[72:80])

    try:                fixed = int(line[80:82])
    except ValueError:  fixed = 0
    if fixed == 1:      a.fixed = True    
    elif fixed == 0:    a.fixed = False
    else:
        print >> sys.stderr, 'Value of fixed is invalid %d, should be 0 or 1. ' % fixed
        sys.exit(1)

    return a

# make output lines:
def returnLine(atom, mode = 'new'):
    """make a bgf line from atom
    mode = 'new' ... the alignment of the atom names is consistent with amber and pbd database
    mode = 'old' ... the alignment of the atom names is consistent with scream and old biogroup stuff
    """
    if atom.fixed: fixed = 1
    else:          fixed = 0 
    if mode == 'old':
        if atom.aTag == 'ATOM  ' and len(atom.aName) <= 4 and atom.aName[0] != 'H':
            aName = ' '+atom.aName
        else:            
            aName = atom.aName
    else:
        if atom.aTag == 'ATOM  ' and len(atom.aName) <= 3:
            aName = ' '+atom.aName
        else:            
            aName = atom.aName
    # previos version:      #line = '%6s%6d %4.4s  %3s %1s %4d %10.5f%10.5f%10.5f %-5s  %1d %1d %8.5f%2d%4d\n'
    #                                                                 -note: last char is unused
    #(a6,        1x,i5,   1x, a5,         1x, a3,1x,      a1,         1x,a5,    3f10.5,                                1x, a5,          i3,              i2,        1x, f8.5,       i2,    f8.5)
    #'%-6s       %6d          %-5s            %-4s        %1s         %5d       %10.5f       %10.5f       %10.5f           %-5s         %3d              %2d            %8.5f       %2d    %4d\n'
    #(atom.aTag, atom.aNo,    atom.aName,     atom.rName, atom.chain, atom.rNo, atom.xyz[0], atom.xyz[1], atom.xyz[2],     atom.ffType, len(atom.bonds), atom.lpair,    atom.charge,fixed, 0)
    line = '%-6.6s%6d %-5.5s %-4.4s%1.1s%5d %10.5f%10.5f%10.5f %-5.5s%3d%2d %8.5f%2d%4d\n'
    line %= (atom.aTag, atom.aNo, aName, atom.rName, atom.chain, atom.rNo, atom.xyz[0], \
             atom.xyz[1], atom.xyz[2], atom.ffType, len(atom.bonds), atom.lpair, \
             atom.charge,fixed,0)
    return line

def returnConect(atom):
    """make a bgf connect line from atom, including the order line, if necessary

    returns list of length one or two or an empty string
    """
    if len(atom.bonds) == 0:
        return ''

    # sort the numbers
    atom.sortBondsList()
            
    conect_line = 'CONECT%6d' % atom.aNo
    order_line  = 'ORDER %6d' % atom.aNo
    print_order = False

    for bonded_atom in atom.bonds:
        token = '%6d' % bonded_atom.aNo
        conect_line += token
    for bonded_order in atom.order:
        if int(bonded_order) == bonded_order:
            token = '%6d' % bonded_order
        else:
            token = '%6.1f' % bonded_order
        order_line += token
        if bonded_order != 1 or atom.aTag == 'HETATM':
            print_order = True

    if print_order:
        return [conect_line+'\n', order_line + '\n']
    else:
        return [conect_line + '\n']

# population functions for files
def getOrder(orderString):
    orderString = orderString.strip()
    if orderString == 'am':
        return 1.3
    elif orderString == 'ar':
        return 1.5
    if orderString == '1.3':
        return 1.3
    elif orderString == '1.5':
        return 1.5
    else:
        return int(orderString)

def readLines(lines):
    """reads the bgf which is passes as list of lines, returns list of atoms, discard advanced features of bgf"""
    AL = []
    i = 0
    firstConectLine = 0
    # read the atom lines
    for i in range(len(lines)):
        line = lines[i]
        if line[0:6] == 'CONECT' or line[0:13] == 'FORMAT CONECT':
            firstConectLine = i
            break
        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
            AL.append( readLine(line) )
    # read the conect lines
    # first make a dictionary of the atom numbers to atoms
    adict = {}
    for atom in AL:
        if atom.aNo in adict:
            sys.exit('Error: This atom number appears twce: %d, %s' % (atom.aNo,atom))
        adict[atom.aNo] = atom
    # now can read the conect lines
    i = firstConectLine
    presentBonds = {}  # store here the number of time each bond was read to check correctness of the bgf
    while i<len(lines):
        if lines[i].startswith('END'):
            break
        if not lines[i].startswith('CONECT'):
            i += 1
            continue
        conect = []
        line = lines[i][6:]
        maxnumber = len(line) / 6 # integer division
        for j in range(maxnumber):
            possiblenumber = line[(6*j):(6*j+6)]
            possiblenumber = possiblenumber.strip()
            if possiblenumber != '':
                conect.append(int(possiblenumber))
        if  i+1 < len(lines) and lines[i+1][0:5] == 'ORDER':
            # reading both conect line and order lines
            orderStrings = lines[i+1].split()
            order = [getOrder(orderi) for orderi in orderStrings[1:]]
            i+=1
            if len(conect) != len(order):
                sys.exit('Error: conect and order lines have different length: %s, %s'%(lines[i],lines[i+1]))
        else:
            order = [1] * len(conect)

        atomi = int(conect[0])
        if atomi not in adict:
            sys.exit('Error: atom number %d on conect line is not in the list of atoms' % atomi)
        atom = adict[atomi]
        for j in range(1,len(conect)):
            atomj = int(conect[j])
            if atomj not in adict:
                sys.exit('Error: atom number %d on conect line is not in the list of atoms' % atomi)
            atom2 = adict[atomj]

            thisBond = (min(atomi,atomj),max(atomi,atomj))
            thisorder = order[j]
            if thisBond in presentBonds:
                if presentBonds[thisBond] == 1:
                    firstBondOrder = atom.getBondOrder(atom2)
                    if thisorder != firstBondOrder:
                        sys.exit('Error: bond %s (index %d) - %s (index %d) has two orders %f and %f ' % (atom,atomi, atom2,atomj,thisorder,firstBondOrder))
                    else:
                        presentBonds[thisBond] = 2
                else:
                    sys.exit('Error: too many bonds between atom %s (index %d) and atom %s (index %d)' % (atom,atomi, atom2,atomj))
            else:
                atom.makeBond(atom2,bondOrder=thisorder)
                presentBonds[thisBond] = 1
        i+=1
    for d in presentBonds:
        if presentBonds[d] != 2:
            sys.exit('Error: bond between atoms of indices %d-%d only occured in the bgf file once' % (d[0],d[1]))

    # renumber atoms to be in range 1..numberOfAtoms
    for i in range(len(AL)):
        AL[i].aNo = i+1
    return AL

def writeLines(AL,mode='auto'):
    """returns list of lines corresponding to the atom list
    need to have a special rule for alignment because of SCREAM
    mode = 'new' ... the alignment of the atom names is consistent with amber and pbd database
    mode = 'old' ... the alignment of the atom names is consistent with scream and old biogroup stuff
    mode = 'auto'... 'old' if one of hte first 24 atoms is amino acid and HN, or HCA, otherwise 'new'
    """
    if mode == 'auto':
        mode = 'new'
        for i in range(min( len(AL), 24 )):
            atom = AL[i]
            if atom.aTag == 'ATOM  ' and (atom.rName in residues.name321) and (atom.aName in ['HN','HCA']):
                mode = 'old'
                break

    atoms.renumberAtomsFrom(AL)
    lines = []
    lines.append('BIOGRF 350\n')
    lines.append('DESCRP atoms\n')
    getlogin = os.getenv('LOGNAME')
    getenv = os.getenv('HOSTNAME')
    gettime = time.strftime('%X %x %Z')
    lines.append('REMARK saved by %s@%s on %s\n' %(getlogin,getenv,gettime))
    lines.append('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,f8.5)\n')
    for a in AL: lines.append(returnLine(a,mode=mode))
    lines.append('FORMAT CONECT (a6,14i6)\n')
    lines.append('FORMAT ORDER (a6,i6,13f6.3)\n')
    for a in AL: 
        conectlines = returnConect(a)
        if conectlines != '':
            lines.extend(conectlines)
    lines.append("END\n")    
    return lines

def readFile(file, debug=False):
    AL = readLines(open(file).readlines())
    if debug: print '%d atoms read from %s.' %(len(AL), file)
    return AL

def writeFile(outfile,AL,debug=False,mode='auto'):
    fout = open(outfile,'w')
    fout.writelines(writeLines(AL,mode=mode))
    fout.close()
    if debug: print 'File %s written.' % outfile

def read(file,debug=False):
    '''the main method for reading bgf files'''
    return readFile(file,debug=debug)
    
def write(outfile,al,debug=False,mode='auto'):
    '''the main method for writing bgf files 
    mode = 'new' ... the alignment of the atom names is consistent with amber and pbd database
    mode = 'old' ... the alignment of the atom names is consistent with scream and old biogroup stuff
    '''
    writeFile(outfile,al,debug=debug,mode=mode)

def selectcatbgf(lines, framename):
    '''lines are lines of a catbgf file
    framename is the string in the 'DESCRP' line
    eg. 'c10116.p.c.s.bgf'
    returns atom list
    '''
    thisname = False
    thislines = []
    for line in lines:
        if line.startswith('DESCRP'):
            s = line.split()[1]
            if s == framename:
                thisname = True
        elif line.startswith('END'):
            thislines.append(line) 
            if thisname:            
                return readLines(thislines)
            else:
                thislines = []
        else:
            thislines.append(line)     
    sys.exit('Error: file %s was not present in the catbgf' % framename)

def extractcatbgf(lines):
    '''made to extract individual files from adams catbgf
    the catbgf file does not need to have all frames with the same connectivity
    lines are lines of a catbgf file
    returns list of atoms lists
    '''
    als = []
    thislines = []
    for line in lines:
        if line.startswith('END'):
            thislines.append(line) 
            als.append(readLines(thislines))
            thislines = []
        else:
            thislines.append(line)     
    return als

def catbgfframenames(lines):
    '''lines are lines of a catbgf file
    returns the list of framenames, ie. strings in the 'DESCRP' line
    if there is more DESCRP lines the last one is used
    if there is no DESCRP line the frame name is ''
    returns list of strings
    '''
    framenames = []
    framename = ''
    for line in lines:
        if line.startswith('DESCRP'):
            framename = line.split()[1]
        elif line.startswith('END'):
            framenames.append(framename)
            framename = ''
    return framenames


