""" reading and outputting pdb file

PDB formatting:
     
        1          2          3            4      5        6          7
01234567890123456789012345678901234567890123456789012345678901234567890123
ATOM      1  N   GLY   225      92.103  42.564  25.823  1.00  0.00      E1
"""

import os,sys,time,copy,numpy,math
import atoms

# population functions for one atom
def readLine(line):
    """store bgf line info into atom and return the atom"""
    a = atoms.Atom()

    if line[0:6] != 'ATOM  ' and line[0:6] != 'HETATM':
        print >> sys.stderr, 'This line is not ATOM/HETATM line: %s' % line
        sys.exit(1)
    a.aTag  = line[0:6]
    a.aNo   = int(line[6:11])
    a.aName = line[12:16].strip()
    a.rName = line[17:20].strip()
    a.chain = line[21]
    a.rNo   = int(line[22:26])
    a.xyz   = numpy.array( [float(line[30:38]), float(line[38:46]), float(line[46:54]) ])

    try:                a.occupancy = float(line[54:60])
    except ValueError:  a.occupancy = 1.0
    try:                a.beta  = float(line[60:66])
    except ValueError:  a.beta = 0.0

    a.lpair = 0
    a.charge= 0.0
    a.fixed = False

    try:                a.ffType = line[76:78]
    except ValueError:  a.ffType = ''

    return a

# make output lines:
def returnLine(atom):
    """make a bgf line from atom"""
    atom.aTag = '%-6s'%atom.aTag

    aName = atom.aName.strip()
    if aName[1:2].islower():
        element = aName[0:2]
    else:
        element = aName[0:1]
    if atom.aTag == 'ATOM  ' and len(aName) <= 3:
        aName = ' '+aName
    try:                   occupancy = atom.occupancy
    except AttributeError: occupancy = 1.0
    try:                   beta = atom.beta
    except AttributeError: beta = 0.0

    line = '%s%5d %-4.4s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s\n'
    line %= (atom.aTag, atom.aNo, aName, atom.rName, atom.chain, atom.rNo, atom.xyz[0], \
             atom.xyz[1], atom.xyz[2], occupancy, beta, element)
    return line

def returnConect(atom, al):
    """make a bgf connect line from atom, including the order line, if necessary

    returns list of length one or two or an empty string
    checks if the connected atoms are present
    """
    if len(atom.bonds) == 0:
        return ''
    # sort the numbers
    atom.sortBondsList()
    conect_line = 'CONECT%5d' % atom.aNo
    for bonded_atom in atom.bonds:
        if bonded_atom not in al:
            sys.exit('Error in writing pdb file: atom %s is bonded to atom %s which is not present in the list of atoms for output' % (atom, bonded_atom))
        token = '%5d' % bonded_atom.aNo
        conect_line += token
    return ['%-70s\n' % conect_line]

# population functions for files

def readLines(lines, alternate = 'A'):
    """reads the bgf which is passes as list of lines, returns list of atoms, discard advanced features of bgf"""
    AL = []
    i = 0
    firstConectLine = -1
    # read the atom lines
    for i in range(len(lines)):
        line = lines[i]
        if line[0:6] == 'CONECT':
            firstConectLine = i
            break
        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
            if (line[16] == ' ') or (line[16] == alternate):
                AL.append( readLine(line) )
    # read the conect lines
    if firstConectLine >= 0:
        i = firstConectLine
        # first make a dictionary of the atom numbers to atoms
        dict = {}
        for atom in AL:
            if atom.aNo in dict:
                print >> sys.stderr, 'Warning: This atom number appears twice: %d, %s' % (atom.aNo,atom)
            dict[atom.aNo] = atom
        # now can read the conect lines
        while i<len(lines):
            if lines[i][0:6] == 'CONECT': 
                # reading only conect line
                conect = []
                maxnumber = (len(lines[i])-6) / 5 # integer division 
                for j in range(maxnumber):
                    possiblenumber = lines[i][(6+5*j):(6+5*(j+1))]
                    possiblenumber = possiblenumber.strip()
                    if possiblenumber != '':
                        conect.append(int(possiblenumber))
                if int(conect[0]) in dict:
                    # we only want the connect lines which occured for the chosen alternative
                    atom = dict[int(conect[0])]
                    for j in range(1,len(conect)):
                        try:
                            atom2 = dict[int(conect[j])]
                            if not atom.isBonded(atom2): 
                                atom.makeBond(atom2)
                        except KeyError:
                            print 'Warning: atom %s number %d has connect line entry bond to atom number %d, but this atom is not present' % (atom,int(conect[0]),int(conect[j]))
            i+=1
    # renumber atoms to be in range 1..numberOfAtoms
    for i in range(len(AL)):
        AL[i].aNo = i+1
    return AL

def writeLines(AL,withConect=True):
    """returns list of lines corresponding to the atom list"""
    lines = []
    getlogin = os.getenv('LOGNAME')
    getenv = os.getenv('HOSTNAME')
    gettime = time.strftime('%X %x %Z')
    lines.append('REMARK saved by %s@%s on %s\n' %(getlogin,getenv,gettime))
    atoms.renumberAtomsFrom(AL)
    for a in AL: lines.append(returnLine(a))
    if withConect:
        for a in AL: 
            conectlines = returnConect(a, AL)
            if conectlines != '':
                lines.extend(conectlines)
    lines.append("END\n")    
    return lines

def readFile(file, debug=False, alternate = 'A'):
    AL = readLines(open(file).readlines(), alternate = alternate)
    if debug: print '%d atoms read from %s.' %(len(AL), file)
    return AL

def writeFile(outfile,AL,withConect=True,debug=False):
    fout = open(outfile,'w')
    fout.writelines(writeLines(AL,withConect=withConect))
    fout.close()
    if debug: print 'File %s written.' % outfile
  
def read(file, alternate = 'A'):
    return readFile(file, alternate = alternate)
    
def write(outfile,al,withConect=True):
    writeFile(outfile,al,withConect=withConect)


pdbDefinition = {'HELIX ': ''' 1 -  6     Record name      "HELIX "
 8 - 10     Integer          serNum       Serial number of the helix.
                                          This starts at 1 and increases
                                          incrementally.
12 - 14     LString(3)       helixID      Helix identifier. In addition
                                          to a serial number, each helix is
                                          given an alphanumeric character
                                          helix identifier.
16 - 18     Residue name     initResName  Name of the initial residue.
20          Character        initChainID  Chain identifier for the chain
                                          containing this helix.
22 - 25     Integer          initSeqNum   Sequence number of the initial
                                          residue.
26          AChar            initICode    Insertion code of the initial
                                          residue.
28 - 30     Residue name     endResName   Name of the terminal residue of
                                          the helix.
32          Character        endChainID   Chain identifier for the chain
                                          containing this helix.
34 - 37     Integer          endSeqNum    Sequence number of the terminal
                                          residue.
38          AChar            endICode     Insertion code of the terminal
                                          residue.
39 - 40     Integer          helixClass   Helix class (see below).
41 - 70     String           comment      Comment about this helix.
72 - 76     Integer          length       Length of this helix.''',

'SHEET ':''' 1 -  6     Record name      "SHEET "
 8 - 10     Integer          strand       Strand number which starts at 1 
                                          for each strand within a sheet 
                                          and increases by one.
12 - 14     LString(3)       sheetID      Sheet identifier.
15 - 16     Integer          numStrands   Number of strands in sheet.
18 - 20     Residue name     initResName  Residue name of initial residue.
22          Character        initChainID  Chain identifier of initial 
                                          residue in strand.
23 - 26     Integer          initSeqNum   Sequence number of initial 
                                          residue in strand.
27          AChar            initICode    Insertion code of initial residue
                                          in strand.
29 - 31     Residue name     endResName   Residue name of terminal residue.
33          Character        endChainID   Chain identifier of terminal
                                          residue.
34 - 37     Integer          endSeqNum    Sequence number of terminal
                                          residue.
38          AChar            endICode     Insertion code of terminal 
                                          residue.
39 - 40     Integer          sense        Sense of strand with respect to
                                          previous strand in the sheet. 0
                                          if first strand, 1 if parallel,
                                          -1 if anti-parallel.
42 - 45     Atom             curAtom      Registration. Atom name in 
                                          current strand.
46 - 48     Residue name     curResName   Registration. Residue name in
                                          current strand.
50          Character        curChainId   Registration. Chain identifier in
                                          current strand.
51 - 54     Integer          curResSeq    Registration. Residue sequence
                                          number in current strand.
55          AChar            curICode     Registration. Insertion code in
                                          current strand.
57 - 60     Atom             prevAtom     Registration. Atom name in
                                          previous strand.
61 - 63     Residue name     prevResName  Registration. Residue name in
                                          previous strand.
65          Character        prevChainId  Registration. Chain identifier in
                                          previous strand.
66 - 69     Integer          prevResSeq   Registration. Residue sequence
                                          number in previous strand.
70          AChar            prevICode    Registration. Insertion code in
                                              previous strand.''',
'SSBOND':''' 1 -  6     Record name      "SSBOND"
 8 - 10     Integer          serNum       Serial number.
12 - 14     LString(3)       "CYS"        Residue name.
16          Character        chainID1     Chain identifier.
18 - 21     Integer          seqNum1      Residue sequence number.
22          AChar            icode1       Insertion code.
26 - 28     LString(3)       "CYS"        Residue name.
30          Character        chainID2     Chain identifier.
32 - 35     Integer          seqNum2      Residue sequence number.
36          AChar            icode2       Insertion code.
60 - 65     SymOP            sym1         Symmetry oper for 1st resid
67 - 72     SymOP            sym2         Symmetry oper for 2nd resid'''}

def parsePdbDefinition():
    dic = {}
    for item in pdbDefinition:
        lines = pdbDefinition[item].split('\n')
        itemdic = {}            
        for line in lines:
            if line[0:2].strip().isdigit():
                ifrom = int(line[0:2].strip()) - 1 
                if line[5:7].strip().isdigit():
                    ito = int(line[5:7].strip())
                else:
                    ito = ifrom + 1
                name = line[29:40].strip()
                itemdic[name] = (ifrom,ito)
        dic[item] = itemdic
    return dic

def parsePdbHeaderLine(line, dic):
    out = {}
    for key in dic:
        ifrom,ito = dic[key]
        out[key] = line[ifrom:ito]
    return out

def parsePdbHeader(lines, whichEntry):
    '''which entry is pdb line name, 'HELIX ', 'SHEET ', SSBOND and has to be defined in my pdbDefinition dictionary'''

    dic = parsePdbDefinition()[whichEntry]
    parsed = []
    for line in lines:
        if line[0:6] == whichEntry:
            parsedLine = parsePdbHeaderLine(line, dic)
            parsed.append(parsedLine)
    return parsed






