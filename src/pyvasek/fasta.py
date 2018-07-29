""" reading and writing fasta files
"""

import os,sys

def readLinesOrKeep(fastafile):
    """filename can be a file name with fasta file or list of lines with a fasta file"""
    if isinstance(fastafile,list):
        return fastafile
    else:
        return open(fastafile).readlines()


def read(fastafile):
    lines = readLinesOrKeep(fastafile)
    seq = ''
    for line in lines:
        if line[0].isalpha():
            seq += line.strip() 

    return seq

def readMultiple(fastafile):
    '''read fasta with multiple sequences'''
    lines = readLinesOrKeep(fastafile)
    seqs = {}
    id = ''
    seq = ''
    for line in lines:
        if line[0].isalpha() or line[0] == '-':
            seq += line.strip()
        else:
            if seq != '':
                assert id != '', 'error in fasta file on line: %s'% line
                seqs[id] = seq
            seq = ''
            id = ''
            if line[0] == '>':
                id = line.split('|')[1]
    if seq != '':
        assert id != '', 'error in fasta file on line, there is no id'
        seqs[id] = seq
    return seqs

def getName(fastafile):
    lines = readLinesOrKeep(fastafile)

    line = lines[0]
    s = line.split()
    gn = s[-3].lower()
    if gn.startswith('gn='):
        name = gn[3:]
    else:
        name = '?'
    return name

        
def getUniprot(fastafile):
    lines = readLinesOrKeep(fastafile)

    line = lines[0]
    s = line.split('|')
    return s[1]

def readFormattedSecondary(formattedfile):
    lines = readLinesOrKeep(formattedfile)
    keep = {}
    prevkey = ''
    for line in lines:
        if line.startswith('KEY:'):
            break
        elif line.startswith('        ') or (len(line.strip()) == 0):
            continue
        s = line.split(':')
        key = s[0].strip().lower()
        if key == prevkey:
            key += '_conf'
        prevkey = key

        text = s[1].strip()
        if key in keep:
            keep[key] += text
        else:
            keep[key] = text
    return keep    

def readRealTM(residueList):
    """as arguement takes residue list from 7 TM bundles (without loops)
    and returns two arrays tmstart and tmend with residues numbers
    """
    rl = residueList
    tmstart = [ rl[0].rNo ]
    tmend = []
    for i in range(1,len(rl)):
        if rl[i].chain != rl[i-1].chain:
            tmstart.append( rl[i].rNo )
            tmend.append( rl[i-1].rNo )
    tmend.append( rl[-1].rNo )

    return tmstart,tmend

def readHeader(fastafile):
    lines = readLinesOrKeep(fastafile)
    line = lines[0]
    if line[0] == '>':
        return line
    else:
        return '>protein'

def writeMFTA(outfile, oldmfta, sequence, structure = '', tmstart = [], tmend =[]):
    """updates mfta file with the new sequence and TM regions
    0. either structure string or tmstart and tmend (arrays with 7 integers each) have to be provided
    1. letter H in the structure string denotes the tm region
    2. from the oldmfta it read the last three colums of the * tm lines
       and header 
    3. sequence can contain '-' characters in which case the letter is skiiped
    """
    lines = readLinesOrKeep(oldmfta)
    header = readHeader(lines)
    if structure != '':
        tmstart = []
        tmend = []
        tmcount = 0
        intm = False
        j = 0
        if len(sequence) != len(structure):
            sys.exit('Error in writeMFTA: length of sequnce has to be the same as length of structure')
        for i in range(len(sequence)):
            if sequence[i] == '-':
                pass
            else:
                j += 1
                if intm:
                    if structure[i] == 'H':
                        pass
                    else:
                        intm = False
                        tmend.append(j-1)
                else:
                    if structure[i] == 'H':
                        intm = True
                        tmstart.append(j)
                    else:
                        pass

    tm = []
    tmcount = 0
    for line in lines:
        if (line[0:1]=='*') and ( line.find('tm') >= 0):
            s = line.split()
            tmcount += 1
            tm.append( '*  %dtm %4d %4d    X    X    X    X    X    X    %s %4s    %s\n' % (tmcount,tmstart[tmcount-1],tmend[tmcount-1],s[-3],s[-2],s[-1]) )

    with open(outfile,'w') as f:
        # header
        f.write(header)
        # seq
        seq = []
        for s in sequence:
            if s != '-':
                seq.append(s)
        sequence = ''.join( seq )
        n = len(sequence)
        d = n / 60
        rem = n % 60
        for i in range(d):
            f.write(sequence[(60*i):(60*(i+1))])
            f.write('\n')        
        if rem > 0:
            f.write(sequence[(60*d):])
            f.write('\n')        
        # tm
        for itm in tm:
            f.write(itm)
        # hcp
        for i in range(1,8):
            f.write('*  %dhpc     X.00     X.00     X.00     X.00     X.00     X.00\n'%i)

def writeMFTA_onlyMT(outfile, infasta, tmstart, tmend, eta = []):
    """makes mfta file with TM entries only from a fasta file and starts and ends of TM regions
    """
    lines = readLinesOrKeep(infasta)
    seq = read(lines)

    with open(outfile,'w') as f:
        for line in lines:
            if len(line.strip()) > 0 and line[0] != '*':
                f.write(line)
        # tm
        for i in range(7):
            tmcount = i + 1
            if len(eta) == 0:
                f.write('*  %dtm %4d %4d    X    X    X    X    X    X    ?    ?    A\n' % (tmcount,tmstart[tmcount-1],tmend[tmcount-1]) )
            else:
                letter = seq[eta[tmcount-1]]
                f.write('*  %dtm %4d %4d    X    X    X    X    X    X    %s %4d    A\n' % (tmcount,tmstart[tmcount-1],tmend[tmcount-1], letter, eta[tmcount-1]) )

def writeMFTA_fromSeq(outfile, seq, header = '>protein', tmstart=[], tmend=[], eta=[]):
    with open(outfile,'w') as f:
        f.write(header.strip())
        f.write('\n')
        n = len(seq)
        i = 0
        while i < n:
            f.write(seq[i:(i+60)])
            f.write('\n')
            i+=60
        # tm
        if tmstart != [] and tmend != []:
          for i in range(7):
            tmcount = i + 1
            if len(eta) == 0:
                f.write('*  %dtm %4d %4d    X    X    X    X    X    X    ?    ?    A\n' % (tmcount,tmstart[tmcount-1],tmend[tmcount-1]) )
            else:
                letter = seq[eta[tmcount-1]]
                f.write('*  %dtm %4d %4d    X    X    X    X    X    X    %s %4d    A\n' % (tmcount,tmstart[tmcount-1],tmend[tmcount-1], letter, eta[tmcount-1]) )

def readMFTATM(mftafile, column = 1):
    """reads the pair of columns of the tm lengths from the mfta file
    default is the first column, or other columns can be requested"""
    lines = readLinesOrKeep(mftafile)
    tmstart = []
    tmend = []
    tmcount = 1
    for line in lines:
        if line[0:1] == '*':
            s=line.split()
            if s[1] == '%dtm'%(tmcount):
                tmcount += 1
                tmstart.append(int(s[2*column]))
                tmend.append(int(s[2*column+1]))
    if len(tmstart) == 7 == len(tmend):
        return tmstart,tmend
    else:
        sys.exit('Error in readMFTA did not manage to obtain 7 tm lengths')

def readMFTAeta(mftafile):
    """reads the eta residue from the mfta file 12th (second to last column) """
    lines = readLinesOrKeep(mftafile)
    tmeta = []
    tmcount = 1
    for line in lines:
        if line[0:1] == '*':
            s=line.split()
            if s[1] == '%dtm'%(tmcount):
                tmcount += 1
                tmeta.append(int(s[11]))
    if len(tmeta) == 7:
        return tmeta
    else:
        sys.exit('Error in readMFTA did not manage to obtain 7 tm eta residues')

def blosum62():
    '''returns dictionary with the blosum62 matrix
       gap represented as '-' '''
    #  Matrix made by matblas from blosum62.iij
    #  * column uses minimum score
    #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
    #  Blocks Database = /data/blocks_5.0/blocks.dat
    #  Cluster Percentage: >= 62
    #  Entropy =   0.6979, Expected =  -0.5209
    blos='''A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  -
        A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
        R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
        N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
        D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
        C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
        Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
        E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
        G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
        H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
        I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
        L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
        K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
        M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
        F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
        P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
        S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
        T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
        W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
        Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
        V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
        B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
        Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
        X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
        - -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1'''
    result = {}
    lines = blos.split('\n')
    lines = [line.strip() for line in lines]
    keys = lines[0].split()
    for line in lines[1:]:
        items = line.split()
        letter = items[0]
        result[letter] = {}
        for i in range(len(keys)):
            result[ letter ][ keys[i] ] = int( items[i+1] )
    return result

def sequenceSimilarity(seq1,seq2):
    ''' gaps in the alignment have to be '-' 
    returns: number of identical residues, number of similar residues, number of residues to compare, total lengths of strings (as floats)
    the two string have to have the same length

    residues are similar if the blosum entry is > 0
    '''
    n = len(seq1)
    if n != len(seq2):
        sys.exit('Error in sequenceSimilarity(): the sequences have different lengths %d and %d' % (n,len(seq2)))
    blos = blosum62()

    idencount = 0
    simcount = 0
    length = 0
    for i in range(n):
        l1 = seq1[i]
        l2 = seq2[i]
        if l1 != '-' or l2 != '-':
            if blos[l1][l2] > 0:
                simcount += 1
            if l1==l2:
                idencount += 1
            length += 1
    return float(idencount), float(simcount), float(length), float(n)

def resno2alnno(seq, resno):
    '''alignment sequence has gaps: -, but we want to know where 
    is some residue number in the alignment
    let's index alignments from 1 as sequences'''
    ai = 0
    ri = 0
    for i in range(len(seq)):
        if seq[i] == '-':
            ai += 1
        else:
            ai += 1
            ri += 1
            if ri == resno:
                return ai
    sys.exit('error in resno2alnno: residue number %s is not in the sequence %s' %(resno,seq))

def alnno2resno(seq, alnno):
    '''alignment sequence has gaps: -, but we want to know what residue number
    is at some alignments position 
    let's index alignments from 1 as sequences'''
    ai = 0
    ri = 0
    for i in range(len(seq)):
        if seq[i] == '-':
            ai += 1
            if ai == alnno:
                sys.exit('error in alnno2resno: alignment position %d corresponds to a gap, so I cannot return residue number')
        else:
            ai += 1
            ri += 1
            if ai == alnno:
                return ri
    sys.exit('error in alnno2resno: alignment item %s is not in the sequence %s' %(resno,seq))

class Gpcr7tm:
    '''class for storing tm start, stop and eta residue numbers
    hydrophobic center of each helix is taken to be residues with Balesteros no=50
    except for helix 3 where hydrophobic center is central D
    rhod tm 3: central D is A117, bal50 is R135 (so the offset is 117-135=-18)

    TM1: G N xxV          N is 50
    TM2: Lxxx D           D is 50
    TM3: E(D) R Y         R is 50
    TM4: W                W is 50
    TM5: Fxx P            P is 50
    TM6: CWx P FF         P is 50
    TM7: N P xxY          P is 50
    x: variable residue 
    '''
    class Helix:
        def __init__(self, itm):
            self.itm = itm   # 1 through 7
            self.start = -1  # residue number 
            self.stop = -1   # residue number
            self.bal50 = -1

        def getHydrophobicCenter(self):
            if self.itm == 3:
                return self.bal50 - 18
            else:
                return self.bal50

        def setBal50fromHydrophobicCenter(self, rNo):
            if self.itm == 3:
                self.bal50 = rNo + 18
            else:
                self.bal50 = rNo

        def getBal50(self, rNo):
            return 50 + rNo - self.bal50

        def getRNo(self, bal50):
            return self.bal50 + bal50 - 50

        def setHelixFromEta(self, start, stop, eta):
            self.start = start
            self.stop = stop
            self.setBal50fromHydrophobicCenter(eta)

    def __init__(self):
        self.tms = {}
        for i in range(1,8):
            self.tms[i] = self.Helix(i)

    def __init__(self, mftafile):
        self.tms = {}
        for i in range(1,8):
            self.tms[i] = self.Helix(i)
        self.setHelicesFromMfta(mftafile)
    
    def __getitem__(self, index):
        return self.tms[index]

    def setHelicesFromEta(self, tmstart, tmend, tmeta):
        for i in range(7):
            self.tms[i+1].setHelixFromEta(tmstart[i], tmend[i], tmeta[i])
    
    def setHelicesFromMfta(self, mftafile):
        tmstart,tmend = readMFTATM(mftafile)
        tmeta = readMFTAeta(mftafile)
        self.setHelicesFromEta(tmstart, tmend, tmeta)

    def bal2rno(self, balles):
        '''convert Ballesteros indexed number to residue number
        can be N1.50 or 1.50 or N/D1.50
        '''
        try :
            s = balles.split('.')
            if len(s) != 2:
                sys.exit('Error in Gpcr7tm.bal2rno: this is not a valid Ballesteros indexed residue: %s'%balles)
            tm = int( s[0][-1] )
            balrno = int( s[1] )
        except (ValueError, IndexError):
            sys.exit('Error in Gpcr7tm.bal2rno: this is not a valid Ballesteros indexed residue: %s'%balles)

        return self.tms[tm].bal50 - 50 + balrno  

    #def rno2bal(self, balles):
    #    '''this is only unique ofr residues in TM; for loops this is not unique'''
    #    ????

    def readcsv(self, csvfile):
        '''the file has one line header and then 7 lines, eg.:
        tm,ballesterosreference,res,tmstart,tmstop
        1,N1.50,55,35,63
        '''
        lines = open(csvfile).readlines()
        for line in lines[1:8]:
            s = line.split(',')
            tm = int(s[0])
            tmbalno = int(s[1].split('.')[-1])
            tmbal = int(s[2])
            self.tms[tm].bal50 = tmbal - tmbalno + 50 
            self.tms[tm].start = int(s[3])
            self.tms[tm].stop = int(s[4])

    def writecsv(self, csvfile):
        with open(csvfile,'w') as f:
            f.write('tm,ballesterosreference,res,tmstart,tmstop\n')
            for i in range(1,8):
                f.write('%d,%d.50,%d,%d,%d\n' % (i, i, self.tms[i].bal50, self.tms[i].start, self.tms[i].stop) )





