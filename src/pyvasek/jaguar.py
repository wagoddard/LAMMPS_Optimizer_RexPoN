""" writing jaguar input .in and reading details from the output .out """

import subprocess, numpy, os, sys, tempfile, shutil
import atoms

def write(fileName, AL, solvant = False, optimization = False, basis='6-311G**++', method='B3LYP', mulliken=False):
    """write input jaguar file for structure in atom list AL""" 
    with open(fileName,'w') as f:
        f.write('&gen\n')
        if solvant:
            f.write('isolv=2\n')
        f.write('basis=%s\n' % basis)
        if optimization:
            f.write('igeopt=1\n')
        if mulliken:
            f.write('mulken=1\n')
        if abs(atoms.charge(AL)) > 0.001:
            f.write('molchg=%d\n' % round(atoms.charge(AL)))
        f.write('dftname=%s\n' % method)
        f.write('iacc=1\n') # 1 ultrafine, 2 fine, 3 coarse
        f.write('&\n')
        f.write('entry_name: Scratch\n')
        f.write('&zmat\n')
        for a in AL:
            if a.fixed:
                format = '%-8s %20.15f# %20.15f# %20.15f#\n'
            else:
                format = '%-8s %20.15f %20.15f %20.15f\n'
            if a.aName.strip()[1:2].islower():
                element = a.aName.strip()[0:2]
            else:
                element = a.aName.strip()[0:1] 
            f.write(format % (element,a.xyz[0],a.xyz[1],a.xyz[2]))        
        f.write('&\n')

def query(jobID):
    output = subprocess.Popen(["jaguar", "jobs", jobID], stdout=subprocess.PIPE).communicate()[0]
    lastline = output.split('\n')[-2]
    status = lastline.split()[2]
    return status

def run(fileName):
    """returns the job ID"""
    output = subprocess.Popen(["jaguar", "run", fileName], stdout=subprocess.PIPE).communicate()[0]
    return output.split()[-1]

def readMulliken(fileName, AL):
    '''bgf stores charge to 5 digit precision 
    if the sum of the charges differs from integer by at most 1e-5 per atoms, then 1e-5 is added or 
    subtracted to first few atoms to make the total charge integer'''

    lines = open(fileName).readlines()
    lineIndex = 0
    for i in range(len(lines)):
        if lines[i].startswith('  Atomic charges from Mulliken'):
            lineIndex = i
    charges = []
    lineIndex = lineIndex + 3
    while not lines[lineIndex - 1].startswith('  sum'):
        s = lines[lineIndex].split()
        charges.extend(s[1:])
        lineIndex = lineIndex + 3
    for i in range(len(AL)):
        AL[i].charge = round(1.0e5*float(charges[i]))/1.0e5
    # making the total charge integer
    charge = atoms.charge(AL)
    dcharge = charge - round(charge)
    if 0 < abs(dcharge) <= 1.0e-5 * len(AL):
        if dcharge > 0:
            sign = 1.0
        else:
            sign = -1.0
        for i in range(int(round(abs(dcharge) / 1.0e-5))):
            AL[i].charge -= sign*1.0e-5
    charge = atoms.charge(AL)
    dcharge = charge - round(charge)
    if abs(dcharge) > 1.0e-10:
        print 'Warning: total charge after reading Mullikan charges from %s is %20.15f (not interger?)' % (fileName,charge)

def readEnergy(fileName):
    """reads jaguar out file: last line starting with ' SCFE' contains the energy
    returns 0 if no such string was found"""
    energy = 0.0
    for line in open(fileName).readlines():
        if line.startswith(" SCFE"):
            s = line.split()
            token = s[-2]
            if token == 'iterations:':
                token = s[-4]
            energy = float(token)
        # some output files contain SCFE but then there is some fatal error
        if line.startswith("  ERROR: fatal error"):
            return 0.0
    return energy * 627.509469  # convert from hartrees to kcal/mol

def readCoordinates(fileName):
    """reads jaguar out file: reads the coordinates after the last '  final geometry:' occurence 
    crashes if no such string is found
    """
    lines = open(fileName).readlines()
    lineIndex = 0
    for i in range(len(lines)):
        if lines[i].startswith('  final geometry:'):
            lineIndex = i
    if not lines[lineIndex].startswith('  final geometry:'):
        print 'Error: the file %s does not contain final geometry' % fileName
        sys.exit(1)
    xyz = []
    lineIndex = lineIndex + 3
    s = lines[lineIndex].split()
    while len(s) > 0:
        x=float(s[1])
        y=float(s[2])
        z=float(s[3])
        xyz.append([x,y,z])
        lineIndex += 1
        s = lines[lineIndex].split()
    return numpy.array(xyz)

def readInputCoordinates(fileName):
    '''reads input coordinates from the out file
    I needed it once for some job checking'''
    lines = open(fileName).readlines()
    lineIndex = 0
    for i in range(len(lines)):
        if lines[i].startswith(' Input geometry:'):
            lineIndex = i
    xyz = []
    lineIndex = lineIndex + 3
    s = lines[lineIndex].split()
    while len(s) > 0:
        x=float(s[1])
        y=float(s[2])
        z=float(s[3])
        xyz.append([x,y,z])
        lineIndex += 1
        s = lines[lineIndex].split()
    return numpy.array(xyz)

def readOutCoordinates(fileName):
    '''reads the coordinates sets from the out file, as the structure is minimized, returns xyzlist and al (only with atom names)'''
    lines = open(fileName).readlines()
    lineIndex = 0
    for i in range(len(lines)):
        if lines[i].startswith(' Input geometry:'):
            lineIndex = i
    xyz = [[]]
    lineIndex = lineIndex + 3
    s = lines[lineIndex].split()
    al = []
    while len(s) > 0:
        al.append(atoms.Atom())
        al[-1].aName = s[0]
        x=float(s[1])
        y=float(s[2])
        z=float(s[3])
        xyz[0].append([x,y,z])
        lineIndex += 1
        s = lines[lineIndex].split()

    while lineIndex < len(lines):        
        if lines[lineIndex].startswith('  new geometry:') or lines[lineIndex].startswith('  final geometry:'):
            xyzitem = []
            lineIndex += 3
            for i in range(lineIndex, lineIndex+len(al)):
                s = lines[i].split()
                x=float(s[1])
                y=float(s[2])
                z=float(s[3])
                xyzitem.append([x,y,z])
            xyz.append(xyzitem)
        lineIndex += 1

    return numpy.array(xyz), al


def qsub(fileName,submit = True):
    """filename has to end with .in"""
    jobname = fileName.split('/')[-1]
    namelist = fileName.split('.')
    namelist[-1] = 'pbs'
    pbsname = '.'.join(namelist)
    f = open(pbsname,'w')
    f.write('#PBS -l nodes=1,walltime=100:00:00 \n')
    f.write('#PBS -N %s \n' % jobname)
    f.write('#PBS -S /bin/tcsh \n')
    f.write('setenv LM_LICENSE_FILE @10.254.1.1 \n')
    f.write('cd $PBS_O_WORKDIR \n')
    f.write('jaguar run -WAIT %s \n' % fileName)
    f.close()
    command = 'qsub %s' % pbsname    
    if submit:
        os.system(command)
    return command

def getMulliken(al):
    """THIS SHOULD NOT RUN ON THE FRONT NODE
    writes jaguar B3LYP / 6-31G** job to compute Mullikan charges
    and reads the charges to the atom al"""
    tmpdir = tempfile.mkdtemp(prefix='tmpJaguar',dir='.')
    os.chdir(tmpdir)
    # 
    write('tempMullikan.in', al, solvant = False, optimization = False, basis='6-31G**', method='B3LYP', mulliken=True)
    command = 'jaguar run -WAIT tempMullikan.in'
    subprocess.call(command.split())
    readMulliken('tempMullikan.out', al)
    #
    os.chdir('..')
    shutil.rmtree(tmpdir)


