""" writing QCHEM input .in and reading details from the output .out """

import subprocess, numpy, os, sys
import atoms, units

def write(fileName, AL, optimization = False, method='B3LYP', basis='6-311++G**', basis2 = False, auxbasis = False, ghostStart = -1, ghostEnd = -1, symOn = False):
    """write input qchem file for structure in atom list AL"""

    # get read filename first 
    s =  fileName.split('.')   
    if s[-1] == 'inp':
        fileNameIn = fileName
    else:
        s.append('inp')
        fileNameIn = '.'.join(s)        

    # qchem input file
    with open(fileNameIn,'w') as f:
        f.write('$molecule \n')
        if abs(atoms.charge(AL)) > 0.001:
            charge = round(atoms.charge(AL))
        else:
            charge = 0
        f.write('%d  1 \n' % charge) # second number is spin multiplicity
        for i in range(len(AL)):
            a = AL[i]
            if a.aName.strip()[1:2].islower():
                element = a.aName.strip()[0:2]
            else:
                element = a.aName.strip()[0:1] 
            if ghostStart <= i < ghostEnd:
                format = ' @%s %20.15f %20.15f %20.15f \n'
            else:
                format = ' %s %20.15f %20.15f %20.15f \n'
            f.write(format % (element,a.xyz[0],a.xyz[1],a.xyz[2]))        
        f.write('$end \n')

        f.write('$rem \n')
        if optimization:
            f.write(' jobtype opt \n')
        else:
            f.write(' jobtype sp \n')

        if method == 'B3LYP':
            f.write(' exchange B3LYP \n')
        elif method == 'XYGJOS':
            f.write(' exchange XYGJOS \n')
        elif method == 'LXYGJOS':
            f.write(' exchange LXYGJOS \n')
            f.write(' omega = 600 \n')
        elif method == 'XYG3RI':
            f.write(' exchange XYG3RI \n')
        elif method == 'MP2':
            f.write(' correlation MP2\n')
        elif method == 'LOCAL_MP2':
            f.write(' correlation LOCAL_MP2\n')
        elif method == 'RIMP2':
            f.write(' correlation RIMP2\n')
        elif method == 'RILMP2':
            f.write(' correlation RILMP2\n')
        else:
            sys.exit('Error: unknown method in qchem.write: %s ' % method)

        f.write(' basis %s \n' % basis)

        if basis2:
            f.write(' dual_basis_energy true\n')
            basis2map = {'CC-PVTZ':'RCC-PVTZ',
                         'CC-PVQZ':'RCC-PVQZ',
                         'AUG-CC-PVDZ':'RACC-PVDZ',
                         'AUG-CC-PVTZ':'RACC-PVTZ',
                         'AUG-CC-PVQZ':'RACC-PVQZ',
                         '6-31G*':'6-31G',
                         '6-311++G**':'6-311G*'}
            f.write(' basis2 %s \n' % basis2map[basis] )

        if auxbasis:
            auxbasismap =  {'CC-PVDZ':'RIMP2-CC-PVDZ',
                            'CC-PVTZ':'RIMP2-CC-PVTZ',
                            'CC-PVQZ':'RIMP2-CC-PVQZ',
                            'AUG-CC-PVDZ':'RIMP2-AUG-CC-PVDZ',
                            'AUG-CC-PVTZ':'RIMP2-AUG-CC-PVTZ',
                            'AUG-CC-PVQZ':'RIMP2-AUG-CC-PVQZ',
                            'G3LARGE':'RIMP2-AUG-CC-PVTZ'}
            f.write(' auxiliary_basis %s \n' % auxbasismap[basis])
        if not symOn:
            f.write(' sym_ignore TRUE \n')
            f.write(' symmetry FALSE \n')
        f.write('$end \n')

def writePbs(fileName, jobList = [], walltime = 100):
    """ writes a suitable pbs file for submitting the job in the queue 
    by default for one job named (parameter fileName), optionally for more jobs in a list jobList""" 
    # get read filenames first 
    s =  fileName.split('.')   
    if s[-1] == 'inp':
        fileNameIn = fileName
    elif s[-1] == 'pbs':
        s[-1] = 'inp'
        fileNameIn = '.'.join(s)
    else:
        s.append('inp')
        fileNameIn = '.'.join(s)        
    s[-1] = 'out'
    fileNameOut = '.'.join(s)
    s[-1] = 'pbs'
    fileNamePbs = '.'.join(s)
    fileNameJob = fileName.split('/')[-1]

    # pbs file
    with open(fileNamePbs,'w') as f:
        f.write('#PBS -l nodes=1:ppn=1,walltime=%d:00:00 \n' % (walltime))
        ijobname = ''.join(fileNameJob.split('.'))
        if ijobname[0].isdigit():
            ijobname = 'x' + ijobname
        f.write('#PBS -N %s \n' % ( ijobname ) )
        f.write('#PBS -S /bin/csh \n')
        f.write('cd $PBS_O_WORKDIR \n')
        # f.write('source /exec/q-chem/qcSer18Aug10/qchem.tcshrc \n')
        if len(jobList) == 0:
            f.write('qchem %s %s \n' % (fileNameIn, fileNameOut) )
        else:
            for iFile in jobList:  
                s =  iFile.split('.')   
                if s[-1] == 'inp':
                    fileNameIn = fileName
                else:
                    s.append('inp')
                    fileNameIn = '.'.join(s)        
                s[-1] = 'out'
                fileNameOut = '.'.join(s)
                f.write('qchem %s %s \n' % (fileNameIn, fileNameOut) )



'''
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
'''

def divideAtomsIntoTwoMolecules(AL):
    """ asumes the atoms are in two molecules: atoms from one molecule are grouped first and 
    the atoms from the second molecule follow
    returns the number of atoms in the first molecule"""
    group = [0 for i in range(len(AL))]
    group[0] = 1
    group[-1] = 2
    unassigned = len(AL) - 2
    maxgroup = 2
    groups = 2
    distances = []
    for i in range(len(AL)-1):
        for j in range(i+1,len(AL)):
            d = AL[i].distance(AL[j])
            distances.append([d,i,j])
    distances.sort()
    for d in distances:
        i = d[1]
        j = d[2]
        if unassigned <= 0 and groups <= 2:
            break
        if group[i] == 0 and group[j] > 0:
            group[i] = group[j]
            unassigned = unassigned - 1
        elif group[i] > 0 and group[j] == 0:
            group[j] = group[i]
            unassigned = unassigned - 1
        elif group[i] == 0 and group[j] == 0:
            maxgroup += 1
            group[j] = maxgroup
            group[i] = maxgroup
            unassigned = unassigned - 2
            groups += 1
        elif group[i] != group[j]: # merging two groups
            groups -= 1
            tempmingroup = min(group[i],group[j])
            tempmaxgroup = max(group[i],group[j])
            for k in range(len(group)):
                if group[k] == tempmaxgroup:
                    group[k] = tempmingroup

    # test if there are only two molecules
    result = 0
    for i in range(len(group)):
        if group[i] != 1:
            result = i
            break
    for i in range(result,len(group)):
        if group[i] != 2:
            print 'Error: could not find just two molecules '
            print group
            print unassigned
            print groups
            sys.exit(1)
    return result

def writeHBond(fileName, AL, method='B3LYP', basis='6-311++G**', basis2 = False, auxbasis = False, pbsFile = True, doCouterPoise = True):
    """fileName does not contain extension"""
    # find two molecules
    n = divideAtomsIntoTwoMolecules(AL)
    # write the files
    #both
    write(fileName+'.both.inp', AL, method=method, basis=basis, basis2=basis2, auxbasis=auxbasis) 
    #part1
    write(fileName+'.part1.inp', AL[0:n], method=method, basis=basis, basis2=basis2, auxbasis=auxbasis) 
    #part2
    write(fileName+'.part2.inp', AL[n:], method=method, basis=basis, basis2=basis2, auxbasis=auxbasis)
    if doCounterPoise: 
        #part1ghost 
        write(fileName+'.part1ghost.inp', AL, method=method, basis=basis, basis2=basis2, auxbasis=auxbasis, ghostStart = n, ghostEnd = len(AL)) 
        #part2ghost
        write(fileName+'.part2ghost.inp', AL, method=method, basis=basis, basis2=basis2, auxbasis=auxbasis, ghostStart = 0, ghostEnd = n) 
        # pbs 
        joblist = [fileName+'.both',fileName+'.part1',fileName+'.part2',fileName+'.part1ghost',fileName+'.part2ghost']
    else:
        # pbs 
        joblist = [fileName+'.both',fileName+'.part1',fileName+'.part2']

    if pbsFile == True:    
        writePbs(fileName, joblist=joblist)
        print 'qsub %s.pbs' % fileName
    return joblist

def readHBondEnergy(fileName):
    """ returns binding energy no BSSE, binding energy with BSSE"""
    both = readEnergy(fileName+'.both.out')
    part1 = readEnergy(fileName+'.part1.out')
    part2 = readEnergy(fileName+'.part2.out')
    part1ghost = readEnergy(fileName+'.part1ghost.out')
    part2ghost = readEnergy(fileName+'.part2ghost.out')
    if (both==0.0) or (part1==0.0) or (part2==0.0) or (part1ghost==0.0) or (part2ghost==0.0):
        print 'Warning: at least one run in hbond computation of %s did not finish' % fileName
        return [0.0,0.0]
    binding = both - part1 - part2
    bindingCP = both - part1ghost - part2ghost
    return [binding,bindingCP]

def readEnergy(fileName):
    """qchem has stupidly formatted total energy differently for every different method"""
    try:
        lines = open(fileName).readlines()
    except IOError:
        print 'Warning: the file %s does not exist ' % fileName
        return 0.0
        
    # test if the run finished ok
    ok = False
    if len(lines) < 11:
        ok = False
    else:
        if lines[-5]=='        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
            ok = True
        elif lines[-4]==' Q-Chem fatal error occurred in module saveScratchFiles.C, line 50:\n' and lines[-10]=='        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
            ok = True
    if not ok:
        print 'Warning: run %s did not finish' % fileName
        return 0.0

    # find out method first
    method = 'undefined'
    for line in lines:
        if line.startswith(' exchange B3LYP'):
            method = 'B3LYP'
        elif line.startswith(' exchange XYGJOS'):
            method = 'XYGJOS'
        elif line.startswith(' exchange LXYGJOS'):
            method = 'LXYGJOS'
        elif line.startswith(' exchange XYG3RI'):
            method = 'XYG3RI'
        elif line.startswith(' correlation RIMP2'):
            method = 'RIMP2'
        elif line.startswith(' correlation RILMP2'):
            method = 'RILMP2'

    # get the energy
    result = 0.0    
    if method == 'B3LYP' :
        for line in lines:
            if line.startswith(' Total energy in the final basis set ='):
                result = float(line.split()[-1])
    elif method == 'XYGJOS' :
        for line in lines:
            if line.startswith(' Total XYGJ-OS energy                    ='):
                result = float(line.split()[-1])
    elif method == 'LXYGJOS' :
        for line in lines:
            if line.startswith(' Total Local XYGJ-OS energy               ='):
                result = float(line.split()[-1])
    elif method == 'XYG3RI' :
        for line in lines:
            if line.startswith(' Total  XYG3-RI     total energy  ='):
                result = float(line.split()[-2])
    elif method == 'RIMP2' :
        for line in lines:
            if line.startswith('        RIMP2         total energy ='):
                result = float(line.split()[-2])
    elif method == 'RILMP2' :
        for line in lines:
            if line.startswith('TOTAL TRIM ENERGY                       ='):
                result = float(line.split()[-2])
    else:
        print 'Error: unknown QM method in qchem output file: %s ' % method
        sys.exit(1)

    return units.hartree2kcalmol( result )



