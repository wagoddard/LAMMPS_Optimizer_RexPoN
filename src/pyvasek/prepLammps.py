"""
Considered pair styles:

lj/charmm/coul/charmm
lj/charmm/coul/long/opt
-- only these two pair style mix

buck/coul/cut
buck/coul/long
-- only the above 4 styles are non hybrid

morse/opt
coul/cut
coul/long 

hbond/dreiding/lj
hbond/dreiding/morse

hybrid/overlay is used to combine these interactions
"""

import sys, numpy, os, time, math
import structure, atoms, periodic

def generateLammpsData(al, ff, fileBase, cell=None,debug=False,exceptResidues=[]):
    """denerates the data file with all the information for Lammps
    only rectangular cells are supported
    so  cell is a list of 
    tree values [xrange, yrange,zrange] with center 0,0,0
    or 
    4x3 numpy array with the cell size in the NAMD notation
    """
    filescript = fileBase+'.in'
    filedata = fileBase+'.data'
    if debug: print 'getting ready to write lammps input and datafiles: %s, %s' %(filescript, filedata)
    
    st = structure.Structure(al,ff,debug=debug,exceptResidues=exceptResidues)
    # find out which angles and types are actually used

    # if the molecule is very flat bin style for neighbor crashes, so use nsq
    useNsqNeighbor = False
    xyz = atoms.getXyzFromList(al)
    xyzmax = numpy.max(xyz,axis=0)
    xyzmin = numpy.min(xyz,axis=0)
    xyzRange = xyzmax-xyzmin
    if min(xyzRange) < 2.5: useNsqNeighbor = True

    # periodicity considerations
    if cell == None:
        isPeriodic = False
        cellBoundary = 's s s'
        xyzmin = xyzmin - ff.cutoffNonBond
        xyzmax = xyzmax + ff.cutoffNonBond
        cell = [xyzmin[0],xyzmax[0],xyzmin[1],xyzmax[1],xyzmin[2],xyzmax[2]]
    else:
        isPeriodic = True
        cellBoundary = 'p p p'
        if len(cell) == 3:
            boxSize = numpy.zeros((4,3))
            boxSize[0,0] = cell[0]
            boxSize[1,1] = cell[1]
            boxSize[2,2] = cell[2]
        elif len(cell) == 4:
            boxSize = cell
        else: 
            sys.exit('Error: prepLammps in generateLammpsData: wrong input for cell dimensions (has to have length 3 or 4, but has %d' % len(cell))
        # need to make sure all the atoms are in the main cell
        # make array with the box size in the NAMD notation
        periodic.putCoordinatesIntoBox(al,boxSize)
        # create list of 6 or 9 numbers xmin,xmax,xy... to specify lammps cell size
        # lammps: "origin" at (xlo,ylo,zlo) 
        # 3 edge vectors starting from the origin given by 
        # A = (xhi-xlo,0,0); 
        # B = (xy,yhi-ylo,0); 
        # C = (xz,yz,zhi-zlo).
        # check that the box has 0's in the right places:
        if boxSize[0,1]!=0 or boxSize[0,2]!=0 or boxSize[1,2]!=0 or boxSize[0,0] <= 0. or boxSize[1,1] <= 0.0 or boxSize[2,2] <= 0.:
            sys.exit('Error: prepLammps: cell vector A is not aligned with +x or B is not in xy plane: %s , but these are required by lammps' % str(boxSize))
        # lammps has pretty weird definition of the periodic cell, just look at hte lammps online page to understand this...
        origin = boxSize[3,:] - (boxSize[0,:]+boxSize[1,:]+boxSize[2,:])/2
        xlo = origin[0]
        ylo = origin[1]
        zlo = origin[2]
        lx = boxSize[0,0]
        ly = boxSize[1,1]
        lz = boxSize[2,2]
        xy = boxSize[1,0]
        xz = boxSize[2,0]
        yz = boxSize[2,1]
        xhi = xlo+lx
        yhi = ylo+ly
        zhi = zlo+lz
        cell = [xlo ,xhi ,ylo ,yhi ,zlo ,zhi ]
        if xy!=0 or xz!=0 or yz!=0:
            cell.append(xy )
            cell.append(xz )
            cell.append(yz )

    # variables describing pair styles:
    #
    # usedPairStyle -- vdw requested in the par file, morse or lj12-6
    # isPeriodic
    # doHydrogenBonds
    #   myHBondType - only defined if doHydrogenBonds is True
    # pairStyleHybrid
    # pairStyleVDW
    #   pairStyleCoul - only defined if usedPairStyle == 'morse'

    # vdwtype
    usedPairStyle = ff.atomTypes[st.sortedAtomLabels[0]].vdwType
    # hbond type    
    if len(st.sortedHydrogenBondTypes) > 0 and ff.doHbonds:
        doHydrogenBonds = True
        if st.sortedHydrogenBondTypes[0].type == 'morse':
            myHBondType = 'hbond/dreiding/morse'
        elif st.sortedHydrogenBondTypes[0].type == 'lj12-10':
            myHBondType = 'hbond/dreiding/lj'
        else:
            sys.exit('Erorr: unknown hydrogen bond type: %s' % (st.sortedHydrogenBondTypes[0].type))
    else:
        doHydrogenBonds = False
    # is this hybrid pair style?
    if (not doHydrogenBonds) and usedPairStyle in ['lj12-6', 'exp6']:
        pairStyleHybrid = False
    else:
        pairStyleHybrid = True

    # prepare names
    if usedPairStyle == 'lj12-6':
        if isPeriodic:
            pairStyleVDW = 'lj/charmm/coul/long/opt'
        else:
            pairStyleVDW = 'lj/charmm/coul/charmm'
    elif usedPairStyle == 'exp6':
        if isPeriodic:
            pairStyleVDW = 'buck/coul/long'
        else:
            pairStyleVDW = 'buck/coul/cut'
    elif usedPairStyle == 'morse':
        pairStyleVDW = 'morse/opt'
        if isPeriodic:
            pairStyleCoul = 'coul/long'
        else:
            pairStyleCoul = 'coul/cut'
    else:
        sys.exit('Erorr: unknown VDW type in prepLammps.py (only lj12-6 and morse available now): %s' % usedPairStyle)

    # PREPROCESS EXPLICIT MIXING
    # hybrid pair does not mix so we need to check if we are doing hydrogen bonds
    
    pair_coeff_lines = []
    if (not pairStyleHybrid) and usedPairStyle == 'lj12-6': # we only have to list offdiagonal VDW when it is requeted in the parameter file
        for i,index in enumerate(st.sortedAtomLabels):
            for j in range(i,len(st.sortedAtomLabels)):
                jindex =  st.sortedAtomLabels[j]
                # first see, if we have to use off-diagonal rule, if not just compute it
                if index+jindex in ff.offDiagVDWTypes:
                    params = ff.offDiagVDWTypes[index+jindex].lammps()
                    pair_coeff_lines.append('pair_coeff %d %d %s # %s-%s\n' % (i+1,j+1,params,index,jindex) )
                elif jindex+index in ff.offDiagVDWTypes:
                    params = ff.offDiagVDWTypes[jindex+index].lammps()
                    pair_coeff_lines.append('pair_coeff %d %d %s # %s-%s\n' % (i+1,j+1,params,index,jindex) ) 
    else: # for hybrid pair style we have to list the offdiagonal VDW always
        if usedPairStyle == 'morse':  # for morse coulomb is not included in the pair style but via hybrid
            pair_coeff_lines.append('pair_coeff * * %s\n' % pairStyleCoul)
        for i,index in enumerate(st.sortedAtomLabels):
            for j in range(i,len(st.sortedAtomLabels)):
                jindex =  st.sortedAtomLabels[j]
                # first see, if we have to use off-diagonal rule, if not just compute it
                if index+jindex in ff.offDiagVDWTypes:
                    params = ff.offDiagVDWTypes[index+jindex].lammps()
                elif jindex+index in ff.offDiagVDWTypes:
                    params = ff.offDiagVDWTypes[jindex+index].lammps()
                else: # computing the diagonal rule here                
                    bi = ff.atomTypes[index]  # these are atom
                    bj = ff.atomTypes[jindex] 
                    params = bi.mixWith(bj,ff.mean).lammps()
                if pairStyleHybrid:
                    pair_coeff_lines.append('pair_coeff %d %d %s %s # %s-%s\n' % (i+1,j+1,pairStyleVDW,params,index,jindex) )    
                else:
                    pair_coeff_lines.append('pair_coeff %d %d %s # %s-%s\n' % (i+1,j+1,params,index,jindex) )    
    if doHydrogenBonds:  
        myHBondPower = st.sortedHydrogenBondTypes[0].power
        for b in st.sortedHydrogenBondTypes:
            bd = st.usedAtomLabels[b.donor]+1
            bh = st.usedAtomLabels[b.hydrogen]+1
            ba = st.usedAtomLabels[b.acceptor]+1
            if bd <= ba:
                flag = 'i'
                bi = bd
                bj = ba
            else:
                flag = 'j'
                bi = ba
                bj = bd
            pair_coeff_lines.append('pair_coeff %d %d %s %d %s %s # donor:%s hydrogen:%s acceptor:%s \n'% (bi, bj, myHBondType, bh, flag, b.lammps(),b.donor,b.hydrogen,b.acceptor))

    # OUTPUTTING
    # output the stuff
    if debug: print 'writing lammps input and datafile'
    # the script file
    with open(filescript,'w') as script:
        script.write('boundary        %s \n' % cellBoundary)
        script.write('units           real \n')
        if useNsqNeighbor: script.write('neighbor        2.0 nsq\n')
        if ff.cutoffNonBond > 15.0:
            # the default number of neighbors in lammps is 2000
            # the number of atoms in radius r is ~0.4*r^3, let's use 50% more
            estimate = int( 0.6*math.pow(ff.cutoffNonBond,3.0) )
            bound = len(al)+10
            neighborone = max(2000, min(estimate, bound) )
            neighborpage = 20 * neighborone 
            script.write('neigh_modify    one %d page %d\n' % (neighborone, neighborpage))
            pass
        #script.write('neigh_modify    delay 10 every 1 check yes one 4000 page 40000\n')
        script.write(' \n')
        script.write('atom_style      full \n')
        
        # BOND STYLE
        if len(st.usedBondTypes) > 0:        script.write('bond_style      harmonic\n')
        else:                             script.write('bond_style      none \n')
        
        # ANGLE STYLE
        if len(st.usedAngleTypes) > 0:
            if st.sortedAngleTypes[0].type == 'harmonic':       
                                          script.write('angle_style     harmonic \n')
            else:                         # cosine angles as in dreiding 
                                          script.write('angle_style     hybrid cosine/periodic cosine/squared \n')
        else:                             script.write('angle_style     none \n')
        
        # TORSION STYLE
        if len(st.usedTorsionTypes) > 0:  script.write('dihedral_style  harmonic \n')
        else:                             script.write('dihedral_style  none \n')
        
        # INVERSION STYLE
        if len(st.usedInversionTypes) > 0:   
            if st.sortedInversionTypes[0].type == 'amber':
                                          script.write('improper_style  cvff \n')
            else:                         # dreiding inversion
                                          script.write('improper_style  umbrella \n')                
        else:                             script.write('improper_style  none \n')
        script.write(' \n')
        
        # PAIR STYLE
        if not pairStyleHybrid:
            if usedPairStyle == 'lj12-6':
                script.write('pair_style      %s %f %f \n' % (pairStyleVDW, ff.splineNonBond, ff.cutoffNonBond))                
                script.write('pair_modify     mix %s \n' % ff.mean)  # note that hte hybrid style does not do mixing
            elif usedPairStyle == 'exp6':
                script.write('pair_style      %s %f \n' % (pairStyleVDW, ff.cutoffNonBond))                
            else:
                sys.exit('Error in prepLammps. One should never get to this part of code, but usedPairStyle is :%s' % usedPairStyle)
        else:
            pairStyleLine = 'pair_style      hybrid/overlay'
            if usedPairStyle == 'lj12-6':
                pairStyleLine += ' %s %f %f' % (pairStyleVDW, ff.splineNonBond, ff.cutoffNonBond)
            elif usedPairStyle == 'morse':
                pairStyleLine += ' %s %f' % (pairStyleVDW, ff.cutoffNonBond)
                pairStyleLine += ' %s %f' % (pairStyleCoul, ff.cutoffNonBond)
            else:
                sys.exit('Erorr: unknown VDW type in prepLammps.py (only lj12-6 and morse available now): %s' % usedPairStyle)
            if doHydrogenBonds:
                pairStyleLine += ' %s %d %f %f %f \n' % (myHBondType, myHBondPower,ff.splineHBond,ff.cutoffHBond,ff.angleHBond)
            pairStyleLine += ' \n'
            script.write(pairStyleLine)
            
        if isPeriodic:
            script.write('kspace_style    pppm 1e-4 \n')
        script.write('dielectric      %f \n' % ff.dielectric)
    
        if ff.special14lj == 1.0 and ff.special14coul == 1.0:
            script.write('special_bonds   dreiding \n')
        elif ff.special14lj == 0.5 and ff.special14coul == 1.0/1.2:
            script.write('special_bonds   amber \n')
        else:             
            print >> sys.stderr, 'special14lj %f and special14coul %f is neither dreiding nor amber ' % (ff.special14lj,ff.special14coul)
            sys.exit(1)
    
        script.write(' \n')
        script.write('read_data       %s \n' % filedata)
        script.write(' \n')
        script.write('#### alternatively read coordinates from the latest (*) restart file \n')
        script.write('# read_restart    restart.*\n')
        script.write(' \n')
        for line in pair_coeff_lines:
            script.write(line) 
    
        # write the run information into the script file
        # define computes
        
        script.write(' \n')
        script.write('variable        step equal step \n')
        script.write('variable        ebond equal ebond \n')
        script.write('variable        eangle equal eangle \n')
        script.write('variable        edihed equal edihed \n')
        script.write('variable        eimp equal eimp \n')
        script.write('variable        emol equal emol \n')
        script.write('variable        ecoul equal ecoul \n')
        
        # the values reported by hbond computes are not correct, so have to use some global variables
        # try 1 after update
        #if doHydrogenBonds:
        #    script.write('compute hbond all pair hbond/dreiding/%s \n' % myHBondType)
        #    script.write('variable hbond equal c_hbond[2] \n')
        #    script.write('variable counthbond equal c_hbond[1] \n')
        #else:
        #    script.write('variable hbond equal 0.0 \n')
        #    script.write('variable counthbond equal 0.0 \n')
        #script.write('variable evdwl equal evdwl-v_hbond \n') # NOTE: thermo evdwl already contains hbond
        # try 2       
        script.write('compute         evdwl all pair %s evdwl \n' % pairStyleVDW)    
        script.write('variable        evdwl equal c_evdwl \n') # NOTE: thermo evdwl already contains hbond
        if doHydrogenBonds:
            script.write('compute         hbondevdwl all pair %s evdwl \n' % myHBondType)
            script.write('variable        hbond equal c_hbondevdwl \n')
            script.write('compute         hbond all pair %s \n' % myHBondType)
            script.write('variable        counthbond equal c_hbond[1] \n')
        else:
            script.write('variable        hbond equal 0.0 \n')
            script.write('variable        counthbond equal 0.0 \n')

        script.write('variable        elong  equal elong \n')
        script.write('variable        epair equal epair \n')
        script.write('variable        pe equal pe \n')
        script.write('variable        ke equal ke \n')
        script.write('variable        etotal equal etotal \n')
        script.write('variable        temp equal temp \n')
        script.write('variable        press equal press \n')
        script.write('variable        fmax equal fmax \n')
        
        script.write(' \n')
        script.write('thermo          100 \n')
        script.write('thermo_style    custom step ebond eangle edihed eimp emol ecoul v_evdwl v_hbond v_counthbond elong epair pe ke etotal temp press fmax\n')
        script.write(' \n')
        # fixed atoms
        # there is 2048 characters per line limit in lammps
        fixed = [a for a in al if a.fixed]
        if len(fixed) > 0:
            groupLines = makeGroupLines(fixed,'freezeatoms')
            script.write(groupLines)
            script.write('fix             freeze freezeatoms setforce 0.0 0.0 0.0')
        script.write('\n')
        # sample commands
        script.write('\n')        
        script.write('#### output the (wrapped) coordinates \n')
        script.write('# dump            1 all custom 100 dump.lammpstrj id type x y z vx vy vz \n')
        script.write('# dump_modify     1 sort id \n')
        script.write('\n')
        script.write('#### minimization \n')
        script.write('# minimize        0.0 0.0 100 10000 \n')
        script.write('\n')
        script.write('#### initialization of velocities \n')
        script.write('# velocity        all create 300.0 4928459 mom yes rot yes dist gaussian\n')
        script.write('\n')
        script.write('#### nve run \n')
        script.write('# timestep        1 \n')
        script.write('# fix             1 all nve \n')
        script.write('# run             1000 \n')
        script.write(' \n')  
        script.write('#### write restart file (* gets replaced by timestep)\n')
        script.write('# write_restart   restart.*\n')
        script.write(' \n')  


    # writing the data file
    with open(filedata,'w') as out:
        # header of the data file
        out.write('Generated by prepLammps.py, %s@%s on %s\n\n' %(os.getenv('LOGNAME'),os.getenv('HOSTNAME'),time.strftime('%X %x %Z')))
        out.write('%d atoms\n' %len(al))
        out.write('%d bonds\n' %len(st.bonds))
        out.write('%d angles\n' %len(st.angles))
        out.write('%d dihedrals\n' %len(st.torsions))
        out.write('%d impropers\n\n' %len(st.inversions))
    
        out.write('%d atom types\n' %len(st.sortedAtomLabels))
        out.write('%d bond types\n' %len(st.sortedBondTypes))
        out.write('%d angle types\n' %len(st.sortedAngleTypes))
        out.write('%d dihedral types\n' %len(st.sortedTorsionTypes))
        out.write('%d improper types\n\n' %len(st.sortedInversionTypes))
    
        out.write('%f %f xlo xhi\n' % (cell[0],cell[1]))
        out.write('%f %f ylo yhi\n' % (cell[2],cell[3]))
        out.write('%f %f zlo zhi\n' % (cell[4],cell[5]))
        if len(cell) == 9:
            out.write('%f %f %f xy xz yz\n' % (cell[6],cell[7],cell[8]))
    
        # atom types
        out.write('\nMasses\n\n')
        for i,b in enumerate(st.sortedAtomLabels):
            bt = ff.atomTypes[b]
            out.write('%d %f # %s\n'%(i+1,bt.mass,b))            
    
        # vdw parameters and hydrogen bonds parameters
        # bond Types
        if not pairStyleHybrid:
            out.write('\nPair Coeffs\n\n')
            for i,b in enumerate(st.sortedAtomLabels):
                bt = ff.atomTypes[b]
                out.write('%d %s # %s\n'%(i+1,bt.lammps(),b))
    
        # bond Types
        if len(st.sortedBondTypes) > 0:
            out.write('\nBond Coeffs\n\n')
            for i,bi in enumerate(st.sortedBondTypes):
                out.write('%d %s\n'%(i+1, bi.lammps()))
     
        # angle Types
        if len(st.sortedAngleTypes) > 0:
            out.write('\nAngle Coeffs\n\n')
            for i,bi in enumerate(st.sortedAngleTypes):
                out.write('%d %s\n'%(i+1, bi.lammps()))
        
        # torsion Types
        if len(st.sortedTorsionTypes) > 0:
            out.write('\nDihedral Coeffs\n\n')
            for i,bi in enumerate(st.sortedTorsionTypes):
                out.write('%d %s\n'%(i+1, bi.lammps()))

        # inversion Types
        if len(st.sortedInversionTypes) > 0:
            out.write('\nImproper Coeffs\n\n')
            for i,bi in enumerate(st.sortedInversionTypes):
                out.write('%d %s\n'%(i+1, bi.lammps()))
    
        # atoms
        out.write('\nAtoms\n\n')
        for a in al:
            atomtypeindex = st.usedAtomLabels[a.ffType]+1
            out.write('%d 0 %d %.10f %.10f %.10f %.10f # %s\n' % (a.aNo,atomtypeindex,a.charge,a.xyz[0],a.xyz[1],a.xyz[2],a))
    
        # bonds
        out.write('\nBonds\n\n')
        for i,bi in enumerate(st.bonds):
            out.write('%d %d %d %d\n'%(i+1, st.usedBondTypes[bi.type]+1, bi.i.aNo, bi.j.aNo))
    
        # angles
        out.write('\nAngles\n\n')
        for i,bi in enumerate(st.angles):
            out.write('%d %d %d %d %d\n'%(i+1, st.usedAngleTypes[bi.type]+1, bi.i.aNo, bi.j.aNo, bi.k.aNo))

        # torsion
        out.write('\nDihedrals\n\n')
        for i,bi in enumerate(st.torsions):
            out.write('%d %d %d %d %d %d\n'%(i+1, st.usedTorsionTypes[bi.type]+1, bi.i.aNo, bi.j.aNo, bi.k.aNo, bi.l.aNo))
    
        # inversions
        out.write('\nImpropers\n\n')
        for i,bi in enumerate(st.inversions):
            if bi.type.type == 'amber':
                out.write('%d %d %d %d %d %d\n'%(i+1, st.usedInversionTypes[bi.type]+1, bi.j.aNo, bi.k.aNo, bi.i.aNo, bi.l.aNo))
            elif bi.type.type == 'dreiding':
                out.write('%d %d %d %d %d %d\n'%(i+1, st.usedInversionTypes[bi.type]+1, bi.i.aNo, bi.j.aNo, bi.k.aNo, bi.l.aNo))
            else:
                print >> sys.stderr, 'Unknown improper type for lammps (probably unimplemented yet) %s ' % bi.type.type
                sys.exit(1)   

    if debug: print 'done writing lammps input and datafiles: %s, %s' %(filescript, filedata)
    return st

#############
# Helping functions

def makeGroupLines(atomlist, name):
    '''defines a group from the atomlist, so that no line of the definition is longer than 2040 characters'''
    al = atomlist[:]
    groupLines = ''
    maxLength = 2030
    while len(al) > 0:
        groupLine = 'group %s id' % name
        while len(al) > 0 and len(groupLine) < maxLength:
            a = al.pop(0)
            groupLine += ' %d'%a.aNo
        groupLines += groupLine + '\n' 
    return groupLines

'''
def extractComputes(lmp):
    names = ['bond','angle','dihedral','improper','ecoul','kspace','evdwl','hbond','pe']
    for name in names:
        try:
            compute = lmp.extract_compute(name,0,0)
        except ValueError:
            compute = 0.0 # is the energy term is not defined, let's say it is 0
        print '  %s: %f' % (name,compute) 
        '''
'''

for name in glob.glob('*.bgf'):
    print name 

    lammpsin = os.tempnam(tempdir, 'lammpsin') 
    lammpsdata = os.tempnam(tempdir, 'lammpsdata')
    
    try:
        AL = vAtom.readBgfFile(name, debug = debug)
    except ValueError:
        print 'Error parsing bgf file %s' %name
        continue
    try:
        vEnergy.generateLammpsData(AL,par,lammpsin,lammpsdata, debug = debug)
    except:
        print 'Error generating structure for %s' % name
        continue
    if debug:
        lmp = lammps.lammps(['-echo','screen'])
    else:
        lmp = lammps.lammps(['-screen','none'])
    lmp.file(lammpsin) 

    vEnergy.extractComputes(lmp)

    #os.remove(lammpsin)
    #os.remove(lammpsdata)
'''        
