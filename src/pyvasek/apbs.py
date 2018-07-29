""" compute interaction energy using apbs
(note uses amber radii in pqr)

apbs input file:
=====================
read
    mol pqr complex.pqr
    mol pqr protein.pqr
    mol pqr ligand.pqr
end
elec name complex_solv
    mg-manual
    dime 129 129 161             # Specifies the number of grid points per processor for grid-based discretization, should be 65, 97, 129, and 161...
    glen 80 80 100               # Specify the mesh domain lengths for multigrid mg-manual calculations.
    gcent 0 0 0                  # Specify the center of the grid based on a molecule's center or absolute coordinates for a mg-manual multigrid calculation.
    mol 1                        # Specify the molecule. IDs are based on the order in which molecules are read by READ mol statements, starting from 1.
    npbe                         # Specifies that the nonlinear (full) Poisson-Boltzmann equation should be solved.
    bcfl sdh                     # "Single Debye-Huckel" boundary condition. 
    # ion charge 1 conc 0.225 radius 2.0   # conc = Mobile ion species concentration (floating point number in M)
    # ion charge -1 conc 0.225 radius 2.0
    pdie 4.0                     # Specify the dielectric constant of the biomolecule. This is usually a value between 2 to 20,
    sdie 80.0                    # Specify the dielectric constant of the solvent.
    srfm mol                     # Specify the model used to construct the dielectric and ion-accessibility coefficients. 
    chgm spl2                    # Specify the method by which the point charges (i.e., Dirac delta functions) by which charges are mapped to the grid
    sdens 10.00                  # Specify the number of grid points per square-angstrom to use in discontinuous surface constructions 
    srad 1.40                    # Specify the radius of the solvent molecules, usually set to 1.4 A for water.
    swin 0.30                    # Specify the size of the support (i.e., the rate of change) for spline-based surface definitions (see srfm). Usually 0.3 A.
    temp 298.15
    calcenergy total
    calcforce no
end
elec name complex_ref
    ...
    sdie 1.0
    ...
end
...
# Solvation energy change due to binding
print energy complex_solv - complex_ref - protein_solv + protein_ref - ligand_solv + ligand_ref end
# Coulomb energy of binding in vacuum
print energy complex_ref - protein_ref - ligand_ref end
# Energy change due to binding
print energy complex_solv - protein_solv - ligand_solv end
quit

"""

import numpy, os, subprocess, sys, tempfile, shutil
from pyvasek import pqr, units, atoms, bgf

stringHeader = """read
   mol pqr %s
   mol pqr %s
   mol pqr %s
end
"""
stringBody = """elec name %s_%s
   mg-manual
   dime %d %d %d
   glen %f %f %f
   gcent 0 0 0 
   mol %d 
   npbe
   bcfl sdh 
   ion charge 1 conc 0.1 radius 2.0 
   ion charge -1 conc 0.1 radius 2.0
   pdie 4.0 
   sdie %f 
   srfm mol 
   chgm spl2 
   sdens 10.00 
   srad 1.40   
   swin 0.30  
   temp 298.15
   calcenergy total
   calcforce no
end
"""
stringEndBind ="""# Energy change due to binding
print elecEnergy complex_solv - protein_solv - ligand_solv end
quit
"""
stringEndSeparate ="""# Solvation energy change due to binding
print elecEnergy complex_solv - complex_ref - protein_solv + protein_ref - ligand_solv + ligand_ref end
# Coulomb energy of binding in vacuum
print elecEnergy complex_ref - protein_ref - ligand_ref end
# Energy change due to binding
print elecEnergy complex_solv - protein_solv - ligand_solv end
# Solvation of the ligand
print elecEnergy ligand_solv - ligand_ref end
quit
"""

def apbsInteraction(al, ligandResidueList, boxSize = [80.,80.,100.], gridSize = [129, 129, 161], separate=True, debug=False, refDielectric = 1.0, amberRadii=True):
    """if separate == True this returns dictionary of three interaction energies
    between the ligandResidue (typically [999] or [444]):
    { dGsolv = charged solvation energy = complex_solv - complex_ref - protein_solv + protein_ref - ligand_solv + ligand_ref
      dGcoul = coulomb interactoin in vaccum = complex_ref - protein_ref - ligand_ref
      dGbind = binding energy = dGsolv + dGcoul = complex_solv - protein_solv - ligand_solv
      dGligsolv = ligand_solv - ligand_ref
    }
    if separate == False, then dictionary with only one energy
    {dGbind} is returned

    note: 
    - nonpolar solvation energy is neglected; apbs has some ways to estimate it too
    - the coordinates of atoms in the atom list are moved to be as far from the box boundary as possible
    - atoms have to by typed in amber force field, because the atomic radii are taken from there
    - if amberRadii == False will use Dreiding vdw radii (using Dreiding does not work yet)
    """
    nBoxSize = numpy.array(boxSize)

    xyz = atoms.getXyzFromList(al)
    xyzmax = numpy.max(xyz,0)
    xyzmin = numpy.min(xyz,0)
    if sum( (xyzmax-xyzmin) >= (nBoxSize-30.)) > 0:
        sys.exit('Error in apbsInteraction: the molecule does not have at leat 15A gap on each side: increase the boxSize (%s)' % (xyzmax-xyzmin))
    offset = (xyzmax + xyzmin)/2.
    xyz -= offset # this already updates atoms
    # atoms.putXyzIntoList(al,xyz) # not necessary because of the previous command

    ligand = [a for a in al if a.rNo in ligandResidueList]
    protein = [a for a in al if a.rNo not in ligandResidueList]
    if len(ligand)==0 or len(protein)==0:
        sys.exit('Error in apbsInteraction: invalid selection of the ligand residues: ligand has %d atoms out of %d' % (len(ligand),len(al)))

    fcom = tempfile.NamedTemporaryFile(mode='w', prefix='tmpApbsComplex', suffix='.pqr', dir='.', delete=False)
    fpro = tempfile.NamedTemporaryFile(mode='w', prefix='tmpApbsProtein', suffix='.pqr', dir='.', delete=False)
    flig = tempfile.NamedTemporaryFile(mode='w', prefix='tmpApbsLigand', suffix='.pqr', dir='.', delete=False)
    fin  = tempfile.NamedTemporaryFile(mode='w', prefix='tmpApbsCommand', suffix='.in', dir='.', delete=False)
    fcom.close()
    fpro.close()
    flig.close()
    fin.close()    
    fcom = fcom.name
    fpro = fpro.name
    flig = flig.name
    fin  = fin.name

    pqr.write(fcom,al,amberRadii)
    pqr.write(fpro,protein,amberRadii)
    pqr.write(flig,ligand,amberRadii)

    with open(fin,'w') as f:
        # header
        f.write(stringHeader % (fcom, fpro, flig)) 
        # body
        entries = ['complex','protein','ligand']
        for i in range(len(entries)):
            f.write(stringBody % (entries[i],'solv',gridSize[0],gridSize[1],gridSize[2],boxSize[0],boxSize[1],boxSize[2],i+1,80.) )
            if separate:
                f.write(stringBody % (entries[i],'ref',gridSize[0],gridSize[1],gridSize[2],boxSize[0],boxSize[1],boxSize[2],i+1,refDielectric) )
        # end
        if separate:
            f.write(stringEndSeparate)
        else:
            f.write(stringEndBind)

    # run apbs
    apbsrun = subprocess.Popen(['apbs', fin], stdout=subprocess.PIPE)
    stdoutstring = apbsrun.communicate()[0]

    dGsolv = 0.
    dGcoul = 0.
    dGbind = 0.
    dGligsolv = 0.

    # parse the output
    lines = stdoutstring.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        if line=='print energy 1 (complex_solv) - 2 (complex_ref) - 3 (protein_solv) + 4 (protein_ref) - 5 (ligand_solv) + 6 (ligand_ref) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGsolv = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 2 (complex_ref) - 4 (protein_ref) - 6 (ligand_ref) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGcoul = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 1 (complex_solv) - 3 (protein_solv) - 5 (ligand_solv) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGbind = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 1 (complex_solv) - 2 (protein_solv) - 3 (ligand_solv) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGbind = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 5 (ligand_solv) - 6 (ligand_ref) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGligsolv = units.joule2calorie( float(s[-2]) )
        i+=1 

    if not debug:
        os.remove(fcom)
        os.remove(fpro)
        os.remove(flig)
        os.remove(fin)
    if separate:
        return {'dGsolv':dGsolv, 'dGcoul':dGcoul, 'dGbind':dGbind, 'dGligsolv': dGligsolv}
    else:
        return {'dGbind':dGbind}

stringBodyMembranePrep = """elec name %s_%s
   mg-dummy
   dime %d %d %d
   glen %f %f %f
   gcent 0 0 0 
   mol %d 
   npbe
   bcfl sdh 
   ion charge 1 conc 0.1 radius 2.0 
   ion charge -1 conc 0.1 radius 2.0
   pdie 4.0 
   sdie %f 
   srfm mol 
   chgm spl2 
   sdens 10.00 
   srad 1.40   
   swin 0.30  
   temp 298.15
   calcenergy total
   calcforce no
   write dielx dx dielx_%s		   # spatially dependent diel. constant
   write diely dx diely_%s		   # out using the Conoly/Molecular surfac
   write dielz dx dielz_%s
   write kappa dx kappa_%s		   # write out the kappa map
end
"""

stringHeaderMembrane = """read
    mol pqr %s
    mol pqr %s
    mol pqr %s

    diel dx dielx_cm.dx diely_cm.dx dielz_cm.dx
    diel dx dielx_pm.dx diely_pm.dx dielz_pm.dx
    diel dx dielx_lm.dx diely_lm.dx dielz_lm.dx

    kappa dx kappa_cm.dx
    kappa dx kappa_pm.dx
    kappa dx kappa_lm.dx
end
"""

stringBodyMembrane = """elec name %s_%s
   mg-manual
   dime %d %d %d
   glen %f %f %f
   gcent 0 0 0 
   mol %d 
   npbe
   bcfl sdh 
   ion charge 1 conc 0.1 radius 2.0 
   ion charge -1 conc 0.1 radius 2.0
   pdie 4.0 
   sdie %f 
   srfm mol 
   chgm spl2 
   sdens 10.00 
   srad 1.40   
   swin 0.30  
   temp 298.15
   calcenergy total
   calcforce no
   usemap diel %d
   usemap kappa %d
end
"""

def apbsInteractionMembrane(al, ligandResidueList, boxSize = [80.,80.,100.], gridSize = [129, 129, 161], debug=False, refDielectric = 1.0, amberRadii=True):
    """returns dictionary of three interaction energies
    between the ligandResidue (typically [999] or [444]):
    { dGsolv = charged solvation energy = complex_solv - complex_ref - protein_solv + protein_ref - ligand_solv + ligand_ref
      dGcoul = coulomb interactoin in vaccum = complex_ref - protein_ref - ligand_ref
      dGbind = binding energy = dGsolv + dGcoul = complex_solv - protein_solv - ligand_solv
      dGligsolv = ligand_solv - ligand_ref
    }

    note: 
    - nonpolar solvation energy is neglected; apbs has some ways to estimate it too
    - the coordinates of atoms in the atom list are moved so that the center of hte protein is at x=0, y=0, the protein has to be aligned to membrane
    - atoms have to by typed in amber force field, because the atomic radii are taken from there
    - if amberRadii == False will use Dreiding vdw radii (using Dreiding does not work yet)
    """
    nBoxSize = numpy.array(boxSize)

    xyz = atoms.getXyzFromList(al)
    xyzmean = numpy.mean(xyz,0)
    xyzmean[2] = 0.0 # do not change the z position
    xyz -= xyzmean # this already updates atoms
    # atoms.putXyzIntoList(al,xyz) # not necessary because of the previous command

    xyzmax = numpy.max(xyz,0)
    xyzmin = numpy.min(xyz,0)
    if sum( xyzmax >= (nBoxSize/2.-15.)) > 0 or sum( xyzmin <= (-nBoxSize/2.+15.)) > 0:
        sys.exit('Error in apbsInteraction: the molecule does not have at leat 15A gap on each side: increase the boxSize (current box %s)' % str(boxSize))

    ligand = [a for a in al if a.rNo in ligandResidueList]
    protein = [a for a in al if a.rNo not in ligandResidueList]
    if len(ligand)==0 or len(protein)==0:
        sys.exit('Error in apbsInteraction: invalid selection of the ligand residues: ligand has %d atoms out of %d' % (len(ligand),len(al)))

    # where should i store temporary files?
    if os.path.exists('/temp1'):
        try:
            getlogin = os.getenv('LOGNAME')            
            tempdirname = '/temp1/'+getlogin
            if not os.path.exists(tempdirname):
                os.mkdir(tempdirname)
            tmpdir = tempfile.mkdtemp(prefix='tmpApbs',dir=tempdirname)
        except:
            tmpdir = tempfile.mkdtemp(prefix='tmpApbs',dir='.')
    else:
        tmpdir = tempfile.mkdtemp(prefix='tmpApbs',dir='.')
    cwd = os.getcwd()
    os.chdir(tmpdir)

    fcom = 'complex.pqr'
    fpro = 'protein.pqr'
    flig = 'ligand.pqr'
    fin  = 'command.in'

    pqr.write(fcom,al, amberRadii=True)
    pqr.write(fpro,protein, amberRadii=True)
    pqr.write(flig,ligand, amberRadii=True)
    if debug:
        bgf.write('complex.bgf',al)
        bgf.write('protein.bgf',protein)
        bgf.write('ligand.bgf',ligand)

    # prepare the first input file -- to outpur dx files
    with open(fin,'w') as f:
        # header
        f.write(stringHeader % (fcom, fpro, flig)) 
        # body
        entries = ['complex','protein','ligand']
        for i in range(len(entries)):
            ei = entries[i]
            el = ei[0]
            f.write(stringBodyMembranePrep % (ei,'solv',gridSize[0],gridSize[1],gridSize[2],boxSize[0],boxSize[1],boxSize[2],i+1,80.,el,el,el,el) )
        # end
        f.write('quit\n')

    # run apbs - dummy run to prepare dx files
    subprocess.call('apbs '+fin, shell=True)  

    # fix the dx files -- i.e. insert the membrane
    for entry in entries:
        el = entry[0]
        subprocess.call('apbs_draw_membrane dielx_%s.dx 4.0 9.0 9.0' % el, shell=True)  

    # prepare new input file 
    with open(fin,'w') as f:
        # header
        f.write(stringHeaderMembrane % (fcom, fpro, flig)) 
        # body
        for i in range(len(entries)):
            ei = entries[i]
            f.write(stringBodyMembrane % (ei,'solv',gridSize[0],gridSize[1],gridSize[2],boxSize[0],boxSize[1],boxSize[2],i+1,80.,i+1,i+1) )
            f.write(stringBody % (ei, 'ref',gridSize[0],gridSize[1],gridSize[2],boxSize[0],boxSize[1],boxSize[2],i+1,refDielectric) )
        # end
        f.write(stringEndSeparate)

    # run apbs actual run
    apbsrun = subprocess.Popen(['apbs', fin], stdout=subprocess.PIPE)
    stdoutstring = apbsrun.communicate()[0]

    dGsolv = 0.
    dGcoul = 0.
    dGbind = 0.
    dGligsolv = 0.

    # parse the output
    lines = stdoutstring.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        if line=='print energy 1 (complex_solv) - 2 (complex_ref) - 3 (protein_solv) + 4 (protein_ref) - 5 (ligand_solv) + 6 (ligand_ref) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGsolv = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 2 (complex_ref) - 4 (protein_ref) - 6 (ligand_ref) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGcoul = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 1 (complex_solv) - 3 (protein_solv) - 5 (ligand_solv) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGbind = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 1 (complex_solv) - 2 (protein_solv) - 3 (ligand_solv) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGbind = units.joule2calorie( float(s[-2]) )
        elif line=='print energy 5 (ligand_solv) - 6 (ligand_ref) end':
            i+=1
            line = lines[i]
            s = line.split()
            dGligsolv = units.joule2calorie( float(s[-2]) )
        i+=1 

    os.chdir(cwd)
    if not debug:
        shutil.rmtree(tmpdir)
    return {'dGsolv':dGsolv, 'dGcoul':dGcoul, 'dGbind':dGbind, 'dGligsolv': dGligsolv}




