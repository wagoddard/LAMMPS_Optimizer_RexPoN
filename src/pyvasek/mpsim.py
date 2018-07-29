'''provides various calls to mpsim'''

import os, subprocess, sys, tempfile, shutil
from pyvasek import bgf

def sgb(al, zeroCharge = True, debug = False):
    '''Call mpsim to evaluate SGB energy
    if zeroCharge = True
        charge of all the atoms in the atom list al is zeroed
        and so the energy is the nonpolar solvation energy
        (which is independent of dielectric, the polar one depends on dielectric)
    '''
    #paramfiles = ['sgb.prr.str.prm.real','sgb.sncorrfnprm.noself','sgb.sncorrfnprm.self','sgb.slr_param','sgb.param']
    #for f in paramfiles:
    #    if os.path.islink(f):
    #        os.unlink(f)

    fcom = tempfile.NamedTemporaryFile(mode='w', prefix='tmpMpsim', suffix='.bgf', dir='.', delete=False)
    fout = tempfile.NamedTemporaryFile(mode='w', prefix='tmpMpsim', suffix='.out', dir='.', delete=False)
    fcom.close()
    fcom.close()    
    fcom = fcom.name
    fout = fout.name
    
    if zeroCharge:
        for a in al:
            a.charge = 0.

    bgf.write(fcom,al)
    subprocess.call('mpsim.pl -f %s -s SGB -a ONEE -p /project/Biogroup/FF/dreiding-0.3LJ.par | grep SGBSolv > %s' % (fcom,fout), shell=True) 

    lines = open(fout).readlines()
    if len(lines) < 1:
        sys.exit('Error in mpsim.sgb: mpsim did not return valid output')
    
    en = float( lines[0].split(':')[-1] )

    if not debug:
        os.remove(fcom)
        os.remove(fout)
        root = fcom.split('.')[0]
        os.remove('%s-done.traj1'%root)
        #for f in paramfiles:
        #    if os.path.islink(f):
        #        os.unlink(f)

    return en

def sgbInteraction(al, residueList = [999], zeroCharge = True):
    lig = [a for a in al if a.rNo in residueList]
    prot = [a for a in al if a.rNo not in residueList]
    return sgb(al,zeroCharge) - sgb(prot,zeroCharge) - sgb(lig,zeroCharge)


