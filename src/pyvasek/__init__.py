'''
vcvicek
'''


''' help to make the module list:
import os
rootdir = os.getenv('VASEK')
list1 = os.listdir(rootdir+'/pyvasek')
list2 = []
for a in list1:
    if a[-3:]=='.py':
        list2.append(a[:-3])
'''
__all__ = ['amberFF', 'atoms',
           'bgf','bonds',
           'constants', 
           'defaults','dx','dreiding3','dreiding4',
           'energy',
           'geometry','gro',
           'fasta','forceField', 
           'jaguar',
           'lammpstrj',
           'mol2','mycsv', 
           'pdb', 'periodic','pqr','prepLammps','prmtop', 'psf',
           'qchem',
           'residues',
           'shortpar','solvate','structure',
           'xst','xyzfile']
