
import sys, numpy, os, time, ctypes, tempfile
import amberFF,shortpar,prepLammps, forceField,defaults
from packages import lammps

names = ['step','ebond','eangle','edihed','eimp','emol','ecoul','evdwl','hbond','counthbond','elong','epair','pe','ke','etotal','temp','press','fmax']
kind = {}
for i in range(len(names)): kind[names[i]] = i

class Energy:
    def __init__(self,type='dreiding',parameterFile='', debug=False, providedForceField = None):
        '''forceField.ForceField instance can be provided as parameter providedForceField or
        it will be read from the file parameterFile'''

        self.debug = debug
        self.lmp = None # lammps engine
        self.structure = None
        
        if providedForceField == None:
            self.ff = forceField.ForceField()
            if type == 'dreiding':
                shortpar.read(self.ff,fileName=parameterFile,debug=self.debug)
            elif type == 'amber':
                if parameterFile == '':
                    amberFF.readAmberDat(self.ff,fileName=defaults.amberDat,debug=self.debug)
                    amberFF.readAmberDat(self.ff,fileName=defaults.amberDatGaff,debug=self.debug)
                else:
                    amberFF.readAmberDat(self.ff,fileName=parameterFile,debug=self.debug)
            else:
                print >> sys.stderr, 'unknown force field  %s  (only dreiding or amber)' % type
                sys.exit(1)
        else:
            self.ff = providedForceField

    def load(self,al,cell=None):
        ''' temporary lammps data and input files are written for the structure in al
        by default new lammps.lammps() instance is created, but an older lammps instance can be reused 
        when passed by the lmp parameter, in this old instance the clear command is issued to allow using of a new structure
        returns the prefix of the temporary names
        '''

        f = tempfile.NamedTemporaryFile(mode='w', prefix='tmpLammps', suffix='.in', dir='.', delete=False)
        f.close()
        filename = f.name[:-3]
        self.structure = prepLammps.generateLammpsData(al, self.ff, filename, cell=cell, debug=self.debug)
        
        if self.debug:
            self.lmp = lammps.lammps(['-echo','screen','-log',filename+'.log'])
        else:
            self.lmp = lammps.lammps(['-screen','none','-log',filename+'.log'])
        
        self.lmp.file(filename+'.in')
        return filename
         
    def eval(self,xyz,minsteps = 0,forceTolerance = 0.0):
        self.setXyz(xyz)
        return self.minimize(steps = minsteps, forceTolerance=forceTolerance)
        
    def get(self):
        energy = numpy.zeros(len(names))
        for i,name in enumerate(names):
            #try:
            #    compute = self.lmp.extract_compute(group+name,0,0)
            #except ValueError:
            #    compute = 0.0 # is the energy term is not defined, let's say it is 0
            # i have redefined the variables names, so that they have to be present
            compute = self.lmp.extract_variable(name,0,0)
            energy[i] = compute
        if self.debug:
            print energy
        return energy 

    def getXyz(self):
        xyz_ctypes = self.lmp.get_coords()
        xyz2 = numpy.frombuffer( buffer(xyz_ctypes) , float) 
                        # there is an extra step to wrap the ctypes array in buffer(), to avoid
                        # some warning numpy() because ctypes produces wrong PEP 3118 code 
        xyz3 = xyz2.reshape(len(xyz2)/3, 3)
        return xyz3
    
    def setXyz(self,xyz):
        self.lmp.command('reset_timestep   0')
        self.lmp.put_coords(xyz.ctypes.data_as(ctypes.c_void_p))
        
    def minimize(self,steps = 0,forceTolerance = 0.0):
        # minimize etol ftol maxiter maxeval
        self.lmp.command('minimize               0.0 %f %d 99999999999' % (forceTolerance, steps))
        return self.get()

def printE(e):
    for name in names:
        print '%s: %f' %(name,e[kind[name]])
        
