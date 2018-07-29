import time, sys
from packages.pizza import vmd as pvmd
import prmtop, bgf, amberFF, forceField, mol2, atoms, energy

class vmd:
    def __init__(self):
        self.pvmd = pvmd.vmd()
        self.pvmd.VMD.delaybeforesend = 0.0
    
    def show(self, al):
        bgf.writeFile('tempvmd.bgf',al)
        self.pvmd.new('tempvmd.bgf',filetype='bgf')
        
    def setAtoms(self,al,append=True):
        xyz = atoms.getXyzFromList(al)
        self.set(xyz,append=append)
        
    def set(self,xyz,append=True):
        if append:
            self.pvmd.__call__('animate dup [molinfo top]')
        # need to break up the data into pieces around 100 atoms each because otherwise the pexpect module stalls
        # pieces of 1000, and 100
        cutlarge = 1000
        cut = 100
        if len(xyz) % cutlarge == 0: 
            nlarge = len(xyz) / cutlarge
        else:
            nlarge = len(xyz) / cutlarge + 1
        for ilarge in range(nlarge):
            firstlarge = ilarge*cutlarge
            lastlarge = min((ilarge+1)*cutlarge, len(xyz))
            
            self.pvmd.__call__('set vmdsel [atomselect top "index >= %d and index < %d"]' % (firstlarge,lastlarge))
            self.pvmd.__call__nowait('$vmdsel set {x y z} {')  # notice: there is not flushing of input here            
            
            if (lastlarge-firstlarge) % cut == 0: 
                n = (lastlarge-firstlarge) / cut
            else:
                n = (lastlarge-firstlarge) / cut + 1
            for i in range(n):
                sys.stdout.write('.')
                #print '....vmd: %d out of %d' % (i,n) 
                first = firstlarge+i*cut
                last = min(firstlarge+(i+1)*cut, len(xyz))
                # now actually updating the data
                cmd = ''
                for ii in range(first,last):
                    cmd += ' { %.2f %.2f %.2f }' %(xyz[ii,0],xyz[ii,1],xyz[ii,2])
                #print first
                #print last
                self.pvmd.__call__nowait(cmd)
            print '.'
            self.pvmd.__call__(' }')            
        self.pvmd.__call__('$vmdsel delete ; unset vmdsel')
        self.pvmd.flush()        

    def minimize(self,en,steps,plot=None):
        step = steps/10
        # plot decreasing energy too
        energies = []
        for i in range(step):
            print 'starting 10 steps of minimization minimization, i = %d' %i
            e = en.minimize(steps = 10)
            xyz = en.getXyz()
            self.set( xyz )
            # plot energy
            energies.append(e[energy.kind['total']])
            if plot is not None:
                plot.plot(energies)    
            print 'finished 10 steps of minimization minimization'
