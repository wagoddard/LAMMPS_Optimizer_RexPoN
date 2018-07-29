#!/usr/bin/env python
'''
Program:              lammpsEnergy.py
Sample command when you have multiple bgf files with exactly same system differing only in
    coordinates:      lammpsEnergy.py -c topology.bgf -n 10 file1.bgf file2.bgf protein*.bgf 
Sample command when you have one bgf file and one xyz file with multiple sets of coordinates:      
                      lammpsEnergy.py -c topology.bgf -f amber -n 10 coordinates.xyz

Description:          Runs minimization in lammps. 

Required inputs: 
  bgf mode: bgf files with the same connectivity to have energy evaluated
    if the topology file is not specified, the first bgf file is used to read 
    connectivitym, charges and force field types
  xyz mode: one xyz file with one or more sets of coordinates and one bgf file with topology

Output:
  if number of minimization steps > 0, the minimized geometry is stored with .m.bgf extension

Optional options:
  -c --topology       topology file (.bgf or .mol2) which is used for reading connectivity
                      charges, and atom types, but not coordinates
  -n --nsteps         number of minimization steps, 
                      default is 0, i.e. just evaluate energy
  -t --threshold      force treshold for minimization; if threshold is met, the
                      minimization stops before reaching min iterations
                      default 0.0
  -s --solvation      yes/no NOT IMPLEMENTED YET
  -f --forcefield     'dreiding3' (default), 'dreiding4' or 'amber'
  -p --parameterfile  in vaclav's "shortpar" format for dreiding
                      or corrected dat format for amber 
                      for default see $VASEK/pyvasek/defaults.py
  -k --keep           keep lammps input files
  -d --debug          print debug info if the option is present
  -h --help

Author:               Vaclav Cvicek
Date:                 4/2011.
'''

import os, sys, string, getopt
from pyvasek import bgf,energy, atoms,defaults,xyzfile,mol2

def usage():
    print __doc__

def printHeader(doPrint = True):
    s = '%15s,' % 'filename'    
    if doPrint:
        for name in energy.names: 
            s += '%10s,' %name
        print s
    else:
        for name in energy.names: 
            s += '%15s,' %name
        return s

def printEnergy(e,fname,doPrint = True):
    s = '%15s,' % fname 
    if doPrint:
        for name in energy.names:
            s += '%10.4f,' % e[energy.kind[name]]
        print s
    else:
        for name in energy.names:
            s += '%15.9f,' % e[energy.kind[name]]
        return s

def addMinToName(bgfname):
    s = bgfname.split('.')
    s.insert(-1,'m')
    return '.'.join(s)
    
def uniRead(filename, debug = False):
    ext = filename.split('.')[-1]
    if ext == 'mol2':
        return mol2.read(filename)    
    else:
        return bgf.read(filename, debug = debug)

def uniWrite(filename,al):
    ext = filename.split('.')[-1]
    if ext == 'mol2':
        mol2.write(filename,al)    
    else:
        bgf.write(filename,al)    

#################### parse input ##############################
def main(sysargv, doPrint = True):  

    # defaults
    topology = ''
    nsteps = 0
    threshold = 0.0
    force = 'dreiding3'
    parameterfile = ''
    debug = False
    keep = False

    # return value
    myreturn = []

    try:
        # options requiring argument have ':' in short string, '=' in long string
        # args contains unparsed options

        #opts, args = getopt.getopt(sys.argv[1:], "c:n:t:f:p:dh", 
        opts, args = getopt.getopt(sysargv, "c:n:t:f:p:dkh", ["topology=","nsteps=","threshold=","forcefield=","parameterfile=",'keep',"debug","help"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-c", "--topology"):
            topology = a
        elif o in ("-n", "--nsteps"):
            nsteps = int(a)
        elif o in ("-t", "--threshold"):
            threshold = float(a)
        elif o in ("-f", "--forcefield"):
            force = a
        elif o in ("-p", "--parameterfile"):
            parameterfile = a
        elif o in ("-k", "--keep"):
            keep = True
        elif o in ("-d", "--debug"):
            debug = True
        else:
            assert False, "unhandled option"
            usage()
            sys.exit()
    # decide if running in bgf mode or xyz mode
    if len(args) == 0:
        print 'ERROR: INPUT REQUIRED!'        
        usage()
        sys.exit()
    elif args[0][-4:] == '.xyz':
        if len(args) != 1:
            print 'ERROR: only one xyz file is read when coordinates are in xyz file, but here there are %d input items on the comand line. ' % len(args)
            sys.exit()            
        inputIsBgf = False
        xyzlist,namelist = xyzfile.read(args[0])
    elif args[0][-4:] == '.bgf':
        inputIsBgf = True
        namelist = args
    elif args[0].split('.')[-1] == 'mol2':
        inputIsBgf = True
        namelist = args
    else:
        print 'ERROR: input file is neither bgf nor mol2 nor xyz: %s ' % args[0]
        sys.exit()            
    
    # finish setting up
    if topology == '':
        topology = namelist[0]

    if force=='dreiding3':
        myforcefield = 'dreiding'
        defaultparameter = defaults.dreiding3
    elif force=='dreiding4':
        myforcefield = 'dreiding'
        defaultparameter = defaults.dreiding4
    elif force=='amber':
        myforcefield = 'amber'
        # if passing empty string, both Amber and Gaff will get loaded in energy.py
        # from the files specified in defaults.py 
        defaultparameter = '' 
    else:
        print 'ERROR: unknown forcefield type: %s:' %force        
        sys.exit()

    if parameterfile == '':
        parameterfile = defaultparameter

    # initialization
    en = energy.Energy(type=myforcefield, parameterFile=parameterfile , debug=debug)
    al = uniRead(topology,debug=debug)

    if inputIsBgf:        
        if topology != namelist[0]:
            file0 = uniRead(namelist[0])
            xyz = atoms.getXyzFromList( file0)
            atoms.putXyzIntoList(al,xyz)
    else:
        xyz = xyzlist[0,:,:]
        atoms.putXyzIntoList(al,xyz)
    nameprefix = en.load(al)
    numAtoms = len(al) # number of atoms

    # first evaluation
    if doPrint:
        printHeader()
    else:
        myreturn.append( printHeader(doPrint = doPrint) )

    index = 1
    e = en.minimize(steps = nsteps, forceTolerance = threshold)
    if doPrint:
        printEnergy(e,namelist[0])
    else:
        myreturn.append( printEnergy(e,namelist[0],doPrint=doPrint) )

    # if nsteps is larger than 0 also store the minimized geometry
    if nsteps > 0:
        xyz = en.getXyz()
        if inputIsBgf:        
            atoms.putXyzIntoList(al,xyz)
            uniWrite(addMinToName(topology),al)
        else:
            xyzlist[0,:,:] = xyz.copy()

    # evaluations
    for i in range(1,len(namelist)):
        a = namelist[i]
        if inputIsBgf:
            al = uniRead(a,debug=debug)
            if len(al) != numAtoms:
                print 'ERROR: file %s has different number of atoms (%d) than the topology file (%d)' % (a,len(al),numAtoms)        
                sys.exit()            
            xyz = atoms.getXyzFromList( al )
        else:
            xyz = xyzlist[i,:,:]
        e = en.eval(xyz,minsteps = nsteps,forceTolerance = threshold)

        if doPrint:
            printEnergy(e,a)
        else:
            myreturn.append( printEnergy(e,a,doPrint=doPrint) )

        if nsteps > 0:
            xyz = en.getXyz()
            if inputIsBgf:        
                atoms.putXyzIntoList(al,xyz)
                uniWrite(addMinToName(a),al)
            else:
                xyzlist[i,:,:] = xyz.copy()
    # in the bgf mode, all the minimized geometries were already written out, but not the xyz mode
    if not inputIsBgf:
        filename = args[0][:-4]+'.m.xyz'
        xyzfile.write(filename,xyzlist,names=namelist,al=al)

    if not keep:
        os.remove(nameprefix+'.in')
        os.remove(nameprefix+'.data')
        os.remove(nameprefix+'.log')
            
    return myreturn

#####################################
##################################
#if __name__ == "__main__":
#    main()


