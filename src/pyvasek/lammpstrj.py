""" reading lammps trj file
"""

import sys,numpy
import atoms,bgf

def read(filename):

  lines = open(filename).readlines()

  # timesteps
  nlines = len(lines)
  timesteps = [int(lines[i+1]) for i in range(nlines-1) if lines[i].startswith('ITEM: TIMESTEP')] 
  ntimes = len( timesteps )

  # how much data do we have
  natoms = int(lines[3])
  items = lines[8].split()
  items = items[2:]
  nitems = len(items)

  # data container
  c = numpy.zeros((ntimes,natoms,nitems))

  itime = 0
  iatom = 0
  readingatoms = False

  for i in range(1,nlines):
    prevline = lines[i-1]
    line = lines[i]
    if prevline.startswith('ITEM: TIMESTEP'):
      itime = int(line)
      indextime = timesteps.index(itime)
    elif prevline.startswith('ITEM: ATOMS'):
      readingatoms = True
      
    if readingatoms:
      #print 'read',itime,iatom,line.strip()
      cdata = numpy.array([float(item) for item in line.split()])    
      c[indextime,iatom,:] = cdata
      #print c[indextime,iatom,:]
      #print cdata
      iatom +=1 

    if iatom == natoms:
      iatom = 0
      readingatoms = False

  return items, timesteps, c

def readQ(filename):
  items, timesteps, c = read(filename)
  if not 'q' in items:
    sys.error('Error: q is not column in the lammps trajectory %s'%filename)

  index = items.index('q')
  cout = c[:,:,index]
  return cout

def readXYZ(filename):
  items, timesteps, c = read(filename)
  if (not 'x' in items) or (not 'y' in items) or (not 'z' in items):
    sys.error('Error: x,y,z are not columns in the lammps trajectory %s (only %s)'%(filename,items))

  indices = [items.index('x'),items.index('y'),items.index('z')]
  cout = c[:,:,indices]
  return cout


