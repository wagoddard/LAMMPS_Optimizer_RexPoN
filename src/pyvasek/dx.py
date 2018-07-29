""" reading and writing opendx file as written by APBS

multigrid data format:

   # Comments
   object 1 class gridpositions counts nx ny nz
   origin xmin ymin zmin
   delta hx 0.0 0.0
   delta 0.0 hy 0.0 
   delta 0.0 0.0 hz
   object 2 class gridconnections counts nx ny nz
   object 3 class array type double rank 0 items n data follows
   u(0,0,0) u(0,0,1) u(0,0,2)
   ...
   u(0,0,nz-3) u(0,0,nz-2) u(0,0,nz-1)
   u(0,1,0) u(0,1,1) u(0,1,2)
   ...
   u(0,1,nz-3) u(0,1,nz-2) u(0,1,nz-1)
   ...
   u(0,ny-1,nz-3) u(0,ny-1,nz-2) u(0,ny-1,nz-1)
   u(1,0,0) u(1,0,1) u(1,0,2)
   ...
   attribute "dep" string "positions"
   object "regular positions regular connections" class field
   component "positions" value 1
   component "connections" value 2
   component "data" value 3

"""

import sys,numpy

class DxData:
    ''' entries:
        n        the number of points 
        origin   contains the origin
        h        is the diagonal of one unit 
        v        values 
        use gridx (resp y/z) to get the grid points  
        '''
    def __init__(self):
        self.n = numpy.array([0,0,0])
        self.origin = numpy.array([0.,0.,0.])
        self.h = numpy.array([0.,0.,0.])
        self.v = numpy.zeros((0,0,0))

    def grid(self, index):
        return numpy.linspace(self.origin[index],self.origin[index]+(self.n[index]-1)*self.h[index],num=self.n[index])

    def gridx(self):
        return self.grid(0)

    def gridy(self):
        return self.grid(1)

    def gridz(self):
        return self.grid(2)

    def potToVolt(self, temperature = 298.15):
        coeff = 1.3806504e-23 * temperature / 1.60217646e-19
        self.v *= coeff

def read(filename):
    dxdata = DxData()
    f = open(filename)

    # comments
    line = f.readline()
    while line.startswith('#'):
        line = f.readline()
    
    # grid n
    s = line.split()
    dxdata.n[2] = int(s[-1])     
    dxdata.n[1] = int(s[-2])            
    dxdata.n[0] = int(s[-3])

    # origin
    line = f.readline()
    s = line.split()
    dxdata.origin[0] = float(s[1])     
    dxdata.origin[1] = float(s[2])            
    dxdata.origin[2] = float(s[3])

    # h
    line = f.readline()
    dxdata.h[0] = float(line.split()[1]) 
    line = f.readline()
    dxdata.h[1] = float(line.split()[2]) 
    line = f.readline()
    dxdata.h[2] = float(line.split()[3]) 

    line = f.readline()
    line = f.readline()

    # reading data
    datapoints = dxdata.n[0]*dxdata.n[1]*dxdata.n[2]
    dxdata.v = numpy.zeros(datapoints)
    i = 0
    while i<datapoints:
        s = f.readline().split()
        for item in s:
            dxdata.v[i] = float(item)
            i += 1
    f.close()
    dxdata.v = dxdata.v.reshape((dxdata.n[0],dxdata.n[1],dxdata.n[2]))
    return dxdata

def write(filename,dxdata):
    f = open(filename,'w')
    f.write('# written by dx.py by vcvicek \n')
    f.write('# \n')
    f.write('object 1 class gridpositions counts %d %d %d\n'%(dxdata.n[0],dxdata.n[1],dxdata.n[2]))
    f.write('origin %12e %12e %12e\n'%(dxdata.origin[0],dxdata.origin[1],dxdata.origin[2]))
    f.write('delta %12e 0.000000e+00 0.000000e+00\n'%dxdata.h[0])
    f.write('delta 0.000000e+00 %12e 0.000000e+00\n'%dxdata.h[1])
    f.write('delta 0.000000e+00 0.000000e+00 %12e\n'%dxdata.h[2])
    f.write('object 2 class gridconnections counts %d %d %d\n'%(dxdata.n[0],dxdata.n[1],dxdata.n[2]))
    f.write('object 3 class array type double rank 0 items %d data follows\n' % (dxdata.n[0]*dxdata.n[1]*dxdata.n[2]))
    counter = 0
    for i in range(dxdata.n[0]):
        for j in range(dxdata.n[1]):
            for k in range(dxdata.n[2]):
                f.write('%12e '%dxdata.v[i,j,k])
                counter += 1
                if counter % 3 == 0:
                    f.write('\n')
            
    if counter % 3 != 0:
        f.write('\n')

    f.write('attribute "dep" string "positions"\n')
    f.write('object "regular positions regular connections" class field\n')
    f.write('component "positions" value 1\n')
    f.write('component "connections" value 2\n')
    f.write('component "data" value 3\n')

    f.close()




