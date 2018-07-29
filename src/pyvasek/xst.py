'''
Formatting of extended trajectory file produced by NAMD:
comments start with #
For the data:

Hi,

Column 1: time step;

Column 2-10: the three unit cell vectors, same meaning as
cellBasisVector[1-3];

Column 11-13: same meaning as cellOrigin;

The last 6 columns are something related to the strain rate of the
unit cell.

Zhu 


Note that the xst file has extra datapoint -- the first line is the
initial size of the system -- i.e. for this one there is no trajectory
snapshot in dcd
'''

import numpy

def readBox(filename):
    '''returns the box sizes from the extended trajectory file xst produced by NAMD

    boxSize is numpy array of shape:

    a_x a_y a_z
    b_x b_y b_z
    c_x c_y c_z
    o_x o_y o_z

    the list is concatenated into extra dimension of the numpy.array:
        numberOfFrames x 4 x 3
    '''
    steps = []
    lines = open(filename).readlines()
    for line in lines:
        s1 = line.strip()
        if s1 != '' and s1[0] != '#':
            s2 = s1.split()
            a_x = float(s2[1])
            a_y = float(s2[2])
            a_z = float(s2[3])
            b_x = float(s2[4])
            b_y = float(s2[5])
            b_z = float(s2[6])
            c_x = float(s2[7])
            c_y = float(s2[8])
            c_z = float(s2[9])
            o_x = float(s2[10])
            o_y = float(s2[11])
            o_z = float(s2[12])
            box = [[a_x,a_y,a_z],[b_x,b_y,b_z],[c_x,c_y,c_z],[o_x,o_y,o_z]]
            steps.append(box)
    return numpy.array( steps )

def readStep(filename):
    '''reads the numbers from the step column in the xst file
    returns list of integers
    '''
    steps = []
    lines = open(filename).readlines()
    for line in lines:
        s1 = line.strip()
        if s1 != '' and s1[0] != '#':
            s2 = s1.split()
            steps.append(int(s2[0]))
    return steps


