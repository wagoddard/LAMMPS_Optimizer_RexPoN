import os, numpy, sys

rootdir = os.getenv('VASEK')
#rootdir = os.path.dirname(os.path.dirname( __file__ ))

if rootdir == None:
    print >> sys.stderr, "Error: Environmental variable VASEK is not set (I don't know where to look for files)."
    sys.exit(1)

amberDat = rootdir + '/amber/parm99sb.dat'
amberDatGaff = rootdir + '/amber/gaff.dat'

dreiding3 = rootdir + '/dreiding/dreiding304.shortpar'
dreiding4 = rootdir + '/dreiding/dreiding402.shortpar'

dreiding3cnv = rootdir + '/dreiding/dreiding304.cnv'


waterBox = rootdir + '/water/water.bgf'
waterBoxSize= numpy.array([[30.0,0.0,0.0],\
                           [0.0,30.0,0.0],\
                           [0.0,0.0,30.0],\
                           [0.0,0.0,0.0]]) # has to be orthogonal for making large box

# reference amino acids
refBgf = {}
refBgf['dreiding3-charmm22'] = rootdir + '/amino/dreiding3-charmm22.bgf'
refBgf['dreiding3-ff99SB'] = rootdir + '/amino/dreiding3-ff99SB.bgf'
refBgf['dreiding4-charmm22'] = rootdir + '/amino/dreiding4-charmm22.bgf'
refBgf['dreiding4-ff99SB'] = rootdir + '/amino/dreiding4-ff99SB.bgf'
refBgf['amber-ff99SB'] = rootdir + '/amino/amber-ff99SB.bgf'

