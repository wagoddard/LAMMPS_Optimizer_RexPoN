""" tools for manipulatin .sph files with spheres from sphgen (part of dock 4.0.1)

Info from 
http://www.csb.yale.edu/userguides/datamanip/dock/DOCK_4.0.1/html/Manual.20.html


Some informative messages are written to a file called outsph. This includes the parameters and files used in the calculation. The spheres themselves are written to outfil. They are arranged in clusters with the cluster having the largest number of spheres appearing first. The sphere cluster file consists of a header followed by a series of sphere clusters. The header is the line

DOCK 3.5 receptor_spheres

followed by a color table. The color table contains color names (format A30) each on a separate line. As sphgen produces no colors, the color table is simply absent. The sphere clusters themselves follow, each of which starts with the line

cluster n number of spheres in cluster i

where n is the cluster number for that sphere cluster, and i is the number of spheres in that cluster. Next, all spheres in that cluster are listed in the format

(I5, 3F10.5, F8.3, I5, I2, I3)

where the values correspond to, respectively,

    The number of the atom with which surface point i (used to generate the sphere) is associated.
    The x, y, and z coordinates of the sphere center.
    The sphere radius.
    The number of the atom with which surface point j (second point used to generate the sphere) is associated.
    The critical cluster to which this sphere belongs.
    The sphere color. The color is simply an index into the color table that was specified in the header. So 1 corresponds to the first color in the header, 2 for the second, etc. 0 corresponds to unlabeled.

The clusters are listed in numerical order from largest cluster found to the smallest. At the end of the clusters is cluster number 0. This is not an actual sphere cluster, but a list of all of the spheres generated whose radii were larger than the minimum radius, before the filtering heuristics (i.e. allowing only one sphere per atom and using a maximum radius cutoff) and clustering were performed. Cluster 0 may be useful as a starting point for users who want to explore a wider range of possible clusters than is provided by the standard sphgen clustering routine. The program cluster takes the full sphere description as input, and allows the user to explore different sphere descriptions of the site. This is particularly useful for Macromolecular Docking ; it is often inefficient to use spheres that fill the entire volume of the "ligand" macromolecule. In addition, only a portion of a cavity in the "receptor" macromolecule may be of interest for docking purposes. If the standard clustered output from sphgen provides a satisfactory description of the ligand molecule or receptor site, running cluster is not necessary.

The program creates three temporary files: temp1.ms, temp2.sph, and temp3.atc. These are used internally by sphgen, and are deleted upon completion of the computation.

"""

import numpy
import pdb,residues,atoms

def getXyzFromSphLine(line):
    x = float(line[5:15]) 
    y = float(line[15:25])
    z = float(line[25:35])
    return numpy.array([x,y,z])

def getDistFromRes(r,xyz):
    dist = 0.0
    count = 0.0
    for a in r.atoms:
        count += 1.
        dist += numpy.linalg.norm(xyz-a.xyz)
    return dist/count

def selectSpheresNearResidues(sphereFilename,outSphereFilename,al,nearResidues, distance):
    '''reads the spheres from the file  sphereFilename  and selects the spheres for which 
    sum of distances from residues (CA atom) is less then 
    the distance parameter
    nearResidues is list of residue numbers (integers)
    '''
    
    rl = residues.makeList(al)
    # select residues of interest
    rl = [r for r in rl if r.rNo in nearResidues]
    # removing non heavy atoms
    for r in rl:
        for a in r.atoms[:]:
            if a.aName.strip() != 'CA':
                r.atoms.remove(a)

    # read the spheres file
    lines = open(sphereFilename).readlines()
    outlines = []
    for line in lines[1:]:
        xyz = getXyzFromSphLine(line)
        dist = 0.0
        #outtext = '' 
        for r in rl:
            di = getDistFromRes(r,xyz)
            dist += di
            #outtext += '%f  ' % di
        #outtext += '%f  ' % dist
        #print outtext
        if dist <= distance:
            outlines.append(line)
    nselected = len(outlines)

    print 'Selected %d out of %s spheres (with distance cutoff %f from reisdues %s)' % (nselected,lines[0].split()[-1],distance,str(nearResidues))

    firstline = 'cluster     0   number of spheres in cluster %5d\n' % nselected
    outlines.insert(0,firstline)

    with open(outSphereFilename,'w') as f:
        f.writelines(outlines)

    return nselected

def sphFileToAtoms(sphereFilename):
    lines = open(sphereFilename).readlines()
    al = []
    for line in lines[1:]:
        if line.strip() != '':
            a = atoms.Atom()
            a.aName = 'Sc'
            a.rNo = 1
            a.rName = 'SPH'
            a.chain = 'S'
            a.xyz = getXyzFromSphLine(line)            
            al.append(a)
    return al

