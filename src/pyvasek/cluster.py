""" clustering algorithms based on Ravi's clustering algorithm new_cv7
Note: the algorithm is worked out for clustering docking poses
"""

import sys
import numpy as np

def rmsdMatrixAndList(manyxyz, align = False):
    if align:
        sys.exit('Error: cluster.rmsdMatrix: align = True is not implemented yet')
    n = manyxyz.shape[0]
    natoms = manyxyz.shape[1]
    oversqrtn = 1./np.sqrt(natoms)
    rmsd = np.zeros( (n,n) )
    rmsdList = []
    for i in range(n):
        for j in range(i+1,n):
            r = np.linalg.norm( manyxyz[i] - manyxyz[j]) * oversqrtn
            rmsd[i,j] = r
            rmsd[j,i] = r
            rmsdList.append( (r,i,j) )
    return rmsd, rmsdList

def cluster(manyxyz, diversity = 2.0, align = False):
    n = manyxyz.shape[0]

    # 1. compute the rmsd between all elements
    rmsdMatrix, rmsdList = rmsdMatrixAndList(manyxyz, align = align)
    rmsdList.sort()

    # 2. go throght the pairs which are close together to possibly start new families
    itemToFamily = [-1]*n # -1 means unassigned
    families = {}         # dictionary[family index] -> list of members 
    nfamilies = 0

    for r,i,j in rmsdList:
        if r <= diversity:
            if itemToFamily[i] < 0:
                # i is not in any family
                if itemToFamily[j] < 0:
                    # j is also not in any family
                    # --> create new family for i,j
                    newfamily = nfamilies
                    nfamilies += 1
                    families[newfamily] = [i,j]
                    itemToFamily[i] = newfamily
                    itemToFamily[j] = newfamily
                else:
                    # j is already in some family
                    # --> check if we can add i
                    family = families[ itemToFamily[j] ]
                    canAddI = True
                    for k in family:
                        if rmsdMatrix[i,k] > diversity:
                            canAddI = False
                            break
                    if canAddI:
                        family.append( i )
                        itemToFamily[i] = itemToFamily[j]
            else: 
                # i is in some family
                if itemToFamily[j] < 0:
                    # but j is not in any family
                    # --> check if we can add j
                    family = families[ itemToFamily[i] ]
                    canAddJ = True
                    for k in family:
                        if rmsdMatrix[j,k] > diversity:
                            canAddJ = False
                            break
                    if canAddJ:
                        family.append( j )
                        itemToFamily[j] = itemToFamily[i]
                else:
                    # j is also in some family
                    # --> check if the two families can be merged
                    ifamily = itemToFamily[i]
                    jfamily = itemToFamily[j]
                    if ifamily == jfamily:
                        # i and j already belong to the same family
                        pass
                    else:
                        # i and j belong to two different families
                        rmsdSelection = rmsdMatrix[ [[f] for f in families[ifamily]], families[jfamily]]
                        if np.max(rmsdSelection) > diversity:
                            # --> cannot merge the two families
                            pass
                        else:
                            # --> we can merge the two families
                            # update the family info for the jfamily
                            for k in families[ jfamily ]:
                                itemToFamily[k] = ifamily 
                            # merge the families
                            families[ifamily].extend( families[jfamily] )
                            del families[jfamily]

    # 3. select family heads (each existing family has at least two members now)
    familyHeads = {}
    for indf in families:
        family = families[indf]
        if len(family) <= 2:
            # for families of size two, choose the family head randomly
            familyHead = family[0]
        else:
            #print
            #print 'FAMILY family', family
            rmsdSelection = rmsdMatrix[[[f] for f in family], family]
            #print rmsdSelection
            #print np.sum(rmsdSelection,axis=0)
            #print np.argmin( np.sum(rmsdSelection,axis=0) ) 
            familyHead = family[  np.argmin( np.sum(rmsdSelection,axis=0) ) ]
            #print familyHead
    
        family.sort()
        familyHeads[familyHead] = family

    # 4. the items, that were not assigned yet, make into families with one member
    for k in range(n):
        if itemToFamily[k] < 0:
            familyHeads[k] = [k]

    #for k in familyHeads:
    #    print 'analyzing family with head  %d'%k
    #    family = familyHeads[k]
    #    for i in family:
    #        for j in family:
    #            r=rmsdMatrix[i,j]
    #            if r>diversity:
    #                print '??? %d - %d is %f' % (i,j,diversity)

    return familyHeads


