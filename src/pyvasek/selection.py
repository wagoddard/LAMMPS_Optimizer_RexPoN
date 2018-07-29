'''
simple selection language for specifying residue and chain:
  example:    134_A,135B,137-139_C
  note: 
    1. middle underscore is not required in 134A but is required in 134_1
    2. chain can be . to mean any chain (default if chain left out)
    3. residues can be range:   134-148
       or list specified by +:  101+103+134-148+150
    4. * means any residue number (default if residues left out)
    5. multiple items are separated with comma
'''

import sys

def residueListFromString(residueString):
    reslist = []
    if residueString == '*':
        reslist.append(-1)
    else:
        resplus = residueString.split('+')
        for resplusi in resplus:
            resminus = resplusi.split('-')
            if len(resminus) == 2:
                istart = int(resminus[0])
                iend = int(resminus[1])+1
                for i in range(istart,iend):
                    reslist.append(i)
            elif len(resminus) == 1:
                reslist.append(int(resminus[0]))
            else:
                sys.exit('Error in selection.residueListFromString() when parsing residue and chain selection rule:%s' %resandchain )
    return reslist

def residueChainSplit(resandchain):
    mystring = resandchain.strip()
    mylist = mystring.split('_')
    if len(mylist) == 2:
        reslist = residueListFromString(mylist[0])
        chain = mylist[1]
    elif len(mylist) == 1:
        if mylist[0] == '*':
            reslist = [-1]
            chain = '.'            
        elif not mylist[0][-1].isdigit():
            reslist = residueListFromString(mylist[0][0:-1])
            chain = mylist[0][-1]
        else:
            # chain was not provided
            reslist = residueListFromString(mylist[0])
            chain = '.'
    else:    
        #defaults:
        reslist = ['-1'] # any residue
        chain = '.' # any chain
    myresids = ['%s_%d'%(chain,r) for r in reslist]
    return myresids
        
def getResidueIDs(selectionString):
    resid = []
    commasplit = selectionString.split(',')
    for resandchain in commasplit:
        resid.extend(residueChainSplit(resandchain))
    resid = list(set(resid))
    resid.sort()
    return resid

def buildResDictionary(al):
    resdict = {}
    for a in al:
        resid = '%s_%d' % (a.chain,a.rNo)
        if resid in resdict:
            resdict[resid].append(a)
        else:
            resdict[resid] = [a]
    return resdict

def get(al,selectionString,separate=False):
    '''returns selection matching the selection rule
    if separate = True, each residue of selection is returned as a separate atom list
    '''
    resids = getResidueIDs(selectionString)
    myresids = buildResDictionary(al)
    sel = []
    for myresid in myresids:
        s = myresid.split('_')
        mychain = s[0]+'_-1'
        myres   = '._'+s[1]    
        myboth   = '._-1'
        if (myboth in resids) or (myresid in resids) or (mychain in resids) or (myres in resids):
            if separate:
                sel.append(myresids[myresid])
            else:
                sel.extend(myresids[myresid])  
    return sel
 
