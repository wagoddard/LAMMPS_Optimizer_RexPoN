""" unit conversion for convenience """

def hartree2kcalmol(energy):
    return energy * 627.509469  # convert from hartrees to kcal/mol

def kcalmol2hartree(energy):
    return energy / 627.509469 

def joule2calorie(energy):
    return energy/4.184

def calorioe2joule(energy):
    return energy*4.184

def bohr2angstrom(x):
    return x * 0.52917721092

def angstrom2bohr(x):
    return x / 0.52917721092

def rydberg2ev(x):
    return x * 13.60569253

def ev2rydberg(x):
    return x / 13.60569253

    




