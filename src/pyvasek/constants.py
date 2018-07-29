""" handy physical constatns and conversions """

import math

permitivityOfVacuum = 8.854187817620e-12
eps0 = permitivityOfVacuum

electronCharge = 1.60217646e-19
e = electronCharge

boltzmannConstant = 1.3806503e-23
kb = boltzmannConstant

avogadrosNumber = 6.02214179e23
Na = avogadrosNumber


kcal = 4184 # joules
kcalPerMol = kcal / Na

coulombFactorKcalMol = 1. / (4. * math.pi * eps0) * e * e * 1e10 / kcalPerMol






