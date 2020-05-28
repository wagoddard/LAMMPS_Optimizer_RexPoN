# LAMMPS_Optimizer_RexPoN

================

The LAMMPS Optimizer is able to optimize any parameters used in a LAMMPS calculation to a given training set.

To run the program:
    Prepare "ffield", "param", "trainset", and the necessary data files.
    Run "python main.py"

Features
--------

- Can choose many different optimization algorithms including SLSQP, Nelder-Mead (gradient free), and a mapping algorithm
- Can constrain the system with bounds and variable constraints
- Can optimize any parameters using the UFFQM, LJCUT, OPLS, MORSE, EXRYD, VAPOX forcefield

Getting Started
---------------

### Installation
- TODO LATER
 

### Calling the optimizer
- Using a python package with an updated version of SCIPY, call 'python main.py'
- Optional arguments for directing to proper input files: '-g' (geometry file, if included, will generate data files), '-t' (trainset file; default: trainset), '-p' (params file; default: param), and '-f' (forcefield file; default: ffield)


### ffield Formatting

Sections include 'OPTIMIZER' (required), 'FFIELD' (required), 'ATOMS', 'BONDS', 'ANGLES', 'DIHEDRALS', 'VDW', 'ELECTRO', 'CONTROL'

1. 'OPTIMIZER' contains the information determing what type of optimization to be performed including the method (options: None, Optimization, or Mapping; default: None), polarization flag (extract the energy due to polarization; default: Off), algorithm (which type of SCIPY optimization algorithm to use; options: SLSQP, Nelder-Mead; default: SLSQP), finish (optimization performed after mapping; options: SLSQP, Nelder-Mead, None; default: None), and valence (used for created of data files, default: 0)
2. 'FFIELD' contains one line with the name of the force field that should be used (options: uffqm, vapox, morse, ljcut, exryd, opls)
3. 'ATOMS' contains a list of the atom types used with format '(atom number) (atom mass) # (atom name)'. Note that any text after a '#' is read as a comment
4. 'BONDS' contains a list of the bond types used with format '(atom 1 number) (atom 2 number) [parameters]'
5. 'ANGLES' contains a list of the angle types used with format '(atom 1 number) (atom 2 number) (atom 3 number) [parameters]'
6. 'DIHEDRALS' contains a list of the dihedral types used with format '(atom 1 number) (atom 2 number) (atom 3 number) (atom 4 number) [parameters]'
7. 'VDW' contains a list of the van der Waal parameters with format '(atom 1 number) (atom 2 number) [parameters]'
8. 'ELECTRO' contains a list of the electrostatic parameters with format '(atom 1 number) (atom 2 number) [parameters]'
9. 'CONTROL' contains a list of other parameters in the LAMMPS input files that can be changed
    - 'units' has default 'real'
    - 'atom_style' has default 'pqeq'
    - 'boundary' has default 'f f f'
    - 'special_bonds' has default 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'


Example ffield

```shell
OPTIMIZER
method      1             # 0 None, 1 Optimize, 2 map
polar       1             # 0 off / 1 on
algorithm   Nelder-Mead   # options: SLSQP, Nelder-Mead
finish      SLSQP         # For map, options: SLSQP, Nelder-Mead, None
valence     0             # For creation of data files, 1 = valence copmlex, 0 = valence simple
FFIELD
uffqm
ATOMS
1 1.0080 # H_
2 12.0110 # C_3
BONDS
2 2 350 1.53 # C_3 C_3
2 1 350 1.09 # C_3 H_
ANGLES
1 2 1 40.0416 109.471 # H_ C_3 H_
1 2 2 49.7123 109.471 # H_ C_3 C_3
2 2 2 57.858 109.471 # X C_3 X
DIHEDRALS
1 2 2 1 0.145389 3 0 0 # H_ C_3 C_3 H_
1 2 2 X 0.167839 3 0 0 # H_ C_3 C_3 X
X 2 2 X 0.458156 3 0 0 # X C_3 C_3 X
VDW
1 1  0.06336  1.87834  3.27603 1.00000 1.89455 0.16379 0.00837 0.03493 0.02619
2 2  0.12734  2.42197  3.82580 1.00000 2.39594 1.12929 0.73792 0.27020 0.10393
1 2  0.12734  2.22197  3.52580 1.00000 2.15594 0.52929 0.35792 0.17020 0.06393
ELECTRO
1 1 4.52800 17.98410  0.37100  0  1.00000  0.37100 2037.20061  0.00000  
2 2 5.34300 10.12600  0.75900  0  1.00000  0.75900  198.84054  0.00000
CONTROL
special_bonds lj 0.0 1.0 1.0 coul 1.0 1.0 1.0
```

### param Formatting

Sections include 'BOUNDS' (contain the parameters to be optimized and any boundary constraints to impose) and 'CONSTRAINTS' (contain relations between parameters)

1. 'BOUNDS' has format '(parameter number) (ffield section) (ffield section line number) (ffield section column number) (lower bound) (upper bound) [step size if performing map]'
2. 'CONSTRAINTS' will take two parameter and force equivalence during the optimization and has format '(parameter number) (BOUNDS parameter number) (BOUNDS parameter number)'

Example param

```shell
BOUNDS
1 DIHEDRALS  1 5  0.0 1.0
2 DIHEDRALS  2 5  0.0 1.0
3 DIHEDRALS  3 5  0.0 1.0
4 VDW        1 3  0.0 1.0
5 VDW        1 5  2.5 3.5
6 ELECTRO    2 5  0.5 2.5
7 ELECTRO    2 8  0.5 2.5 
CONSTRAINTS
1 6 7
```


### trainset Formatting

Sections include 'ENERGY' (relative of absolute energies of structures) and 'CHARGE' (charges of atoms in structures)

1. 'ENERGY' has format '(weight) [(+ or -) (structure name)/(structure weight) ...] (energy)' and should be terminated with 'ENDENERGY'
2. 'CHARGE' has format '(structure) (weight) (atom ID number) (charge)' and should be terminated with 'ENDCHARGE'

Example trainset

```shell
ENERGY
#weight +/ structure/weight (... +/- structure/weight) energy
10.0    +  CH/1    -  CH2     200.5
1.0     +  CH2                  5.4
ENDENERGY
CHARGE
# molec_name weight atom_ID charge  
CH          1.0     1      -0.25
CH          1.0     7       0.15
CH          1.0     8       0.10
ENDCHARGE
```

### data files Formatting and Generation

- Data files should be named 'data_lammps_(structure name)'. 
- These files can be generated using the tool 'precompile_uffqm.py' (for uffqm forcefield) or including the '-g' option when calling the optimizer.
- The '-g' command will take any geometry file (concatenated .bgf files) and generate all of the proper data files depending on the force field name in 'ffield'.
- The 'precompile_uffqm.py' file will also generated 'temp.ffield' which can be useful in writing the 'ffield' file. If there are many atom/bond/angle/dihedral types in the data files, the 'temp.ffield' file will order all of the types properly. 

Support
-------

- Contact Saber Naserifar (naseri@caltech.edu)

License
-------

- TODO later

