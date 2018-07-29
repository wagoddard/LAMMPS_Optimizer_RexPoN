# Write the data.in file for the LAMMPS optimizer
# Assuming the force field is EXRYDC

def write_in(ffield, data_file):
    ''' Write the in file '''

    dict = ffield.get_ffield()

    ctrl = dict['CONTROL']

    # Retrieve commands
    units         = ctrl['units']
    atom_style    = ctrl['atom_style']
    boundary      = ctrl['boundary']
    special_bonds = ctrl['special_bonds']

    g = open('data.in', 'w')

    g.write('# EXRYDC - generated by optimizer\n')
    g.write('units            ' + units + '\n')
    g.write('atom_style       ' + atom_style + '\n')
    g.write('boundary         ' + boundary + '\n') 
    g.write('#dielectric      1\n')                                   # Commented out
    g.write('box tilt large\n')
    g.write('special_bonds    lj 0.0 0.0 1.0 coul 1.0 1.0 1.0\n')     # Change later
    g.write('\n')
    g.write('read_data        ' + data_file + '\n')
    g.write('\n')
    g.write('pair_style       hybrid/overlay coul/pqeqgauss 0.00 12.0 exrydc 0.0 12.0\n')  # Change later
    #g.write('pair_style       hybrid/overlay exrydc 0.0 12.5\n')  # Change later
    #g.write('pair_style       hybrid/overlay coul/pqeqgauss 0.00 12.50 \n')
    g.write('\n')
    g.write('#bond_style       harmonic\n')
    g.write('#angle_style      harmonic\n')
    g.write('#dihedral_style   dreiding\n')
    g.write('#improper_style   none\n')
    g.write('#kspace_style     none\n')
    g.write('\n')
    g.write('#           i j                  Xi       Ji       Rc       p  Qc       Rs          K2      K4\n')
    g.write('pair_coeff  * *  coul/pqeqgauss  0.00000  0.00000  0.00000  0  0.00000  0.00000     0.0000  0.00000  # dummy\n')
    for i in range(len(dict['ELECTRO'])):
        temp = dict['ELECTRO'][i]
        g.write('pair_coeff ' + str(temp[0]) + ' ' + str(temp[1]) + ' coul/pqeqgauss ' + ' '.join(map(str, temp[2:])) + '\n')
    g.write('#                       D R0 L alpha a0 a1 a2 a3 a4 a5\n')
    g.write('pair_coeff  * *  exrydc  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n')
    for i in range(len(dict['VDW'])):
        temp = dict['VDW'][i]
        g.write('pair_coeff ' + str(temp[0]) + ' ' + str(temp[1]) + ' exrydc ' + ' '.join(map(str, temp[2:])) + '\n')
    g.write('\n')
    g.write('pair_modify pair exrydc mix geometric\n')
    g.write('neigh_modify one 100000 page 1000000\n')
    g.write('# CHECK MIXING RULE FOR EXRYDC\n')
    g.write('\n')
    """
    for i in range(len(dict['BONDS'])):
        temp = dict['BONDS'][i]
        g.write('bond_coeff ' + str(i+1) + ' ' + ' '.join(map(str, temp[2:])) + '\n')
    for i in range(len(dict['ANGLES'])):
        temp = dict['ANGLES'][i]
        g.write('angle_coeff ' + str(i+1) + ' ' + ' '.join(map(str, temp[3:])) + '\n')
    for i in range(len(dict['DIHEDRALS'])):
        temp = dict['DIHEDRALS'][i]
        g.write('dihedral_coeff ' + str(i+1) + ' ' + ' '.join(map(str, temp[4:])) + '\n')
    """ 
    # PUT HYDROGEN BOND STYLE LATER
    g.write('\n')
    #g.write('fix              pqeq  all pqeq  1 0.0 20.0 1.0e-6\n')
    g.write('fix  pqeq all pqeq method 2 nevery 1 charge 0.0 tolerance 1.0e-6 damp 1.0\n')
    g.write('\n')
    g.write('compute          pqeq all pair coul/pqeqgauss\n')
    g.write('variable epqeq equal c_pqeq\n')

    g.write('thermo_style     custom etotal pe evdwl v_epqeq \n')

    g.write('compute 2 all property/atom q\n')
    g.write('compute 3 all property/atom type\n')
    g.write('compute 4 all property/atom id\n')
    g.write('variable charge equal q\n')
    g.write('variable typei equal type\n')
    g.write('variable idi equal id\n')
    g.write('compute fox all property/atom fx\n')
    g.write('compute foy all property/atom fy\n')
    g.write('compute foz all property/atom fz\n')

    g.write('variable etota equal etotal\n')
    g.write('compute eng all pe\n')

    g.write('run 0\n')
    g.write('\n')

    g.close()





