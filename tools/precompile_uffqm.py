# Parse a geo file and use tpascal script to convert to data_lammps_ files

import os
import sys
import glob
import logging


def data_file_complete(geo):
    ''' Dreiding - Take a geo file of .bgf's and convert to data files '''

    logging.info('Converting geo to data_files, using uffqm...')

    if os.path.isfile(geo):
       print '(1) Reading Geo...'
    else:
       print 'No File Exists'
       sys.quit()

    f = open(geo, 'r')
    lines = f.readlines()
    f.close()

    # --------------------------------------------------

    # Parse the geo file by BIOGRF
    out_array = []
    for line in lines:
        if line.startswith("BIOGRF"):
            out_array.append(line)
        else:
            out_array[-1] += line

    # -------------------------------------------

    # Create temp .bgf files

    out_name = []
    for i in range(len(out_array)):
        s = out_array[i]
        name = s[s.find('DESCRP') + len('DESCRP'):].split()[0]
        out_filename = name
        out_name.append(out_filename)

        # TODO: Check that file does not already exist

        g = open(out_filename, 'w')
        g.write(out_array[i])
        g.close()

    # --------------------------------------------

    # Create data.lammps files and generate atom dictionary

    # TODO: use a function rather than repeating the same code

    atom_lst  = []
    bond_lst  = []
    angl_lst  = []
    dihe_lst  = []
    impr_lst  = []

    lst       = ['Created', 'Masses', 'Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers']
    lst2      = ['Bond', 'Angle', 'Pair', 'Dihedral', 'Improper']

    atom_dic  = {}
    bond_dic  = {}
    angl_dic  = {}
    dihe_dic  = {}
    impr_dic  = {}

    for out_filename in out_name:
        atom_dic[out_filename] = {}
        bond_dic[out_filename] = {}
        angl_dic[out_filename] = {}
        dihe_dic[out_filename] = {}
        impr_dic[out_filename] = {}

        # Run Tod's script to generate data.lammps files
        os.system('/ul/joppenhe/lammps_optimizer/lib/createLammpsInput.pl -b %s -f \"DREIDING\"' % out_filename)
        os.remove(out_filename)

        g = open('data.lammps', 'r')
        lines2 = g.readlines()
        g.close()

        # Run through each file and create dictionaries and lists to renumber everything
        for j in range(len(lines2)):
            temp = lines2[j].split()

            if temp == ['Masses']:
                k = 1
                temp2 = lines2[j+k].split()
                while j+k<len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                    if len(temp2) != 0 and ' '.join(temp2[1:]) not in atom_lst:
                        atom_lst.append(' '.join(temp2[1:]))
                    if len(temp2) != 0 and temp2[0] not in atom_dic[out_filename]:
                        atom_dic[out_filename][temp2[0]] = atom_lst.index(' '.join(temp2[1:])) + 1
                    k += 1
                    temp2 = lines2[j+k].split()
            elif temp == ['Bond', 'Coeffs']:
                k = 1
                temp2 = lines2[j+k].split()
                while j+k<len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                    if len(temp2) != 0 and ' '.join(temp2[1:]) not in bond_lst:
                        bond_lst.append(' '.join(temp2[1:]))
                    if len(temp2) != 0 and temp2[0] not in bond_dic[out_filename]:
                        bond_dic[out_filename][temp2[0]] = bond_lst.index(' '.join(temp2[1:])) + 1
                    k += 1
                    temp2 = lines2[j+k].split()
            elif temp == ['Angle', 'Coeffs']:
                k = 1
                temp2 = lines2[j+k].split()
                while j+k<len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                    if len(temp2) != 0 and ' '.join(temp2[1:]) not in angl_lst:
                        angl_lst.append(' '.join(temp2[1:]))
                    if len(temp2) != 0 and temp2[0] not in angl_dic[out_filename]:
                        angl_dic[out_filename][temp2[0]] = angl_lst.index(' '.join(temp2[1:])) + 1
                    k += 1
                    temp2 = lines2[j+k].split()
            elif temp == ['Dihedral', 'Coeffs']:
                k = 1
                temp2 = lines2[j+k].split()
                while j+k<len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                    if len(temp2) != 0 and ' '.join(temp2[1:]) not in dihe_lst:
                        dihe_lst.append(' '.join(temp2[1:]))
                    if len(temp2) != 0 and temp2[0] not in dihe_dic[out_filename]:
                        dihe_dic[out_filename][temp2[0]] = dihe_lst.index(' '.join(temp2[1:])) + 1
                    k += 1
                    temp2 = lines2[j+k].split()
            elif temp == ['Improper', 'Coeffs']:
                k = 1
                temp2 = lines2[j+k].split()
                while j+k<len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                    if len(temp2) != 0 and ' '.join(temp2[1:]) not in impr_lst:
                        impr_lst.append(' '.join(temp2[1:]))
                    if len(temp2) != 0 and temp2[0] not in impr_dic[out_filename]:
                        impr_dic[out_filename][temp2[0]] = impr_lst.index(' '.join(temp2[1:])) + 1
                    k += 1
                    temp2 = lines2[j+k].split()

        # I use 'cp' rather than 'mv' for debugging purposes - it is deleted later on
        os.system('cp data.lammps data_lammps_%s' % out_filename)

    # For debugging purposes, to check if all dictionaries and lists have been made
    logging.debug('atom_lst: %s' % str(atom_lst))
    logging.debug('atom_dic: %s' % str(atom_dic))
    logging.debug('bond_lst: %s' % str(bond_lst))
    logging.debug('bond_dic: %s' % str(bond_dic))
    logging.debug('angl_lst: %s' % str(angl_lst))
    logging.debug('angl_dic: %s' % str(angl_dic))
    logging.debug('dihe_lst: %s' % str(dihe_lst))
    logging.debug('dihe_dic: %s' % str(dihe_dic))
    logging.debug('impr_lst: %s' % str(impr_lst))
    logging.debug('impr_dic: %s' % str(impr_dic))

    # --------------------------------------------------

    # Rewrite data.lammps file without lst2 sections and with new numbers

    dirs = glob.glob("data_lammps_*")

    # Loop over every data_lammps_ file
    for i in range(len(dirs)):

        out_filename = '_'.join(dirs[i].split('_')[2:])

        g = open('data_lammps_%s' % out_filename, 'r')
        lines2 = g.readlines()
        g.close()

        g = open('data_lammps_%s' % out_filename, 'w')

        # Loop over every line and rewrite
        for j in range(len(lines2)):
            temp = lines2[j].split()
            if len(temp) > 0:
                if temp[0] == 'Created':
                    k = 1
                    g.write(lines2[j])
                    temp2 = lines2[j+k].split()
                    while j+k < len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                        if len(lines2[j+k].split()) == 3:
                            if lines2[j+k].split()[1] == 'atom':
                                g.write(' '*11 + str(len(atom_lst)) + ' atom types\n')
                            elif lines2[j+k].split()[1] == 'bond':
                                g.write(' '*11 + str(len(bond_lst)) + ' bond types\n')
                            elif lines2[j+k].split()[1] == 'angle':
                                g.write(' '*11 + str(len(angl_lst)) + ' angle types\n')
                            elif lines2[j+k].split()[1] == 'dihedral':
                                g.write(' '*11 + str(len(dihe_lst)) + ' dihedral types\n')
                            elif lines2[j+k].split()[1] == 'improper':
                                g.write(' '*11 + str(len(impr_lst)) + ' improper types\n')
                            else: 
                                g.write(lines2[j+k])
                        else:
                            g.write(lines2[j+k])
                        k += 1
                        if j+k < len(lines2):
                            temp2 = lines2[j+k].split()
                elif temp[0] == 'Masses':
                    g.write('Masses\n\n')
                    for k in range(len(atom_lst)):
                        g.write('    ' + str(k+1) + ' ' + atom_lst[k] + '\n')
                    g.write('\n')
                elif temp[0] == 'Atoms':
                    k = 1
                    g.write('Atoms\n')
                    temp2 = lines2[j+k].split()
                    while j+k < len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                        if len(temp2) == 0:
                            g.write('\n')
                        else:
                            g.write('    ' + temp2[0] + ' ' + temp2[1] + ' ' + str(atom_dic[out_filename][str(temp2[2])]) + ' ' + ' '.join(temp2[3:]) + '\n')
                        k += 1
                        if j + k < len(lines2):
                            temp2 = lines2[j+k].split()
                elif len(temp) == 1 and temp[0] == 'Bonds':
                    k = 1
                    g.write('Bonds\n')
                    temp2 = lines2[j+k].split()
                    while j+k < len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                        if len(temp2) == 0:
                            g.write('\n')
                        else:
                            g.write('    ' + temp2[0] + ' ' + str(bond_dic[out_filename][str(temp2[1])]) + ' ' + ' '.join(temp2[2:]) + '\n')
                        k += 1
                        if j + k < len(lines2):
                            temp2 = lines2[j+k].split()
                elif len(temp) == 1 and temp[0] == 'Angles':
                    k = 1
                    g.write('Angles\n')
                    temp2 = lines2[j+k].split()
                    while j+k < len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                        if len(temp2) == 0:
                            g.write('\n')
                        else:
                            g.write('    ' + temp2[0] + ' ' + str(angl_dic[out_filename][str(temp2[1])]) + ' ' + ' '.join(temp2[2:]) + '\n')
                        k += 1
                        if j + k < len(lines2):
                            temp2 = lines2[j+k].split()
                elif len(temp) == 1 and temp[0] == 'Dihedrals':
                    k = 1
                    g.write('Dihedrals\n')
                    temp2 = lines2[j+k].split()
                    while j+k < len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                        if len(temp2) == 0:
                            g.write('\n')
                        else:
                            g.write('    ' + temp2[0] + ' ' + str(dihe_dic[out_filename][str(temp2[1])]) + ' ' + ' '.join(temp2[2:]) + '\n')
                        k += 1
                        if j + k < len(lines2):
                            temp2 = lines2[j+k].split()
                elif len(temp) == 1 and temp[0] == 'Impropers':
                    k = 1
                    g.write('Impropers\n')
                    temp2 = lines2[j+k].split()
                    while j+k < len(lines2) and (len(temp2) == 0 or temp2[0] not in (lst + lst2)):
                        if len(temp2) == 0:
                            g.write('\n')
                        else:
                            g.write('    ' + temp2[0] + ' ' + str(impr_dic[out_filename][str(temp2[1])]) + ' ' + ' '.join(temp2[2:]) + '\n')
                        k += 1
                        if j + k < len(lines2):
                            temp2 = lines2[j+k].split()

        g.close()

    # ------------------------------------

    # Write dictionaries to a file
    # This is important for getting FFIELD correct

    g = open('temp.ffield', 'w')
    # Atoms
    g.write('ATOMS\n')
    for a in atom_lst:
        g.write(a + '\n')
    # Bonds        
    g.write('BONDS\n')
    for b in bond_lst:
        g.write(b + '\n')
    # Angles
    g.write('ANGLES\n')
    for a in angl_lst:
        g.write(a + '\n')
    # Dihedrals
    g.write('DIHEDRALS\n')
    for d in dihe_lst:
        g.write(d + '\n')
    # Impropers
    g.write('IMPROPERS\n')
    for i in impr_lst:
        g.write(i + '\n')
    g.close()


    # ----------------------------------------------------------


    # Clean up files
    os.remove('data.lammps')
    #os.remove('in.lammps')
    #os.remove('in.lammps_singlepoint')
    #os.remove('lammps.lammps.pbs')

    logging.info('...conversion success.')

    return

if __name__=="__main__":
    if len(sys.argv) != 2:
        print "Usage: python precompile_uffqm.py [geo]"
        print " geo: a concatenation of .bgf files to compile into LAMMPS data files with proper numbering and format for UFFQM"
        print " outputs: data_lammps_* for each structure, temp.ffield for building ffield file with proper numbering"
        sys.exit()

    geo = sys.argv[1]

    data_file_complete(geo)




