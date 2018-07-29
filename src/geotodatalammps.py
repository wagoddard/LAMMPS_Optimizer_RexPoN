# Parse a geo file and use tpascal script to convert to data_lammps_ files

import os
import sys
import glob
from pyvasek import bgf, periodic
import atomass
import logging
import string


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


def data_file_simple(geo):
    ''' rexpon '''

    logging.info('Converting geo to data_files, using rexpon...')

    gap = 50.0 # add to control file
    if os.path.isfile(geo):
       print '(1) Reading Geo...'
    else:
       print 'No File Exists'
       sys.quit()

    f = open(geo, 'r')
    lines = f.readlines()
    f.close()

    # Parse the geo file by BIOGRF
    out_array = []
    for line in lines:
        if line.startswith("BIOGRF"):
            out_array.append(line)
        else:
            out_array[-1] += line

    print '(2) Writing data files...'
    for i in range(len(out_array)):
        s = out_array[i]
        name = s[s.find('DESCRP')+len('DESCRP'):].split()[0]
        out_filename = 'data_lammps_%s'%name
        # Write a temporary .bgf file
        h = open('temp', 'w')
        h.write(s)
        h.close()

        lines = open('temp','r').readlines()
        al = bgf.read('temp')
        boxSize = periodic.readBgf(lines)

        xcoords = []; ycoords = []; zcoords = []
        natoms = 0
        types = []
        for a in al:
          xcoords.append(a.xyz[0])
          ycoords.append(a.xyz[1])
          zcoords.append(a.xyz[2])
          type = a.ffType.split('_')[0].strip()[:2]
          if type not in types:
            types.append(type)
          natoms += 1  
        types.sort()

        if boxSize != None:
            origin = boxSize[3,:] - (boxSize[0,:] + boxSize[1,:] + boxSize[2,:]) / 2
            xlo = origin[0]
            ylo = origin[1]
            zlo = origin[2]
            lx = boxSize[0,0]
            ly = boxSize[1,1]
            lz = boxSize[2,2]
            xy = boxSize[1,0]
            xz = boxSize[2,0]
            yz = boxSize[2,1]
            xhi = xlo + lx
            yhi = ylo + ly
            zhi = zlo + lz
            cell = [xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz]

        else:
            # min, max of coords for cell bboundaries
            minx = min(xcoords); maxx = max(xcoords)
            miny = min(ycoords); maxy = max(ycoords)
            minz = min(zcoords); maxz = max(zcoords)
            xlo = minx - gap
            ylo = miny - gap
            zlo = minz - gap
            xhi = maxx + gap
            yhi = maxy + gap
            zhi = maxz + gap
            yx = xz = yz = 0.0
            cell = [xlo, xhi, ylo, yhi, zlo, zhi, yx, xz, yz]

        # write the data file
        g = open(out_filename, 'w')
        g.write('lammps data file \n\n')
        g.write('%i atoms \n' % natoms)
        g.write('%i atom types \n' % len(types))
        g.write('\n')
        g.write('%f %f xlo xhi\n' % (cell[0], cell[1]))
        g.write('%f %f ylo yhi\n' % (cell[2], cell[3]))
        g.write('%f %f zlo zhi\n\n' % (cell[4], cell[5]))
        g.write('%f %f %f xy xz yz\n\n' % (cell[6], cell[7], cell[8]))
        g.write('Masses\n')
        g.write('\n')
        for i in range(len(types)):
          g.write('%i %s # %s\n' % (i+1,atomass.mass[types[i]],types[i]))
        g.write('\n')
        g.write('Atoms\n')
        g.write('\n')

        ID = 0
        mID = 1
        for a in al:
          ID += 1
          # molecID only for water system!!!!
          if ( ID%3 == 0.0):
            molecID = mID
            mID += 1
          else:
            molecID = mID
          # -----
          type = a.ffType.split('_')[0].strip()[:2]
          typeID = types.index(type)+1
          g.write('%i %i %i %f %f %f %f %f %f %f\n'%
                  (ID,molecID,typeID,0.0,a.xyz[0],a.xyz[1],a.xyz[2],0.0,0.0,0.0))
        g.close()

    logging.info('...conversion success.')

    return

def data_file_libattery(geo):
    print ' just for Li battery project '
    
    logging.info('Converting geo (.xyz)  to data_files...')

    gap = 0.0 # add to control file
    if os.path.isfile(geo):
       print '(1) Reading Geo...'
    else:
       print 'No File Exists'
       sys.quit()

    f = open(geo, 'r')
    lines = f.readlines()
    f.close()

    natoms = lines[0].split()[0]
    if  lines[0]:
      natoms = int(lines[0].split()[0])
    else:
      print 'geo file does not have correct format ...'
      sys.quit()


    # Parse the geo file by BIOGRF
    out_array = []
    for i in range(len(lines)):
        if i%(natoms+2) == 0:
            out_array.append(lines[i])
        else:
            out_array[-1] += lines[i]

    print '(2) converting each .xyz to .bgf using jaguar babel and then to data_lammps ...'
    #os.system('rm frame_*.bgf)
    for i in range(len(out_array)):
        s = out_array[i]
        h = open('temp', 'w')
        h.write(s)
        h.close()
        os.system('jaguar babel -ixyz temp -obgf frame_%i.bgf'%(i+1))
        os.system('rm temp')
        
        al = bgf.read('frame_%i.bgf'%(i+1))

        out_filename = 'data_lammps_%i'%(i+1)

        xcoords = []; ycoords = []; zcoords = []
        natoms = 0
        types = []

        # get info from each .bgf file
        for a in al:
          xcoords.append(a.xyz[0])
          ycoords.append(a.xyz[1])
          zcoords.append(a.xyz[2])
          #type = a.ffType.split('_')[0].strip()[:2]
          type = a.ffType.strip()[:]
          if type == 'H_HB': type = 'H_'

          if type not in types:
            types.append(type)
          natoms += 1
        types.sort()

        for j in range(len(types)):
          if types[j] == 'H_HB':
            types[j] = 'H_'
          if not (types[j] == 'C_3' or types[j] == 'H_' or \
             types[j] == 'Li' or types[j] == 'N_2' or \
             types[j] == 'N_3' or types[j] == 'O_2' or \
             types[j] == 'O_3') :
             print 'unknown atom type %s:  '%types[j]
             sys.exit()

        # box info should be known and same for all frames
        minx = 0.0;  maxx = 18.0703
        miny = 0.0;  maxy = 18.0703
        minz = 0.0;  maxz = 18.0703

        g = open(out_filename, 'w')
        g.write('lammps data file \n\n')
        g.write('%i atoms \n'%natoms)
        g.write('%i atom types \n'%len(types))
        g.write('\n')
        g.write('%f %f xlo xhi\n'%(minx-gap,maxx+gap))
        g.write('%f %f ylo yhi\n'%(miny-gap,maxy+gap))
        g.write('%f %f zlo zhi\n\n'%(minz-gap,maxz+gap))
        g.write('Masses\n')
        g.write('\n')
        for k in range(len(types)):
          g.write('%i %s # %s\n'%(k+1,atomass.mass[types[k]],types[k]))
        g.write('\n')
        g.write('Atoms\n')
        g.write('\n')
        ID = 0
        for a in al:
          ID += 1
          #type = a.ffType.split('_')[0].strip()[:2]
          type = a.ffType.strip()[:]
          if type == 'H_HB': type = 'H_'

          typeID = types.index(type)+1
          g.write('%i %i %i %f %f %f %f %f %f %f\n'%
                  (ID,1,typeID,0.0,a.xyz[0],a.xyz[1],a.xyz[2],0.0,0.0,0.0))
        g.close()
        os.system('rm frame_%i.bgf'%(i+1))
        if ( (i+1)%100!=0 and i!=0 ):
          os.system('rm data_lammps_%i'%(i+1))

    logging.info('...conversion success.')

    return


