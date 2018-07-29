######################################################
# Optimizer for LAMMPS                              ##
# Based on optimization functions in scipy          ##
# Inputs: data.lammps, trainset, ffield, params     ##
# Outputs: Optimal parameters (history)             ##
######################################################

############## TO-DO ###############################################
# All To-do sections marked by many #####                         ##
# 1 - Implement in parallel rather than serial                    ##
####################################################################


# LAMMPS
from lammps import lammps                                        # Python-LAMMPS interface

# SYSTEM PACKAGES
import sys                                                       # Command-line inputs, exit()
import glob                                                      # Find files in the working directory
import logging                                                   # Keep a log to debug the optimizer
import time                                                      # Used to get the time of optimization
from ctypes import *                                             # LAMMPS extract functions use ctypes
import re                                                        # Needed to modify some strings

# MATH
from scipy.optimize import minimize, brute                       # Optimization functions, Assuming current scipy version
import numpy                                                     # Used for some arrays

# SRC CODE
import forcefield                                                # Used to store force field data
import trainset                                                  # Used to store trainset data
import param                                                     # Used to store parameter data
import geotodatalammps                                           # Used to convert geo file to LAMMPS data files

# DATA FILE GENERATORS
import write_uffqm                                               # Used to write the data.in files
import write_rexpon                                               # ...
import write_morse                                               # ...
import write_ljcut                                               # ...
import write_exryd, write_exrydc                                 # ...
import write_pprlg                                               # ...
import write_opls                                                # ...



def prepare_eng_table(ff, polar):
    ''' Prepare energy_table.txt '''

    if not polar:
        if ff.get_name().lower() == 'uffqm':
            logging.info('Writing energy_table.txt output...')
            h = open('energy_table.txt', 'w')
            h.write("%-30s  ETOTAL   PE        KE       EVDWL   ELECTRO   EBOND    EANGLE   EDIHED\n" % ('#str'))
        elif ff.get_name().lower() == 'rexpon':
            logging.info('Writing energy_table.txt output...')
            h = open('energy_table.txt', 'w')
            h.write('%-30s  Ep       Eb        Evdw     Ehb     Elp       Eangle   Etors    Eelectro\n' % ('#str'))
        elif ff.get_name().lower() == 'opls':
            logging.info('Writing energy_table.txt output...')
            h = open('energy_table.txt', 'w')
            h.write("%-30s  ETOTAL   PE        KE       EVDWL   ELECTRO   EBOND    EANGLE   EDIHED\n" % ('#str'))
        elif ff.get_name().lower() == 'exrydc':
            logging.info('Writing energy_table.txt output...')
            h = open('energy_table.txt', 'w')
            h.write("%-30s  ETOTAL  PE  ELECTRO\n" % ('#str'))
        elif ff.get_name().lower() == 'pprlg':
            logging.info('Writing energy_table.txt output...')
            h = open('energy_table.txt', 'w')
            h.write("%-30s  ETOTAL  PE  ELECTRO\n" % ('#str'))
        else:
            logging.warning('Neither uffqm, rexpon, or opls supplied for energy_table.txt')
            h = open('energy_table.txt', 'w')
            h.write('Warning, neither uffqm, rexpon, or opls supplied\n')
    else:
        logging.warning('Polar Calculation - Cannot write energy_table.txt')
        h = open('energy_table.txt', 'w')
        h.write('Warning, cannot write energy_table.txt for polar calculation\n')

    return h



def write_in_file(ff, data, polar):
    ''' Write the proper data.in file depending on the forcefield type '''

    ff_name = ff.get_name().lower()
    if ff_name == 'uffqm':                                   # Write the data.in file
        write_uffqm.write_in(ff, data, polar)             #  ... uffqm
    elif ff_name == 'rexpon':
        write_rexpon.write_in(ff, data)                    #  ... rexpon
    elif ff_name == 'morse':
        write_morse.write_in(ff, data)                    #  ... morse
    elif ff_name == 'ljcut':
        write_ljcut.write_in(ff, data)                    #  ... lj/cut
    elif ff_name == 'opls':                                 
        write_opls.write_in(ff, data)                     #  ... OPLS
    elif ff_name == 'exryd':
        write_exryd.write_in(ff, data)
    elif ff_name == 'exrydc':
        write_exrydc.write_in(ff, data)
    elif ff_name == 'pprlg':
        write_pprlg.write_in(ff, data)
    else:
        print "Forcefield name does not have a data.in writing function"
        logging.error("Forcefield name does not have a data.in writing function")
        sys.exit()

    return


def write_history(param, err):
    ''' Write the history file. '''

    to_write = ""
    for item in param:
        to_write += '%0.6f\t' % item
    to_write += str(err) + '\n'

    g = open("history", 'a')                                     # Create a history file of all errors
    g.write(to_write)
    g.close()

    return



def compute_error(train, energ_dict, charg_dict, force_dict, ID_dict):
    ''' Compute the error between energy dictionary and trainset. '''

    error_lst = {}                                               # Store the error for each training set line
    error_lst['energy'] = []
    error_lst['charge'] = []
    error_lst['force'] = []
    lammp_lst = []                                               #   ... for the LAMMPS energy
    train_lst = []                                               #   ... for the QM energy

    # Calculate error for energies
    for i in range(len(train.get_energy())):
        weight = float(train.get_energy()[i][0])
        struct = train.get_energy()[i][2::2]
        sign   = train.get_energy()[i][1:-1:2]
        tr_err = float(train.get_energy()[i][-1])

        lammps_en = 0.0

        for j in range(len(struct)):
            sign[j] = -1 if sign[j] == '-' else 1

            # Calculate the LAMMPS energy factoring weighting and sign 
            try:
                lammps_en += energ_dict[struct[j].split('/')[0]] * sign[j] / float(struct[j].split('/')[-1])
            except Exception as e:
                print "Error matching trainset to structures, ", e
                logging.error("Error matching trainset to structures, %s" % e)
                sys.exit()

        #err = ((lammps_en - tr_err) / (weight * tr_err)) ** 2             # Calculate the error
        err = ((lammps_en - tr_err) / (weight)) ** 2             # Calculate the error

        lammp_lst.append(lammps_en)                              # Save for fort.99
        train_lst.append(tr_err)                                 # Save for fort.99
        error_lst['energy'].append(err)                          # Save for fort.99 and output of function

    # Calculate error for charges
    for i in range(len(train.get_charge())):
        struct = train.get_charge()[i][0]
        weight = float(train.get_charge()[i][1])
        at_no  = int(train.get_charge()[i][2])
        tr_chg = float(train.get_charge()[i][3])

        try:
            lammps_chg = charg_dict[struct][at_no-1]
        except Exception as e:
            print "Error matching trainset charge to structures, ", e
            logging.error("Error matching trainset charge to structures, %s" %e)
            sys.exit()

        err = ((lammps_chg - tr_chg) / (weight)) ** 2

        lammp_lst.append(lammps_chg)
        train_lst.append(tr_chg)
        error_lst['charge'].append(err)

    # Calculate error for forces
    for i in range(len(train.get_force())):
        struct = train.get_force()[i][0]
        weight = float(train.get_force()[i][1])
        at_no  = int(train.get_force()[i][2]) 
        at_indx = ID_dict[struct].index(at_no)
        direct = train.get_force()[i][3]
        tr_for = float(train.get_force()[i][4])

        dir = -1
        if direct.lower() == 'x':
            dir = 0
        elif direct.lower() == 'y':
            dir = 1
        elif direct.lower() == 'z':
            dir = 2
        elif direct.lower() == 't':
            dir = 3
        else:
            print 'ERROR: Improper direction given for FORCES (x,y,z,t)'
            logging.error('Improper direction given for FORCES (x,y,z,t)')
            sys.exit()

        try:
            #lammps_chg = force_dict[struct][dir][at_no-1]
            lammps_chg = force_dict[struct][dir][at_indx]
        except Exception as e:
            print "Error matching trainset force to structures, ", e
            logging.error("Error matching trainset force to structures, %s" % e)
            sys.exit()

        err = ((lammps_chg - tr_for) / (weight)) ** 2

        lammp_lst.append(lammps_chg)
        train_lst.append(tr_for)
        error_lst['force'].append(err)

    return lammp_lst, error_lst, train_lst



def write_fort_99(train, lammp_lst, train_lst, error_lst):
    ''' Write the fort.99 file. '''

    logging.info('Writing fort.99 output...')

    g = open("fort.99", 'w')
    g.write('#str' + ' ' * (46) + 'LAMMPS   QM     WEIGHT  ERROR   TOTAL ERROR\n')

    a = len(train.get_energy())
    b = len(train.get_charge())
    c = len(train.get_force())

    for i in range(a):
        structure = ' '.join(train.get_energy()[i][0:-1])
        write_lmp = lammp_lst[i]
        write_trn = train_lst[i]
        weight    = train.get_energy()[i][0]
        write_err = error_lst['energy'][i]
        write_ter = sum(error_lst['energy'][0:i+1])

        g.write('%s %-40s %8.4f %8.4f %-5s %8.4f %8.4f\n' % ('ENERGY', structure, write_lmp, write_trn, weight, write_err, write_ter))

    tot_en_err = sum(error_lst['energy'])

    for i in range(b):
        structure = ' '.join(train.get_charge()[i][0:-1])
        write_lmp = lammp_lst[i+a]
        write_trn = train_lst[i+a]
        weight    = train.get_charge()[i][1]
        write_err = error_lst['charge'][i]
        write_ter = sum(error_lst['charge'][0:i+1]) + tot_en_err

        g.write('%s %-40s %8.4f %8.4f %-5s %8.4f %8.4f\n' % ('CHARGE', structure, write_lmp, write_trn, weight, write_err, write_ter))

    tot_en_err += sum(error_lst['charge'])

    for i in range(c):
        structure = ' '.join(train.get_force()[i][0:-1])
        write_lmp = lammp_lst[i+a+b]
        write_trn = train_lst[i+a+b]
        weight    = train.get_force()[i][1]
        write_err = error_lst['force'][i]
        write_ter = sum(error_lst['force'][0:i+1]) + tot_en_err
   
        g.write('%s %-40s %8.4f %8.4f %-5s %8.4f %8.4f\n' % ('FORCE ', structure, write_lmp, write_trn, weight, write_err, write_ter))

    g.close()

    logging.info('...writing success.')

    return



def write_eng_table(lmp, h, name, data):
    ''' Write the energy table given a ffield type. '''

    if name == 'uffqm':
        etotal  = lmp.extract_variable("etota", 0, 0)         # Total energy
        pe      = lmp.extract_compute("eng", 0, 0)               # Potential energy
        ke      = lmp.extract_variable("ke2", 0, 0)              # Kinetic energy
        evdwl   = lmp.extract_variable("evdw", 0, 0)             # Van der Waal energy
        electro = lmp.extract_variable("epqeq", 0, 0)            # Electrostatic energy
        ebond   = lmp.extract_variable("ebon", 0, 0)             # Bond energy
        eangle  = lmp.extract_variable("eangl", 0, 0)            # Angle energy
        edihed  = lmp.extract_variable("edihe", 0, 0)            # Dihedral energy

        #phi     = lmp.extract_variable("phi2", 0, 1)            # List of dihedral angles
        #phi     = lmp.extract_compute("1", 2, 1)[1]             #   ... extract a single dihedral

        h.write('%-30s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (data.split("data_lammps_")[1], etotal, pe, ke, evdwl, electro, ebond, eangle, edihed))

    elif name == 'rexpon':
        pe      = lmp.extract_compute("eng", 0, 0)               # Potential energy
        evdwl   = lmp.extract_variable("ew", 0, 0)               # Van der Waal energy
        electro = lmp.extract_variable("epqeq", 0, 0)            # Electrostatic energy
        ebond   = lmp.extract_variable("eb", 0, 0)               # Bond energy
        eangle  = lmp.extract_variable("ev", 0, 0)               # Angle energy
        etors   = lmp.extract_variable("et", 0, 0)               # Torsion energy
        ehydb   = lmp.extract_variable("ehb", 0, 0)              # Hydrogen bond energy
        elpb    = lmp.extract_variable("elp", 0, 0)              # ?????????????

        h.write('%-30s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (data.split("data_lammps_")[1], pe, ebond, evdwl, ehydb, elpb, eangle, etors, electro))         

    elif name == 'opls':
        etotal  = lmp.extract_variable("etota", 0, 0)            # Total energy
        pe      = lmp.extract_compute("eng", 0, 0)               # Potential energy
        ke      = lmp.extract_variable("ke2", 0, 0)              # Kinetic energy
        evdwl   = lmp.extract_variable("evdw", 0, 0)             # Van der Waal energy
        electro = lmp.extract_variable("epqeq", 0, 0)            # Electrostatic energy
        ebond   = lmp.extract_variable("ebon", 0, 0)             # Bond energy
        eangle  = lmp.extract_variable("eangl", 0, 0)            # Angle energy
        edihed  = lmp.extract_variable("edihe", 0, 0)            # Dihedral energy

        h.write('%-30s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (data.split("data_lammps_")[1], etotal, pe, ke, evdwl, electro, ebond, eangle, edihed))

    elif name == 'opls':
        etotal  = lmp.extract_variable("etota", 0, 0)            # Total energy
        pe      = lmp.extract_compute("eng", 0, 0)               # Potential energy
        ke      = lmp.extract_variable("ke2", 0, 0)              # Kinetic energy
        evdwl   = lmp.extract_variable("evdw", 0, 0)             # Van der Waal energy
        electro = lmp.extract_variable("epqeq", 0, 0)            # Electrostatic energy
        ebond   = lmp.extract_variable("ebon", 0, 0)             # Bond energy
        eangle  = lmp.extract_variable("eangl", 0, 0)            # Angle energy
        edihed  = lmp.extract_variable("edihe", 0, 0)            # Dihedral energy

        h.write('%-30s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (data.split("data_lammps_")[1], etotal, pe, ke, evdwl, electro, ebond, eangle, edihed))

    elif name == 'exrydc':
        etotal  = lmp.extract_variable("etota", 0, 0)            # Total energy
        pe      = lmp.extract_compute("eng", 0, 0)               # Potential energy
        electro = lmp.extract_variable("epqeq", 0, 0)            # Electrostatic energy

        h.write('%-30s %8.4f %8.4f %8.4f\n' % (data.split("data_lammps_")[1], etotal,pe,electro))

    elif name == 'pprlg':
        #etotal  = lmp.extract_variable("etota", 0, 0)            # Total energy
        pe      = lmp.extract_compute("eng", 0, 0)               # Potential energy
        #electro = lmp.extract_variable("epqeq", 0, 0)            # Electrostatic energy

        #h.write('%-30s %8.4f %8.4f %8.4f\n' % (data.split("data_lammps_")[1], etotal,pe,electro))
        h.write('%-30s %8.4f\n' % (data.split("data_lammps_")[1],pe))

    else:
        h.write('Warning: Neither uffqm, rexpon, or opls specified\n')
        logging.warning('Neither uffqm, rexpon, or opls specified for energy table.')


    return



def write_chg_table(lmp, m, ff, qm, data):
    ''' Write the charge table '''

    # Note the QM list starts at index 0, but trainset starts at 1

    for i in range(len(qm)):
        write_nm = data.split("data_lammps_")[1]
        if qm[i][0] == write_nm:
            write_qm = float(qm[i][3])
            try:
                write_ff = float(ff[qm[i][0]][int(qm[i][2])-1])
            except Exception as e:
                print 'Charge error matching trainset to data files: ' + str(qm[i][0])
                logging.error('Charge error matching trainset to data files: %s' % str(qm[i][0]))
                sys.exit()
            m.write('%-30s %8.4f %8.4f\n' % (write_nm, write_qm, write_ff))


    return



def write_force_table(lmp, n, ff, qm, data):
    ''' Write the force table '''


    for i in range(len(qm)):
        write_nm = data.split("data_lammps_")[1]

        if qm[i][0] == write_nm:
            write_qm = float(qm[i][4])

            direct = qm[i][3]
            dir = -1
            if direct.lower() == 'x':
                dir = 0
            elif direct.lower() == 'y':
                dir = 1
            elif direct.lower() == 'z':
                dir = 2
            elif direct.lower() == 't':
                dir = 3
            else:
                print 'ERROR: Improper direction given for FORCES (x,y,z,t)'
                logging.error('Improper direction given for FORCES (x,y,z,t)')
                sys.exit()

            try:
                write_ff = float(ff[qm[i][0]][dir][int(qm[i][2])-1])
            except Exception as e:
                print 'Force error matching trainset to data files: ' + str(qm[i][0])
                logging.error('Force error matching trainset to data files: %s' % str(qm[i][0]))
                sys.exit()
 
            n.write('%-30s %8.4f %8.4f\n' % (write_nm, write_qm, write_ff))

    return



def write_out(ff, param, param_loc):
    ''' Write results of the optimization to a force field to fort.4 '''

    logging.info("Writing fort.4 file...")

    ff_temp = ff.get_ffield()

    for k in range(len(param)):                                  # Get force field with final parameters
        ff_temp[param_loc[k][0]][param_loc[k][1] - 1][param_loc[k][2] - 1] = str(param[k])
        ff.set_ffield(ff_temp)

    g = open('fort.4', 'w')

    g.write('OPTIMIZER\n')
    for key, value in ff.get_optimizer().iteritems():
        g.write(key + ' ' + value + '\n')

    g.write('FFIELD\n')
    g.write(ff.get_name() + '\n')

    for x in ff.get_ffield():                                    # Write sections in no particular order
        g.write(x + '\n')
        if x != 'CONTROL':
            for y in range(len(ff.get_ffield()[x])):
                g.write(' '.join(ff.get_ffield()[x][y]) + '\n')
        if x == 'CONTROL':
            for key, value in ff.get_ffield()[x].iteritems():
                g.write(key + ' ' + value + '\n')

    g.close()    

    logging.info("...writing success.")

    return



def readData():
    ''' Generate a list of all data_lammps files in the working directory. '''

    logging.info("Reading data files...")

    data = glob.glob("data_lammps_*")                            # File stored as list

    if len(data) == 0:
        print "Error: no data files"
        logging.error("Error: no data files")
        sys.exit()

    logging.debug('readData(): %s' % str(data))
    logging.info("...reading success.")

    return data



def optimizer(param, data, train, ff, param_loc, polar, res_print=False, debugging=True, disp_screen=False, cite=False):
    ''' Optimizer function that runs lammps based on current values
        and outputs energy to be minimized. param are the current
        values of the parameters being optimized.
    '''

    #param[-1] = param[-2]

    if debugging:
        logging.debug("Optimizer data: %s\nEntering optimizer function..." % str(data))

    # ---------------------------------------------------------------------------------------------

    # Update the force field with current optimizer parameters
    ff_temp = ff.get_ffield()
    try:
        float(param)
        a = 0
    except (ValueError, TypeError) as e:
        a = 1
    if a == 0:
        param = [param]
    for k in range(len(param)):
        ff_temp[param_loc[k][0]][param_loc[k][1] - 1][param_loc[k][2] - 1] = str(param[k])
        ff.set_ffield(ff_temp)

    # ---------------------------------------------------------------------------------------------

    if res_print: 
        h = prepare_eng_table(ff, polar)                         # Prepare energy_table.txt

        logging.info('Writing charge_table.txt output...')       # Prepare charge_table.txt
        m = open('charge_table.txt', 'w')
        m.write('%-30s QM-Charge FF-Charge\n' % ('#data'))

        logging.info('Writing force_table.txt output...')        # Prepare force_table.txt
        n = open('force_table.txt', 'w')
        n.write('%-30s QM-force FF-Force\n' % ('#data')) 

    # ---------------------------------------------------------------------------------------------

    energ_dict = {}                                              # Store a dictionary of energies of structures
    charg_dict = {}                                              # Store a dictionary of charges for atoms in structures
    force_dict = {}                                              # Store a dictionary of forces for each component of an atom (x,y,z,tot)
    ID_dict = {}

    # For each structure
    for s in range(len(data)):
        cmdarg = []
        if not cite:                                             # Do not create log.cite
            cmdarg.append('-nocite')
        if not disp_screen:                                      # Do not display output to terminal screen
            cmdarg = cmdarg + ['-screen','none']

        write_in_file(ff, data[s], polar)                        # Write data.in file
        lmp = lammps(cmdargs=cmdarg)                             # Create a new instance of LAMMPS
        lmp.file('data.in')                                      # Call LAMMPS
 
        with open('log.lammps', 'r') as f:
            for line in f:
                if line == 'WARNING: Warning: core-shell dist more than 0.250 A (../fix_pqeq.cpp:737)\n':
                    write_history(param, None)
                    return 77777777 #None                        # Arbitrary large number

        # Sabers new code
        energy     = lmp.extract_compute("eng",0,0)              # Extract energy
        charge     = lmp.extract_compute("2", 1, 1)              # Extract charge
        natoms     = lmp.get_natoms()                            # Extract number of atoms
        types      = lmp.extract_compute("3", 1, 1)              # Extract types
        IDs        = lmp.extract_compute("4", 1, 1)              # Extract types
        fox        = lmp.extract_compute("fox", 1, 1)
        foy        = lmp.extract_compute("foy", 1, 1)
        foz        = lmp.extract_compute("foz", 1, 1)

        charge_lst = [charge[i] for i in xrange(natoms)]         # Convert C pointer to double array to python list
        type_lst   = [types[i] for i in xrange(natoms)]          # Convert C pointer to double array to python list
        ID_lst     = [int(IDs[i]) for i in xrange(natoms)]
        fox_lst    = [fox[i] for i in xrange(natoms)]
        foy_lst    = [foy[i] for i in xrange(natoms)]
        foz_lst    = [foz[i] for i in xrange(natoms)]
        fot_lst    = [(fox[i] ** 2 + foy[i] ** 2 + foz[i] ** 2) ** (0.5) for i in xrange(natoms)]

        force_lst = [fox_lst, foy_lst, foz_lst, fot_lst]

        if polar:
            energy_p = lmp.extract_compute("5", 0, 0)
            energy = energy - energy_p


        energ_dict[data[s].split("data_lammps_")[1]] = energy    # Store the energies into the dictionary
        charg_dict[data[s].split("data_lammps_")[1]] = charge_lst# Store the charges into the dictionary
        force_dict[data[s].split("data_lammps_")[1]] = force_lst # Store the forces into the dictionary
        ID_dict[data[s].split("data_lammps_")[1]] = ID_lst       # Store the atom IDs into the dictionary


        if res_print and not polar:                                            # Write the energy_table.txt
            write_eng_table(lmp, h, ff.get_name().lower(), data[s])            # ... only if not a polar calculations
            write_chg_table(lmp, m, charg_dict, train.get_charge(), data[s])
            write_force_table(lmp, n, force_dict, train.get_force(), data[s])

        lmp.close()

    if res_print:
        h.close()
        m.close()
        n.close()
        logging.info('...energy_table.txt success\n...charge_table.txt success\n...force_table.txt success')

    # ---------------------------------------------------------------------------------------------

    lammp_lst, error_lst, train_lst = compute_error(train, energ_dict, charg_dict, force_dict,ID_dict)    # Compure the error

    sum_error = 0
    for key in error_lst:
        sum_error += sum(error_lst[key])

    write_history(param, sum_error)                              # Write a history file

    if res_print:                                                # Write results to fort.99
        write_fort_99(train, lammp_lst, train_lst, error_lst)

    if debugging:
        logging.debug("Dictionaries...\n...train_energy: %s\n...train_charge: %s\n...train_force: %s\n...energ_dict: %s\n...charg_dict: %s\n...force_dict: %s\n...exiting optimizer" % (str(train.get_energy()), str(train.get_charge()), str(train.get_force()), str(energ_dict), str(charg_dict), str(force_dict)))
    

    return sum_error                                             # Return sum of output energies (not too sure about the error handling here)



def run(ffield, params, train, geo):
    ''' An optimizer for LAMMPS files. '''

    #Create a log file for the optimizer
    logging.basicConfig(filename='log.opt', filemode='w', level=logging.INFO)
    lvl_set = {0:'NOTSET',10:'DEBUG',20:'INFO',30:'WARNING',40:'ERROR',50:'CRITICAL'}
    lvl = logging.getLogger().getEffectiveLevel()
    if lvl in lvl_set:
        logging.info('Logging level - string: %s, integer: %d' % (lvl_set[lvl], lvl))
    else:
        logging.info('Logging level not set.')

    # --------------------------------------------------------------------------------------

    train_read = trainset.Trainset(train)                        # Read in the trainset file
    ff_read    = forcefield.ForceField(ffield)                   # Read in the force field using an object
    param_read = param.Param(params)                             # Read in the param file

    #if geo != None:                                              # Create data files if geo file supplied with proper valence
    #    if ff_read.get_optimizer()['valence'] == '0':
    #        geotodatalammps.data_file_complete(geo)
    #    elif ff_Read.get_optimizer()['valence'] == '1':
    #        geotodatalammps.data_file_simple(geo)

    if geo != None:                                              # Create data files if geo file supplied
        name = ff_read.get_name()
        if name == 'uffqm':
            geotodatalammps.data_file_complete(geo)
            #geotodatalammps.data_file_simple(geo)
        elif name == 'rexpon':
            geotodatalammps.data_file_simple(geo)
        elif name == 'morse':
            geotodatalammps.data_file_complete(geo)
        elif name == 'ljcut':
            geotodatalammps.data_file_simple(geo)
        elif name == 'exryd':
            geotodatalammps.data_file_simple(geo)
        elif name == 'exrydc':
            #geotodatalammps.data_file_libattery(geo)
            geotodatalammps.data_file_simple(geo)
        elif name == 'pprlg':
            #geotodatalammps.data_file_libattery(geo)
            geotodatalammps.data_file_simple(geo)
        else:
            logging.error('ERROR: no geo file supplied, and neither uffqm, rexpon, morse, exryd, or lj/cut selected.')
            print 'ERROR: no geo file supplied and neither uffqm, rexpon, morse, exryd, or lj/cut selected!'
            sys.exit()

    data_read  = readData()                                      # Read in the data files

    logging.debug('%s' % str(ff_read.get_ffield()))

    # --------------------------------------------------------------------------------------

    init = []                                                    # Pull the initial values of the parameters from the force field
    try:
        for k in range(len(param_read.get_bounds())):
            init.append(float(ff_read.get_ffield()[param_read.get_bounds()[k][0]][param_read.get_bounds()[k][1] - 1][param_read.get_bounds()[k][2] - 1]))
    except Exception as e:
        print "Error reading param file: ", e
        logging.error("Error reading param file: %s" % e)
        sys.exit()
    logging.info('Initial parameters read.')

    # --------------------------------------------------------------------------------------

    to_write = ""
    for i in range(len(init)):
        to_write += 'param[%d]\t' % i
    to_write += 'error\n'
    g = open('history', 'w')                                     # Set up the history file
    g.write(to_write)
    g.close()
    logging.info('history file created.')

    # --------------------------------------------------------------------------------------

    method    = ff_read.get_optimizer()['method']
    polar     = False if ff_read.get_optimizer()['polar'] == '0' else True
    algorithm = ff_read.get_optimizer()['algorithm']
    finish    = ff_read.get_optimizer()['finish']

    if method == '0':                                   # Run without optimization
        logging.info('Optimization flag turned off.')

        res = optimizer(init, data_read, train_read, ff_read, param_read.get_bounds(), polar, res_print=True, debugging=True, disp_screen=True)

    elif method == '1':                                 # Call optimizer on ffield parameters

        logging.info('Optimization flag turned on, preparing to optimize.')

        bnds = ()                                                # Get the bounds from params
        for i in range(len(param_read.get_bounds())):
            bnds = bnds + ((param_read.get_bounds()[i][-2], param_read.get_bounds()[i][-1]),)

        con = ()                                                 # Set constraints for parameters
        for i in range(len(param_read.get_constraints())):
            con = con + ({'type':'eq','fun':lambda x:x[param_read.get_constraints()[i][0]]-x[param_read.get_constraints()[i][1]]},)

        # Run the optimization
        start_time = time.time()
        if algorithm == 'slsqp':
            res = minimize(optimizer, init, args=(data_read, train_read, ff_read, param_read.get_bounds(), polar), method='SLSQP', bounds=bnds, constraints=con)
        elif algorithm == 'nelder-mead':
            res = minimize(optimizer, init, args=(data_read, train_read, ff_read, param_read.get_bounds(), True), method='Nelder-Mead', options={'eps':1e-9})#, options={'fatol':0.01})
        else:
            print 'Error: Algorithm not supplied. Options \'None\', \'SLSQP\', and \'Nelder-Mead\''
            logging.error('Error: Algorithm not supplied. Options \'None\', \'SLSQP\', and \'Nelder-Mead\'')
            sys.exit()
        final_time = time.time()

        write_out(ff_read, res.x, param_read.get_bounds())       # Create the fort.4 result force field

        # Create a fort.99 type file using optimized values
        optimizer(res.x, data_read, train_read, ff_read, param_read.get_bounds(), polar, res_print=True, debugging=False)

        logging.info("Optimal parameters: %s" % res)
        logging.info('Time of optimization: %0.3f' % (final_time-start_time))

    elif method == '2':                                          # Perform a map
        logging.info('Optimization flag turned on, preparing to optimize.')

        rrange = ()                                              # Get the bounds from params
        for i in range(len(param_read.get_bounds())):
            lower  = param_read.get_bounds()[i][-3]
            upper  = param_read.get_bounds()[i][-2]
            steps  = param_read.get_bounds()[i][-1]
            rrange = rrange + ((lower,upper,steps),)

        # Run the optimization
        start_time = time.time()
        if finish == 'none':
            res = brute(optimizer, rrange, args=(data_read, train_read, ff_read, param_read.get_bounds(), polar), finish=None)
        elif finish == 'slsqp':
            res = brute(optimizer, rrange, args=(data_read, train_read, ff_read, param_read.get_bounds(), polar), finish='SLSQP')
        elif finish == 'nelder-mead':
            res = brute(optimizer, rrange, args=(data_read, train_read, ff_read, param_read.get_bounds(), polar), finish='Nelder-Mead')
        else:
            print 'Error: Finish not supplied. Options \'None\', \'SLSQP\', and \'Nelder-Mead\''
            logging.error('Error: Finish not supplied. Options \'None\', \'SLSQP\', and \'Nelder-Mead\'')
            sys.exit()
        final_time = time.time()

        write_out(ff_read, res.x, param_read.get_bounds())       # Create the fort.4 result force field

        # Create a fort.99 type file using optimized values
        optimizer(res.x, data_read, train_read, ff_read, param_read.get_bounds(), polar, res_print=True, debugging=False)

        logging.info("Optimal parameters: %s" % res)

    else:
        print 'ERROR: No valid method supplied'
        logging.error('ERROR: No valid method supplied')
        sys.exit()

    return res


