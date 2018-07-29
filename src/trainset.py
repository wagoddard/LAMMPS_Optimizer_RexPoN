# Class for the trainset file

import sys
import logging

class Trainset(object):
    ''' Read in all of the energies from trainset.in '''


    def __init__(self, train):
        """ Assume train is the trainset file to read """

        self.energy     = []                                     # The energy of a given structure, relative to other structures
        self.charge     = []                                     # The charge of an atom in a structure
        self.force      = []                                     # The atomic forces on (x,y,z,t) directions of particular atom

        self.geometry   = []                                     # Eventually
        self.cell_param = []                                     # Eventually
        self.pressure   = []                                     # Eventually


        self.headers = ["ENERGY", "CHARGE", "FORCES", "GEOMETRY", "CELL PARAMETER", "PRESSURE"]
        self.ffield  = {}


        logging.info("Reading trainset file...")

        try:
            f = open(train, 'r')
            lines = f.readlines()
            f.close()
        except IOError:
            logging.error("Trainset does not exist")
            sys.exit()


        type = None
        for i in range(len(lines)):
            head = lines[i].split()[0]

            # All sections stored as lists
            if head == 'ENERGY':
                j = 1
                while i+j<len(lines) and lines[i + j].split()[0] != 'ENDENERGY':
                    if lines[i + j][0] != '#':
                        self.energy.append([])
                        for k in range(len(lines[i + j].split())):
                            self.energy[-1].append(lines[i + j].split()[k])
                    j += 1
                if lines[i + j].split()[0] != 'ENDENERGY':
                    print 'ERROR: Section ENERGY not ended.'
                    logging.error('Section ENERGY not ended.')
                    sys.exit()
            if head == 'CHARGE':
                j = 1
                while i+j<len(lines) and lines[i + j].split()[0] != 'ENDCHARGE':
                    if lines[i + j][0] != '#':
                        self.charge.append([])
                        for k in range(len(lines[i + j].split())):
                            self.charge[-1].append(lines[i + j].split()[k])
                    j += 1
                if lines[i + j].split()[0] != 'ENDCHARGE':
                    print 'ERROR: Section CHARGE not ended.'
                    logging.error('Section CHARGE not ended.')
                    sys.exit()
            if head == 'FORCES':
                j = 1
                while i+j<len(lines) and lines[i + j].split()[0] != 'ENDFORCES':
                    if lines[i + j][0] != '#':
                        self.force.append([])
                        for k in range(len(lines[i + j].split())):
                            self.force[-1].append(lines[i + j].split()[k])
                    j += 1
                if lines[i + j].split()[0] != 'ENDFORCES':
                    print 'ERROR: Section FORCES not ended.'
                    logging.error('Section FORCES not ended.')
                    sys.exit()

        logging.debug('readTrainset(): %s' % str(self.energy))
        logging.info("...reading success.")


    def get_energy(self):
        """ Return the energy array """
        return self.energy

    def get_charge(self):
        """ Return the charge array """
        return self.charge

    def get_force(self):
        """ Return the force array """
        return self.force

    def get_geometry(self):
        """ Return the geometry array """
        return self.geometry

    def get_cell_param(self):
        """ Return the cell parameter array """
        return self.cell_param

    def get_pressure(self):
        """ Return the pressure array """
        return self.pressure





