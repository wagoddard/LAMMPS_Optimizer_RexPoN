# Class for the parameter file

import sys
import logging

class Param(object):
    ''' Read in all of the energies from trainset.in '''


    def __init__(self, param):
        """ Assume train is the trainset file to read """

        self.bounds     = []                                     # The bounds for a given variable
        self.cons       = []                                     # Constraints relating many variables

        self.headers = ["BOUNDS","CONSTRAINTS"]

        logging.info("Reading param file...")

        try:
            f = open(param, 'r')
            lines = f.readlines()
            f.close()
        except IOError:
            logging.error("Params does not exist")
            sys.exit()

        # Read the parameter file into a list of lists
        # Assuming the propert syntax of pair_style [int], pair_type [int], param [int], range_low [float], range_high [float]
        try:
            for i in range(len(lines)):
                head = lines[i].split()[0]
  
                # All sections stored as lists
                if head == "BOUNDS":
                    j = 1
                    while i+j < len(lines) and lines[i + j].split()[0] not in self.headers:
                        if lines[i + j].split()[0][0] != '#' and len(lines[i+j].split()) == 6:
                            self.bounds.append([lines[i+j].split()[1], int(lines[i+j].split()[2]), int(lines[i+j].split()[3]), float(lines[i+j].split()[4]), float(lines[i+j].split()[5])])
                        elif lines[i + j].split()[0][0] != '#' and len(lines[i+j].split()) == 7:
                            self.bounds.append([lines[i+j].split()[1], int(lines[i+j].split()[2]), int(lines[i+j].split()[3]), float(lines[i+j].split()[4]), float(lines[i+j].split()[5]), float(lines[i+j].split()[6])])
                        j += 1
                elif head == "CONSTRAINTS":
                    j = 1
                    while i+j < len(lines) and lines[i + j].split()[0] not in self.headers:
                        if lines[i + j].split()[0][0] != '#':
                            self.cons.append([int(lines[i+j].split()[1])-1, int(lines[i+j].split()[2])-1])
                        j += 1
        except Exception as e:
            logging.error("Error reading params: %s" % str(e))
            print "Error reading params, see log.opt"
            sys.exit()


        logging.debug('bounds: %s' % str(self.bounds))
        logging.debug('cons: %s' % str(self.cons))
        logging.info("...reading success.")


    def get_bounds(self):
        """ Return the bounds array """
        return self.bounds

    def get_constraints(self):
        """ Return the constraints array """
        return self.cons 

