# ForceField class that should be used while reading in ffield in the optimizer

import sys
import logging

class ForceField(object):
    """ A force field object that consists of pair_styles, mol. weights, and pair_types """


    def __init__(self, ffield):
        """ Assume ffield is the force field file to read """

        # Read the force field file to be parsed
        f = open(ffield, 'r')
        lines = f.readlines()
        f.close()

        self.headers = ["GENERAL", "ATOMS", "BONDS", "ANGLES",   # The expected sections of ffield 
            "DIHEDRALS", "HYDROGEN_BOND", "VDW", "ELECTRO", "CONTROL", "OPTIMIZER", "FFIELD"]
        self.ffield  = {}                                        # The force field is stored in this dictionary
        self.name    = None                                      # Type of force field (e.g. rexpon, uffqm, ...)
        #self.opt     = 1                                         # Optimization flag is turned on by default
        self.optimizer = {'method':'0','polar':'0','algorithm':'slsqp','finish':'none','valence':'0'}

        for i in self.headers:
            self.ffield[i] = []

        # Get the forcefield name first
        for i in range(len(lines)):
            head = lines[i].split()[0]
            if head == "FFIELD":
                self.name = lines[i + 1].split()[0]
                break
        if self.name == None:
            logging.error('ERROR: FFIELD type not defined!\n')
            print 'ERROR: FFIELD type not defined!\n'
            sys.exit()


        # Default setting for CONTROL file
        if self.name == 'uffqm':
            self.set_ctrl_uffqm()
        elif self.name == 'rexpon':
            self.set_ctrl_rexpon()
        elif self.name == 'opls':                # TEMPORARY
            self.set_ctrl_uffqm()
        elif self.name == 'morse':
            self.set_ctrl_morse()
        elif self.name == 'exryd':
            self.set_ctrl_exryd()
        elif self.name == 'exrydc':
            self.set_ctrl_exrydc()
        elif self.name == 'ljcut':
            self.set_ctrl_ljcut()
        elif self.name == 'pprlg':
            self.set_ctrl_ljcut()

        type = None
        for i in range(len(lines)):
            head = lines[i].split()[0]
  
            # Optimization flag
            #if head == "OPT":
                ##if lines[i + 1].split()[0] == '0':
                ##    self.opt = False
                #self.opt = int(lines[i+1].split()[0])

            # All sections, excluding CONTROL, stored as lists
            if head in self.headers and head not in ['CONTROL','OPTIMIZER','FFIELD']:
                self.ffield[head] = []
                j = 1
                while i+j < len(lines) and lines[i + j].split()[0] not in self.headers:
                    if lines[i + j].split()[0][0] != '#':
                        self.ffield[head].append([])
                        for k in lines[i + j].split():
                            self.ffield[head][-1].append(k)
                    j += 1
 
            # CONTROL section stored as dictionary
            if head == 'CONTROL':
                j = 1
                while i+j < len(lines) and lines[i + j].split()[0] not in self.headers:
                    if lines[i + j].split()[0][0] != '#':
                        self.ffield['CONTROL'][lines[i+j].split()[0]] = ' '.join(lines[i+j].split()[1:])
                    j += 1

            # OPTIMIZER section stored as dictionary
            if head == 'OPTIMIZER':
                j = 1
                while i+j < len(lines) and lines[i + j].split()[0] not in self.headers:
                    if lines[i + j].split()[0][0] != '#':
                        if lines[i + j].split()[0] in self.optimizer:
                            self.optimizer[lines[i + j].split()[0]] = lines[i + j].split()[1].lower()
                        else:
                            print 'OPTIMIZER variable %s not defined - continuing program' % lines[i + j].split()[0]
                            logging.warning('OPTIMIZER variable %s not defined - continuing program' % lines[i + j].split()[0])
                    j += 1

    def set_ctrl_uffqm(self):
        ''' Set the defaults for the control section assuming UFFQM. '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'
        self.ffield['CONTROL']['special_bonds'] = 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'

    def set_ctrl_rexpon(self):
        ''' Set the defaults for the control section assuming rexpon. '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'


    def set_ctrl_morse(self):
        ''' Set the defaults for the control section assuming MORSE. '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'
        self.ffield['CONTROL']['special_bonds'] = 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'


    def set_ctrl_exryd(self):
        ''' Set the defaults for the control section assuming EXRYD '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'
        self.ffield['CONTROL']['special_bonds'] = 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'

    def set_ctrl_exrydc(self):
        ''' Set the defaults for the control section assuming EXRYDC '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'
        self.ffield['CONTROL']['special_bonds'] = 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'


    def set_ctrl_ljcut(self):
        ''' Set the defaults for the control section assuming LJCUT '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'
        self.ffield['CONTROL']['special_bonds'] = 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'

    def set_ctrl_pprlg(self):
        ''' Set the defaults for the control section assuming PPRLG '''

        self.ffield['CONTROL'] = {}
        self.ffield['CONTROL']['units'] = 'real'
        self.ffield['CONTROL']['atom_style'] = 'pqeq'
        self.ffield['CONTROL']['boundary'] = 'f f f'
        self.ffield['CONTROL']['special_bonds'] = 'lj 0.0 0.0 1.0 coul 1.0 1.0 1.0'

    def get_ffield(self):
        ''' Return the array of pair styles '''
        return self.ffield

    def get_name(self):
        ''' Return the name/type of the force field '''
        return self.name

    #def get_opt(self):
    #    ''' Return the optimization variable '''
    #    return self.opt

    def get_optimizer(self):
        ''' Return the optimizer parameters '''
        return self.optimizer

    def set_ffield(self, ff):
        ''' Set a new ff '''
        self.ffield = ff


