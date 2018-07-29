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


import sys                                                       # Command-line inputs
from src import optimizer                                        # The optimizer has been moved to a separate file
from optparse import OptionParser                                # Command line arguments

def main():

    helpString ="""
        An optimizer for LAMMPS files, see README for more information
        Usage: python optimizer.py -f [ffield] -p [params] -t [trainset] -g [geo]
    """

    args = sys.argv[1:]
    numArgs = len(args)
    if numArgs == 0:
        print helpString
        #sys.exit()


    # parse the command line
    parser = OptionParser()
    parser.add_option("-f",dest="F",help="ffield and control file",default="ffield")
    parser.add_option("-p",dest="P",help="parameter list to optimize ",default="params")
    parser.add_option("-t",dest="T",help="trainset file",default="trainset")
    parser.add_option("-g",dest="G",help="geometry file with bgf format")
    (options,args) = parser.parse_args()

    # set input files
    ffield=options.F
    params=options.P
    train=options.T
    geo=options.G

    # perform optimization
    results = optimizer.run(ffield, params, train, geo)

    print results

if __name__ == "__main__":
    main()

