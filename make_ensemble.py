#!/usr/bin/env python

import os, sys, getopt
# Get parameters from ens_dart_param and the default cm1 namelist
from ens_dart_param import *
from namelist_utils import read_namelist, write_namelist
cm1nml = read_namelist(os.path.join(dir_dom, 'namelist.input'))

"""
This is a quickly-reformatted version of the make_ensemble script
from version 1 of this code that makes it compatible with the new
namelist version.  This script will make an ensemble (in serial) with
a number of members as specified in the namelist variable "Ne" or manually
specified from the command line.  The script will copy in all of the necessary
files.  No perturbation of initial conditions is enabled at this time.
LEM 09/2015

"""


(opts,args) = getopt.getopt(sys.argv[1:],'n:')
for o,a in opts:
   if o == '-n':
      Ne = int(a)



def main():
    DOMDIR = dir_dom
    os.chdir(DOMDIR)
    """
    if os.path.isdir(DOMDIR+'/mems'):
        remove = raw_input('Remove contents of current ensemble directory (0 or 1)?')
        if remove:
            os.system('rm -rf {:s}/mems/*'.format(DOMDIR))
        else:
            exit(1)
    else:
        os.system('mkdir mems')  
    """
    print "############### MAKING ENSEMBLE ###############"
    print "############## NUMBER OF MEMS: %d #############" % Ne
    # Now populate the ensemble directory
    print "Populating ensemble directory (mems)"
    os.chdir(DOMDIR+'/mems')

    # Check our initialization method to see what we need for
    # each ensemble member
    if cm1nml['param2']['isnd'] == 7:
        # Will need an input sounding file.  Check for this
        if not os.path.exists('{:s}/input_sounding'.format(DOMDIR)):
            print("ERROR: CM1 namelist has isnd=7, but no input_sounding file found in {:s}".format(DOMDIR))
            exit(1)


    # Loop through each member
    print "######################################"
    print "Making ensemble member sub-directories"
    print "######################################"
    print ""
    for k in range(Ne):
    #for k in [0]:
        k = k+1 # Don't start at 0
        print "***** MEMBER {:d} *****".format(k)
        os.chdir(DOMDIR)
        # Run script to generate the namelist
        write_namelist(cm1nml, 'namelist.input')
        # That will be written as namelist.input
        # hold onto that while making a directory for it
        os.chdir('{:s}/mems'.format(DOMDIR))
        os.system('mkdir {:s}/mems/m{:d}'.format(DOMDIR,k))
        os.chdir('{:s}/mems/m{:d}'.format(DOMDIR,k))
        # Copy in the namelist
        os.system('cp {:s}/namelist.input namelist.input'.format(DOMDIR))
        # See if we need an input sounding file
        if cm1nml['param2']['isnd'] == 7:
            os.system('cp {:s}/input_sounding .'.format(DOMDIR))
            

        # Link in executable
        os.system('ln -sf {:s}/LANDUSE.TBL .'.format(dir_src_model))
        os.system('ln -sf {:s}/cm1.exe .'.format(dir_src_model))

        # Go back up a directory to start the next one
        os.chdir('%s/mems' % DOMDIR)

    print "### DONE MAKING DIRECTORIES ###"
    print ""

    os.chdir(DOMDIR)



if __name__ == '__main__':
    main()
