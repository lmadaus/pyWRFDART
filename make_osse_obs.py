#!/usr/bin/env python
from __future__ import division, print_function
import os, sys
import re
from netCDF4 import Dataset
sys.path.append('/glade/p/work/lmadaus/cm1/kffc_ensemble')
from ens_dart_param import *

"""
Script to operate on output from a "truth" run, found in the members directory,
to extract an observation sequence at the specified time.
"""

GENERATE_IDEAL_OBS = True


def main():

    build_obs_structure(0)

    #if GENERATE_IDEAL_OBS:
    #    generate_ideal_obs(intime=0)

def build_obs_structure(intime, gridspace=8):
    """ Function to specify what variables we want and how dense they should be """
    from make_namelist_dart import set_namelist_sectors, write_namelist
    from write_cm1_namelist import set_namelist_defaults
    # Generate the DART namelist structure so we
    # can query (and modify) it
    dartnml = set_namelist_sectors()
    cm1nml = set_namelist_defaults()
    # Parse out the observations to assimilate
    use_obs = [s.strip()[1:-1] for s in dartnml['obs_kind']['assimilate_these_obs_types'].split(',')]
    # I'm just testing this for now
    use_obs = [use_obs[1]]
    print(use_obs)

    # Figure out the grid structure
    dx = cm1nml['param1']['dx']
    dy = cm1nml['param1']['dy']
    nx = cm1nml['param0']['nx']
    ny = cm1nml['param0']['ny']

    # Now figure out how many obs we'll need
    gridspace *= 1000
    obx = range(0,int(nx*dx),gridspace)
    oby = range(0,int(ny*dx),gridspace)
    total_obs = len(obx)*len(oby)*len(use_obs)
    print("Total number of obs:", total_obs)
    
    




def generate_ideal_obs(intime):
    """
    Run the DART sequence to generate "truth" observations from a
    run within the truth directory
    """
    # Multiply intime by 60 for seconds
    intime *= 60
    intime = int(intime)

    # Check to see that the truth run is where it's supposed to be
    truthdir = os.path.join(dir_mems, 'truth')
    if not os.path.exists(truthdir):
        raise ProceduralError('Unable to find directory{:s}/truth'.format(dir_mems),\
                                'generate_ideal_obs')
    # Search the restart files for the appropriate time
    rst_filenames = [f for f in os.listdir(truthdir) if 'rst' in f and f.endswith('.nc')]
    if len(rst_filenames) < 1:
        raise ProceduralError('Unable to find any rst files in {:s}'.format(truthdir),\
                                'generate_ideal_obs')
    
    # Match these files to the intime
    rst_files = [Dataset(f,'r') for f in rst_filenames]
    file_times = [int(r.variables['time'][0]) for r in rst_files]
    for f in rst_files:
        f.close()
    if intime not in file_times:
        raise ProceduralError('Could not find time {:d} in any rst files in {:s}'.format(intime,truthdir),\
                                'generate_ideal_obs')
    # Get this file
    rst_file = rst_filenames[file_times.index(intime)]
    print("Found restart file at time {:d}: {:s}".format(intime, rst_file))

    # Write the dart namelist for this time
    os.system('./write_namelist_dart.py -d {:d}'.format(intime/60))
    os.system('cp input.nml {:s}/'.format(dir_obs))

    # Go into the obs directory
    os.chdir(dir_obs)

    # Link in the executables we may need and run in order
    exe_sequence = ['preprocess','create_obs_sequence',\
        'create_fixed_network_sequence','perfect_model_obs']
    for exe in exe_sequence:
        if not os.path.exists(exe):
            os.system('ln -sf {:s} .'.format(os.path.join(dir_src_dart,exe)))

        # Special instructions for each here
        if exe == 'create_obs_sequence':
            # Here goes the code that tells
            # which obs to generate

            pass
        elif exe == 'create_fixed_network_sequence':
            # Here we have more code to tell
            # create_fixed_network_sequence about
            # the time for the obs
            pass

        os.system('./{:s}'.format(exe))





class ProceduralError(Exception):
    """
    A simple class to throw a text exception
    """
    def __init__(self, text, function):
        self.text = text
        self.function = function
    def __str__(self):
        return repr(self.text)
            


if __name__ == '__main__':
    main()
