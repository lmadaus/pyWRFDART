#!/usr/bin/env python
from __future__ import division, print_function
import os, sys
import re
from datetime import datetime, timedelta
from netCDF4 import Dataset
sys.path.append('/home/disk/pvort/nobackup/lmadaus/cm1/DOMAINS/kdvn_ensemble')
from ens_dart_param import *

"""
Script to operate on output from a "truth" run, found in the members directory,
to extract an observation sequence at the specified time.
"""

GENERATE_IDEAL_OBS = True


error_var = {'LAND_SFC_TEMPERATURE' : 1.0,
             'LAND_SFC_U_WIND_COMPONENT' : 1.0,
             'LAND_SFC_V_WIND_COMPONENT' : 1.0,
             'LAND_SFC_PRESSURE' : 100.0,
             'LAND_SFC_SPECIFIC_HUMIDITY' : 0.001}

def main():

    #build_obs_structure(0)

    if GENERATE_IDEAL_OBS:
        generate_ideal_obs(intime=120)

def build_obs_structure(intime, rst_file, gridspace=16):
    """ Function to specify what variables we want and how dense they should be """
    #from make_namelist_dart import set_namelist_sectors, write_namelist
    #from write_cm1_namelist import set_namelist_defaults
    from namelist_utils import read_namelist, write_dart_namelist

    # Generate the DART namelist structure so we
    # can query (and modify) it
    dartnml = read_namelist('input.nml')
    cm1nml = read_namelist('namelist.input')
    # Build the epoch from the cm1 namelist
    param11 = cm1nml['param11']
    startdate = datetime(param11['year'], param11['month'],\
                     param11['day'], param11['hour'],\
                     param11['minute'], param11['second'])
    
    # Parse out the observations to assimilate
    use_obs = dartnml['obs_kind_nml']['assimilate_these_obs_types']
    # I'm just testing this for now
    #use_obs = use_obs[0:2]
    use_obs = ['LAND_SFC_PRESSURE']
    print(use_obs)
    # Put time as datetime
    intime = startdate + timedelta(seconds=intime)
    epoch = datetime(1601,1,1,0)
    print("Indate:", intime)

    # Figure out the grid structure
    dx = cm1nml['param1']['dx']
    dy = cm1nml['param1']['dy']
    nx = cm1nml['param0']['nx']
    ny = cm1nml['param0']['ny']

    # Now figure out how many obs we'll need
    gridspace *= 1000
    obx = range(int(dx/2.),int(nx*dx),gridspace)
    oby = range(int(dy/2.),int(ny*dx),gridspace)
    #obx = [16000]
    #oby = [16000]
    total_obs = len(obx)*len(oby)*len(use_obs)
    print("Total number of obs:", total_obs)
   
    # Write a new dart namelist for the current time
    delt = intime - epoch
    #dartnml['perfect_model_obs_nml']['first_obs_days'] = str(delt.days)
    #dartnml['perfect_model_obs_nml']['first_obs_seconds'] = str(delt.seconds)
    #dartnml['perfect_model_obs_nml']['last_obs_days'] = str(delt.days)
    #dartnml['perfect_model_obs_nml']['last_obs_seconds'] = str(delt.seconds)
    #dartnml['perfect_model_obs_nml']['init_time_days'] = str(delt.days)
    #dartnml['perfect_model_obs_nml']['init_time_seconds'] = str(delt.seconds)
    dartnml['perfect_model_obs_nml']['obs_seq_in_file_name'] = 'obs_seq.in'
    dartnml['perfect_model_obs_nml']['restart_in_file_name'] = rst_file

    # Make sure we have the right io pattern
    dartnml['io_filenames_nml']['restart_in_stub'] = 'notinuse'
    dartnml['io_filenames_nml']['overwrite_input'] = False
    dartnml['io_filenames_nml']['rpointer'] = True
    dartnml['io_filenames_nml']['rpointer_file'] = 'input_filelist.txt'

    # Write the pointer file
    with open('input_filelist.txt', 'w') as pointerfile:
        pointerfile.write(rst_file)
    os.system('mv input_filelist.txt {:s}/'.format(dir_obs))

    # Write the modified namelist
    write_dart_namelist(dartnml)

    # Set up the input file
    with open('obs_seq_input.txt','w') as infile:
        infile.write(str(total_obs)+'\n') # Total num obs
        infile.write('0\n') # 0 copies
        infile.write('0\n') # 0 QC values
        for obtype in use_obs:
            for x in obx:
                for y in oby:
                    infile.write('0\n')
                    infile.write(obtype+'\n') # Identify ob type
                    infile.write('0\n') # Specify location
                    infile.write(str(x)+'\n') # X coordinate
                    infile.write(str(y)+'\n') # Y coordinate
                    infile.write('0\n') # Z coordinate
                    infile.write('{:d} {:d} {:d} {:d} {:d} {:d}\n'.format(\
                        intime.year, intime.month,\
                        intime.day, intime.hour,\
                        intime.minute, intime.second))
                    infile.write(str(error_var[obtype])+'\n') # Error variance
        infile.write('obs_seq.in\n')
    # Move these files to the obs dir
    if os.path.exists(os.path.join(dir_obs, 'input.nml')):
        os.system('rm -f {:s}'.format(os.path.join(dir_obs, 'input.nml')))
    if os.path.exists(os.path.join(dir_obs, 'obs_seq_input.txt')):
        os.system('rm -f {:s}'.format(os.path.join(dir_obs, 'obs_seq_input.txt')))
    os.system('mv obs_seq_input.txt {:s}'.format(dir_obs))
    os.system('cp input.nml {:s}'.format(dir_obs))

    return None 

def generate_ideal_obs(intime):
    """
    Run the DART sequence to generate "truth" observations from a
    run within the truth directory
    """
    # Multiply intime by 60 for seconds
    intime *= 60
    intime = int(intime)

    # Check to see that the truth run is where it's supposed to be
    truthdir = os.path.join(dir_members, 'truth')
    if not os.path.exists(truthdir):
        raise ProceduralError('Unable to find directory{:s}/truth'.format(dir_members),\
                                'generate_ideal_obs')
    # Search the restart files for the appropriate time
    rst_filenames = [f for f in os.listdir(truthdir) if 'rst' in f and f.endswith('.nc')]
    if len(rst_filenames) < 1:
        raise ProceduralError('Unable to find any rst files in {:s}'.format(truthdir),\
                                'generate_ideal_obs')
    
    # Match these files to the intime
    rst_files = [Dataset(os.path.join(truthdir,f),'r') for f in rst_filenames]
    file_times = [int(r.variables['time'][0]) for r in rst_files]
    for f in rst_files:
        f.close()
    if intime not in file_times:
        raise ProceduralError('Could not find time {:d} in any rst files in {:s}'.format(intime,truthdir),\
                                'generate_ideal_obs')
    # Get this file
    rst_file = rst_filenames[file_times.index(intime)]
    os.system('rm -f {:s}/cm1out_rst*.nc'.format(dir_obs))
    os.system('ln -sf {:s}/{:s} {:s}/{:s}'.format(truthdir, rst_file, dir_obs, rst_file))
    print("Found restart file at time {:d}: {:s}".format(intime, rst_file))
    
    # Build the obs structure file based on this          
    build_obs_structure(intime, rst_file)

    #os.system('./write_namelist_dart.py -d {:d}'.format(intime/60))
    #os.system('cp input.nml {:s}/'.format(dir_obs))
    
    
    
    
    # Go into the obs directory
    os.chdir(dir_obs)

    # Cleanup
    os.system('rm -f True_State.nc obs_seq.in obs_seq.out')
    if not os.path.exists('cm1out_rst_000001.nc'):
        os.system('ln -sf {:s} .'.format(os.path.join(truthdir, 'cm1out_rst_000001.nc')))
    if not os.path.exists('namelist.input'):
        os.system('cp {:s}/namelist.input .'.format(dir_dom))
    # Link in the executables we may need and run in order
    exe_sequence = ['create_obs_sequence','perfect_model_obs']
    for exe in exe_sequence:
        if not os.path.exists(exe):
            os.system('ln -sf {:s} .'.format(os.path.join(dir_src_dart,exe)))

        # Special instructions for each here
        if exe == 'create_obs_sequence':
            # Here goes the code that tells
            # which obs to generate
            os.system('./{:s} < obs_seq_input.txt'.format(exe))

        else:
            os.system('./{:s}'.format(exe))


    # Check to be sure it worked
    if not os.path.exists('obs_seq.out'):
        raise ProceduralError("perfect_model_obs sequence failed! No obs_seq.out file.",'')
    else:
        print("Success!  Made obs_seq.out file.")
        os.system('mv obs_seq.out {:d}_obs_seq.prior'.format(intime))



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
