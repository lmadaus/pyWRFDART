#!/usr/bin/env python
import os, sys
os.chdir('/glade/u/home/lmadaus/DOMAINS/relampago')
sys.path.append('/glade/u/home/lmadaus/DOMAINS/relampago')
from namelist_utils import read_namelist, write_namelist, update_time_wrf
from WRF_dart_param import *
from datetime import datetime, timedelta
from time import sleep

# Number of perturbations to generate
nperts = 250


def main():
    """
    Runs through the sequence to generate perturbations
    """

    setup_perturbation_making()
    subfile = 'make_pert.csh'
    make_wrfda_submit_script(subfile)
    # Now loop through
    for pnum in range(nperts):
        run_wrfda_for_perts(pnum, subfile)


def setup_perturbation_making():
    """
    Preliminary work for making pertubation generator work
    Writes namelists and copies an existing wrfinput_d01 file
    into the "fg" template file.  This file's time is currently
    expected to match the "date_end" time in the WRF_dart_param file.
    If it does not, will need to adjust so that the wrfvar18/analysis_date
    namelist option below matched the date of the wrfinput_d01/fg file.

    """
    # Get the test date as the start date of this experiment
    testdate = datetime.strptime(date_end, '%Y%m%d%H')

    # Read in the namelist
    nmld = read_namelist('TEMPLATE_namelist.input')
    # Update the namelist with the relevant wrfvar sections
    use_keys = [k for k in wrf_namelist.keys() if k.startswith('wrfvar')]
    for k in use_keys:
        nmld[k] = wrf_namelist[k]

    # Set the right time
    nmld['wrfvar18']['analysis_date'] = testdate.strftime('%Y-%m-%d_%H:%M:%S.0000')
    # Set seeds here in case we want to use
    nmld['wrfvar11']['seed_array1'] = int(testdate.strftime('%Y%m%d%H'))
    nmld['wrfvar11']['seed_array2'] = 1 
    # Setting this to False ignores the seeds above and just does random
    # perts every time
    nmld['wrfvar5']['put_rand_seed'] = False
    # CV3 uses stored fixed background error covariances in wrfda
    nmld['wrfvar7']['cv_options'] = 3
    # Scaling for perturbations (variance, horizontal length, vertical length)
    nmld['wrfvar7']['as1'] = [0.15,1.0,0.75] # Streamfunction scaling (variance, horizontal, vertical)
    nmld['wrfvar7']['as2'] = [0.15,1.0,0.75] # Unbalanced Velocity potential
    nmld['wrfvar7']['as3'] = [0.15,1.0,0.75] # Unbalanced temperature
    nmld['wrfvar7']['as4'] = [0.15,1.0,0.75] # Pseudo relative humidity
    nmld['wrfvar7']['as5'] = [0.15,1.0,0.75] # Unbalanced surface pressure
    # This tells us just to perturb the state...no assimilation
    nmld['wrfvar17']['analysis_type'] = 'RANDOMCV'

    # Update the times and set the max_dom to be 1 (only perturbing outer domain)
    nmld = update_time_wrf(nmld, testdate, testdate)
    nmld['domains']['max_dom'] = 1
    # Write out the namelist
    write_namelist(nmld, 'namelist.input')
    # Copy in our template wrfinput_d01 file to the required "fg" file
    # wrfdart will be looking for
    os.system('cp wrfinput_d01 fg')
    return



def make_wrfda_submit_script(subfile_name):
    """
    Makes a queue submission script for running da_wrfvar.exe
    """
    with open(subfile_name, 'w') as outfile:
        outfile.write('#!/bin/csh\n')
        outfile.write('#PBS -N {:s}\n'.format(subfile_name[:-4]))
        outfile.write('#PBS -A {:s}\n'.format(NCAR_GAU_ACCOUNT))
        outfile.write('#PBS -l walltime=00:02:00\n') # Should be fast...you can change time here
        outfile.write('#PBS -q {:s}\n'.format(ADVANCE_QUEUE_MEMBER))
        #outfile.write('#PBS -j oe\n')
        #outfile.write('#PBS -m abe\n')
        #outfile.write('#PBS -M {:s}@ucar.edu\n'.format(getuser()))
        outfile.write('#PBS -l select=1:ncpus={:d}:mpiprocs={:d}\n'.format(numprocs_per_node, numprocs_per_node))
        outfile.write("#==================================================================\n")
        outfile.write("\n")
        outfile.write('# Change into appropriate directory\n')
        outfile.write('cd {:s}\n'.format(dir_wrf_dom))
        outfile.write('# Run da_wrfvar.exe\n')
        outfile.write('mpiexec_mpt ./da_wrfvar.exe\n')

    return







def run_wrfda_for_perts(pertnum, subfile_name):
    """
    Uses pre-defined queue submission script to run wrfda,
    take difference of perturbed and original wrfinput_d01,
    and storing those perturbations in the pert bank directory
    """

    print("Making perturbation member {:d}".format(pertnum+1))
    # Make sure this file does not exist before starting
    if os.path.exists('wrfvar_output'):
        os.system('rm wrfvar_output')
    # Submit the script to run da_wrfvar
    os.system('qsub < {:s}'.format(subfile_name))

    # Wait until the output file has been generated
    while not os.path.exists('wrfvar_output'):
        sleep(5)
    # Clean up all this other output
    os.system('rm buddy_check cost_fn grad_fn jo check_max_iv qcstat_conv_01')
    os.system('rm gts_omb_*')
    os.system('rm rej_obs_co*')
    os.system('rm unpert_obs*')
    os.system('rm new_diff.nc')

    # ncdiff will give the perturbations
    os.system('ncdiff wrfvar_output wrfinput_d01 new_diff.nc')

    # Only need to retain the needed variables
    os.system('ncks -v T,U,V,QVAPOR,MU new_diff.nc {:s}/pert_bank_mem_{:d}.nc'.format(pert_bank_path,pertnum+1))

if __name__ == '__main__':
    main()
