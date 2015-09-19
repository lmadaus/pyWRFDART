#!/usr/bin/env python

import os, sys
from datetime import datetime, timedelta
import re
from netCDF4 import Dataset
sys.path.append('/home/disk/pvort/nobackup/lmadaus/cm1/DOMAINS/kdnr_ensemble')
from ens_dart_param import *


"""
# New master script that will run all member-specific tasks on the queue
# Leaving the run of "filter" as the only remaining task

# Currently, this script is set up to do cycling without forecast runs 


# At present, supplying the member number and the start time comes from
# a dumped text file supplied by the master controls script


# Changing these variables will change which segments of the sequence
# to process.  This lets you fix problems with just one part of the sequence
# without having to run the whole thing.  Set all to True to run everything.
# Note that even if set to true, dart_to_cm1 and cm1_to_dart won't necessarily
# run if the starttime or endtime of the cycle is not in the list of
# assimilation times
"""

PRE_CLEAN     = True
DART_TO_MODEL = False 
RUN_MODEL     = True
POST_MODEL    = True
MODEL_TO_DART = False



def main():
    # MAIN PROCESSING SEQUENCE BEGINS HERE
    if PRE_CLEAN:
        preclean()

    # Get assimilation times so we know if this is starting post-assimilation,
    # running to an assimilation, or none of the above
    print "Getting assim times"
    assim_times = get_assim_times()
    print "ASSIM TIMES:",assim_times

    # Get the relevant information for this member
    print "Getting model info"
    memnum, start, end = get_basic_info() 
    print "START TIME:", start
    print "END TIME:", end
    # Compute cycle_len
    # Get any additional forecast length
    fcst_end = end
    if fcst_len < 0:
        fcst_end = exp_length
    else:
        fcst_end += fcst_len
    cycle_len = int(end - start)
    #forecast_len = int(fcst_end - end)

    # Check here if the flag_direct_netcdf_io is True.  If so, set the
    # DART_TO_WRF and WRF_TO_DART sequences to False
    if flag_direct_netcdf_io:
        DART_TO_MODEL = False
        MODEL_TO_DART = False

    # Check to see if we need to run dart_to_cm1
    if DART_TO_MODEL:
        print "Checking for dart_to_cm1"
        if start in assim_times:
            run_dart_to_cm1(memnum, startdate)


   
    # Run MODEL
    if RUN_MODEL:
        print "Running CM1"
        if start in assim_times:
            # Prepare cm1 to start from a restart
            #restart_name = 'cm1out_rst_000001.nc'
            cm1_prep(memnum,start,fcst_end)
        else:
            cm1_prep(memnum,start,fcst_end)
        run_cm1(memnum)

    # Do post-MODEL cleanup
    if POST_MODEL:
        print "Doing post-MODEL cleanup"
        #if end in assim_times:
        post_model_cleanup(memnum, start, end, fcst_end)
        #else:
        #    post_model_cleanup(memnum, start, end, fcst_end)

    # If we end at an assimilation time, run cm1_to_dart
    if MODEL_TO_DART:
        print "Checking for cm1_to_dart"
        if end in assim_times:
            run_cm1_to_dart(memnum, enddate)


    # MAIN PROCESSING SEQUENCE ENDS

def preclean():
    # Delete a few unneeded files
    os.system('rm dart_log.out')
    os.system('rm dart_log.nml')
    os.system('rm rsl.*')
    os.system('rm m*_run_member.py.p*')

def error_handler(error_msg, error_function):
    # Function to handle what to do with error messages
    print "!!!!!! ERROR in function: %s !!!!!!" % error_function
    print error_msg
    exit(1)


def get_assim_times():
    """ Function to look at ens_dart_param and return possible assimilation times 

        REQUIRES:
            In ens_dart_param:
                exp_length, cycle_len, N_assim, assim_start, assim_interval

        RETURNS:
                list of assimilation times"""
    # Get the cycling length
    # Get the length between cycles in minutes
    # Get the total number of assimilations
    # Find out which cycle starts the assimilation
    # and build the list of assimilation times from there
    # List indices start at 0

    # Make a list of all cycles
    # This list is in minutes since simulation start
    cycle_times = range(0, int(exp_length)+int(cycle_len), int(cycle_len))

    # Now whittle this down based on the above
    # Start where requested
    if assim_start < 0:
        cycle_times = []
    else:
        cycle_times = cycle_times[assim_start:N_assim:assim_interval]


    # Return the list of assimilation times
    return cycle_times



def get_basic_info():
    """ Function to figure out which member number we are processing based on the directory the
        script is running in.  Queue submission **MUST** include the -wd {path} command where
        {path} is the directory of the member being run.  Restart time to use is from environmental
        variable

        REQUIRES:
            Queue submission:
                qsub command **MUST** use -wd option with the member's directory specified

        RETURNS:
            memnum, start, end """

    # Try to read the environment variables
    start = int(os.environ['STARTTIME'])
    end = int(os.environ['ENDTIME'])


    # Try to find the member based on the directory we are in
    curdir = os.getcwd()
    mem = int(re.search('mems/m(\d{1,3})', curdir).groups()[0])

    return mem, start, end



def run_dart_to_cm1(mem,start):
    """ If starting post-DART, run dart_to_wrf to populate the
     wrfinput file. """

    # Refresh with new WRF_dart_param variables
    #from WRF_dart_param import *
    print "Running dart_to_wrf"

    # We know the member number, so mv the new dart state vector
    # for that member number into the working directory
    if os.path.exists('%s/wrfdart/filter_ic_new.%04d' % (dir_wrf_dom, mem)):
        os.system('cp %s/wrfdart/filter_ic_new.%04d dart_wrf_vector' % (dir_wrf_dom, mem))
    else:
        error_handler('Could not find update ic file wrfdart/filter_ic_new.%04d' % mem, \
                      'run_dart_to_wrf')

    # Write a new dart namelist for this time
    nml_good = write_dart_namelist(mem,start)
    if not nml_good:
        error_handler('Trouble writing input.nml for member %d' % mem, 'run_dart_to_wrf')

    # Loop through each domain to be sure that wrfinput files are present
    for dn in range(max_dom+1)[1:]:
        os.system('cp %s/archive_bdy/%s_wrfinput_d01 wrfinput_d01' % (dir_wrf_dom,start.strftime('%Y%m%d%H')))
        if not os.path.exists('wrfinput_d%02d' % dn):
            error_handler('Could not find wrfinput_d%02d for mem %d' % (dn,mem), \
                          'run_dart_to_wrf')

    # Link in the dart_to_wrf executable
    os.system('ln -sf %s/dart_to_wrf .' % dir_src_dart)

    # Run the executable for the current member
    if os.path.exists('dtw.out'):
        os.system('rm dtw.out')
    os.system('./dart_to_wrf >> dtw.out')

    # Check for errors in the dart_to_wrf_sequence
    darterror = False
    logfile = open('dtw.out','r')
    for line in logfile:
        if re.search('error',str(line)):
            darterror = True

    # Error out if errors found
    if darterror:
        error_handler('Error creating wrfinput files from filter_ic_new.%04d' % (mem), \
                      'run_dart_to_wrf')
    else:
        # If no errors, just clean up the directory
        os.system('rm -f dart_wrf_vector')
        os.system('rm -f dart_to_wrf')
        os.system('rm input.nml')

    # End of routine




def cm1_prep(mem,start,end,restart_name='cm1out_rst_000001.nc'):
    """ Function to write new namelist.input and check to be
        sure files are in order 

    """
    # If this is a restart, be sure we have the right file
    if start != 0:
        # Try loading the restart file with the correct name
        if not os.path.exists(restart_name):
            error_handler('unable to find file {:s} in member {:d} directory'.format(restart_name, mem), 'cm1_prep')
        # Make sure that the restart file is indeed at the correct time
        rstnc = Dataset(restart_name, 'r')
        time_sec = int(rstnc.variables['time'][0])
        rstnc.close()
        if time_sec != start*60:
            error_handler('restart file {:s} for member {:d} does not match current cycle time: {:d} min.'.format(restart_name, mem, start), 'cm1_prep')
    
    # Get the files needed to set up namelist
    if os.path.exists('namelist.input'):
        os.system('rm namelist.input')
    if not os.path.exists('write_cm1_namelist.py'):
       os.system('ln -sf {:s}/write_cm1_namelist.py .'.format(dir_dom))
    if not os.path.exists('ens_dart_param.py'):
        os.system('ln -sf {:s}/ens_dart_param.py'.format(dir_dom))

    # Find out how long the run should be 
    dtime = int(end-start)

    # Write the namelist
    if start != 0:
        os.system('./write_cm1_namelist.py -r {:s} -l {:d}'.format(restart_name, dtime))
    else:
        os.system('./write_cm1_namelist.py -l {:d}'.format(dtime))
    if not os.path.exists('namelist.input'):
        error_handler('Unable to find namelist.input','wrf_prep')

    # End of wrf_prep



def run_cm1(mem):
    """ Function that actually runs CM1
        REQURES:
            From ens_dart_param:
                mpi_run_command, mpi_numprocs_flag, dir_src_model """

    # All we do here is run CM1
    #mp_numprocs_member = 16 
    if not os.path.exists('cm1.exe'):
        os.system('ln -sf {:s}/cm1.exe'.format(dir_src_model))
    #os.system('mpirun -np %d wrf.exe' % mp_numprocs_member)
    os.system('%s %s ./cm1.exe' % (mpi_run_command,mpi_numprocs_flag))



def post_model_cleanup(mem,start,end,fcst_end):
    """ Function to move output files to appropriate places, handle
    calculating tendency if desired and process auxilliary outputs """
    print "Beginning post-model cleanup"

    # First, verify that CM1 finished correctly
    if start != 0:
        # Find the second restart name
        rst_files = [f for f in os.listdir('.') if\
                     f.startswith('cm1out_rst') and f.endswith('.nc')]
        rst_files.sort()
        # It's the seocnd file
        restart_name = rst_files[1]
        if not os.path.exists(restart_name):
            error_handler('Restart file {:s} not found for member {:d}'.format(restart_name, mem),'post_model_cleanup')
        # Check to be sure this is at the right time
        with Dataset(restart_name,'r') as rstfile:
            rst_time = rstfile.variables['time'][0]
            if rst_time != end * 60:
                error_handler('Restart file {:s} not at correct time {:d} min. for member {:d}'.format(restart_name, int(end), mem),'post_model_cleanup')

    # Also check for out file
    if not os.path.exists('cm1out.nc'):
        error_handler('Output file cm1out.nc not found for member {:d}'.format(mem), 'post_model_cleanup')


    print "CM1 restart file found.  Success!"    
    # Check the rsl.error files to be doubly sure
    # INSERT CODE HERE

    # Now run through a litany of possible post-processing things to be done

    # Check if we're computing tendency
    if flag_compute_tendency:
        print "Computing altimeter tendency"
        if not os.path.exists('wrf_tendency'):
            os.system('ln -sf %s/wrf_tendency .' % dir_utils)
        # Write a new dart namelist for this time
        nml_good = write_dart_namelist(mem,end)
        if not nml_good:
            error_handler('Trouble writing input.nml for member %d' % mem, 'post_model_cleanup')

        for dn in range(int(max_dom)+1)[1:]:
            # Don't overwrite the wrfinput_d01 file
            if dn > 1:
                os.system('cp wrfinput_d01 wrfinput_d01_orig')
                os.system('rm wrfinput_d01')
                os.system('ln -sf wrfinput_d%02d wrfinput_d01' % dn)
            os.system('ln -sf wrfout_d%02d_%s wrfout_d01' % (dn, end.strftime('%Y-%m-%d_%H:%M:%S')))
            os.system('./wrf_tendency')
            os.system('rm wrfout_d01 wrf_tendency')
            if dn > 1:
                os.system('rm wrfinput_d01')
                os.system('mv wrfinput_d01_orig wrfinput_d01')


    # Actual file turnaround goes here
    print "Copying over files for restart."
    # Save the previous file in case something goes wrong
    if start != 0:
        os.system('mv cm1out_rst_000001.nc prev_cm1out_rst_000001.nc')
        # THIS IS WHERE THE ACTUAL SWITCHOVER GOES
        print("Copying over restart:", restart_name)
        os.system('cp {:s} cm1out_rst_000001.nc'.format(restart_name))

    # Now remove all other restart files
    rst_file = [f for f in os.listdir('.') if 'cm1out_rst' in f and not f.endswith('000001.nc')]
    for f in rst_file:
        os.system('rm -f {:s}'.format(f))


    # Check to see if we keep the regular out files
    if not flag_keep_outs:
        print "Removing old out files."
        os.system('rm -f cm1out.nc')

    # Archive the forecast if fcst_len != 0
    if fcst_len:
        os.system('mv cm1out.nc cm1out_m{:d}_{:06d}.nc'.format(mem,start))



    # End of post_wrf_cleanup


def run_cm1_to_dart(mem,end):
    """ Run wrf_to_dart to prepare for assimilation"""
    print "Beginning wrf_to_dart."
    # Check to be sure that the date in wrfinput_d01 matches the end time
    for dn in range(max_dom+1)[1:]:
        # Use ncdump to get the time from the wrfinput_d?? file
        ncout = os.popen('ncdump -v Times wrfinput_d%02d' % dn)
        timeline = ncout.readlines()[-2]
        wrfin_timestring = re.search('"(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})"',timeline).groups()[0]
        wrfin_time = datetime.strptime(wrfin_timestring,'%Y-%m-%d_%H:%M:%S')

        # Make sure the wrfinput file is at the right time
        if wrfin_time != end:
            error_handler('wrfinput_d%02d time (%s) does not match assim time (%s) for mem %d' \
                      % (dn,wrfin_time.strftime('%Y%m%d%H'), end.strftime('%Y%m%d%H'), mem), \
                      'run_wrf_to_dart')   

    # Write a new dart namelist for this time
    nml_good = write_dart_namelist(mem,end)
    if not nml_good:
        error_handler('Trouble writing input.nml for member %d' % mem, 'run_wrf_to_dart')

    # Link in the dart_to_wrf executable
    os.system('ln -sf %s/wrf_to_dart .' % dir_src_dart)

    # Run the executable for the current member
    if os.path.exists('wtd.out'):
        os.system('rm wtd.out')
    os.system('./wrf_to_dart >> wtd.out')

    # Check for errors in the wrf_to_dart sequence
    darterror = False
    logfile = open('wtd.out','r')
    for line in logfile:
        if re.search('error',str(line)):
            darterror = True

    # Error out if errors found
    if darterror or not os.path.exists('dart_wrf_vector'):
        error_handler('Error creating dart_wrf_vector from wrfinput on mem %d' % (mem), \
                      'run_wrf_to_dart')
    else:
        # If no errors, just clean up the directory
        os.system('mv -f dart_wrf_vector %s/wrfdart/filter_ic_old.%04d' % (dir_wrf_dom,mem))
        os.system('rm -f wrf_to_dart')
        os.system('rm input.nml')

    # End of function run_wrf_to_dart


def write_dart_namelist(mem,dtime):
    """ Function that makes a new input.nml file in the current working directory 

        RETURNS:
            True if successful"""
 
    # Remove the old file if it exists
    if os.path.exists('input.nml'):
        os.system('rm input.nml')


    # Link in the WRF_dart_param file
    if not os.path.exists('ens_dart_param.py'):
        os.system('ln -sf {:s}/ens_dart_param.py'.format(dir_dom))
    if not os.path.exists('ens_dart_obtypes.py'):
        os.system('ln -sf {:s}/ens_dart_obtypes.py'.format(dir_dom))
    
    # Link in make_namelist_dart
    if not os.path.exists('make_namelist_dart.py'):
        os.system('ln -sf {:s}/make_namelist_dart.py .'.format(dir_dom))

    # Run make_namelist_dart
    os.system('./make_namelist_dart.py -d %s' % dtime.strftime('%Y%m%d%H'))

    # Check to be sure it was successful
    if not os.path.exists('input.nml'):
        return False
    else:
        return True


if __name__ == '__main__':
    main()

