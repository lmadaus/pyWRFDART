#!/usr/bin/env python

import os, sys
from datetime import datetime, timedelta
import re
sys.path.append('/home/disk/pvort/nobackup/lmadaus/WRF/DOMAINS/ens_july27')
from WRF_dart_param import *


"""
# New master script that will run all member-specific tasks on the queue
# Leaving the run of "filter" as the only remaining task

# Currently, this script is set up to do cycling without forecast runs 


# At present, supplying the member number and the start time comes from
# a dumped text file supplied by the master controls script


# Changing these variables will change which segments of the sequence
# to process.  This lets you fix problems with just one part of the sequence
# without having to run the whole thing.  Set all to True to run everything.
# Note that even if set to true, dart_to_wrf and wrf_to_dart won't necessarily
# run if the starttime or endtime of the cycle is not in the list of
# assimilation times
"""

PRE_CLEAN   = True
DART_TO_WRF = False 
PROCESS_BCS = True
RUN_WRF     = True
POST_WRF    = False
WRF_TO_DART = True



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
    memnum, startdate, enddate = get_basic_info() 
    print "STARTDATE:", startdate

    # Check here if the flag_direct_netcdf_io is True.  If so, set the
    # DART_TO_WRF and WRF_TO_DART sequences to False
    if flag_direct_netcdf_io:
        DART_TO_WRF = False
        WRF_TO_DART = False

    # Check to see if we need to run dart_to_wrf
    if DART_TO_WRF:
        print "Checking for dart_to_wrf"
        if startdate in assim_times:
            run_dart_to_wrf(memnum, startdate)

    # Check to see if we need to perturb boundary conditions
    # And do a final check of things
    if flag_pert_bcs:
        if PROCESS_BCS:
            print "Perturbing boundary conditions"
            perturb_bcs(memnum, startdate, enddate)
        if RUN_WRF:
            print "Prepping for WRF"
            wrf_prep(memnum, startdate, enddate)
            print "Updating boundary conditions"
            update_bcs(memnum,startdate)
    else:
        if PROCESS_BCS:
            print "Prepping for WRF"
            os.system('cp %s/archive_bdy/%s_wrfbdy_d01 wrfbdy_d01' % (dir_wrf_dom,startdate.strftime('%Y%m%d%H')))
            wrf_prep(memnum, startdate, enddate)
            print "Updating boundary conditions"
            update_bcs(memnum,startdate)
   
    # Run WRF
    if RUN_WRF:
        print "Running WRF"
        run_wrf(memnum)

    # Do post-WRF cleanup
    if POST_WRF:
        print "Doing post-WRF cleanup"
        post_wrf_cleanup(memnum,startdate,enddate)

    # If we end at an assimilation time, run wrf_to_dart
    if WRF_TO_DART:
        print "Checking for wrf_to_dart"
        if enddate in assim_times:
            run_wrf_to_dart(memnum, enddate)


    # MAIN PROCESSING SEQUENCE ENDS

def preclean():
    # Delete a few unneeded files
    os.system('rm dart_log.out')
    os.system('rm dart_log.nml')
    os.system('rm rsl.*')
    os.system('rm m*_run_member.py.p*')
    os.system('rm fg buddy_check cost_fn grad_fn jo')
    os.system('rm dart_log.*')
    os.system('rm met_em.*')
    os.system('rm gts_omb_*')
    os.system('rm rej_obs_co*')
    os.system('rm unpert_obs*')

def error_handler(error_msg, error_function):
    # Function to handle what to do with error messages
    print "!!!!!! ERROR in function: %s !!!!!!" % error_function
    print error_msg
    exit(1)


def get_assim_times():
    """ Function to look at WRF_dart_param and return possible assimilation times 

        REQUIRES:
            In WRF_dart_param:
                date_start, date_end, fct_len, N_assim, assim_start

        RETURNS:
                list of assimilation times"""
    # Get the cycling start and end times from WRF_dart_param
    try:
        cycling_start = datetime.strptime(date_start, '%Y%m%d%H')
        cycling_end = datetime.strptime(date_end, '%Y%m%d%H')
    except:
        error_handler("Could not process date_start or date_end from WRF_dart_param", 'get_assim_times')


    # Get the length between cycles in minutes
    cycling_interval = fct_len

    # Get the total number of assimilations
    total_assim_cycles = N_assim

    # Find out which cycle starts the assimilation
    # and build the list of assimilation times from there
    first_assimilation = assim_start
    # List indices start at 0, so a first_assimilation of '1' should actually be '0'
    first_assimilation = first_assimilation - 1

    # Make a list of all cycles
    cycle_times = []
    curdate = cycling_start
    while curdate <= cycling_end:
        cycle_times.append(curdate)
        curdate = curdate + timedelta(minutes=cycling_interval)

    # Now simply parse that list down using the values from WRF_dart_param
    if first_assimilation > len(cycle_times):
        assim_times = []
    else:
        assim_times = cycle_times[first_assimilation:(first_assimilation+total_assim_cycles)]

    # Return the list of assimilation times
    return assim_times



def get_basic_info():
    """ Function to figure out which member number we are processing based on the directory the
        script is running in.  Queue submission **MUST** include the -wd {path} command where
        {path} is the directory of the member being run.  Start and end dates are pulled from
        environmental variables

        REQUIRES:
            Environment variables:
                STARTDATE,ENDDATE
            Queue submission:
                qsub command **MUST** use -wd option with the member's directory specified

        RETURNS:
            memnum, startdate, enddate (in Python form)"""

    # Try to read the environment variables
    start = int(os.environ['STARTDATE'])
    end = int(os.environ['ENDDATE'])

    startdt = datetime.strptime(str(start),'%Y%m%d%H')
    enddt = datetime.strptime(str(end),'%Y%m%d%H')
    

    # Try to find the member based on the directory we are in
    curdir = os.getcwd()
    mem = int(re.search('mems/m(\d{1,3})', curdir).groups()[0])

    return mem, startdt, enddt



def run_dart_to_wrf(mem,start):
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


def perturb_bcs(mem,start,end):
    """ Function to run the BC perturbation sequence.  Uses fixed-covariance
        perturbations through the WRF-VAR framework.  Requires a met_em.d01 file
        to be present in dir_wrf_dom valid at the end time of the current cycle
        for the perturbation to work.  The actual re-compuation of the boundary
        conditions uses the program pert_wrf_bcs from DART"""
  
    print "Starting to perturb boundary conditions"

    # Go ahead and format times here for simplicity
    datem = start.strftime('%Y%m%d%H')
    datef = end.strftime('%Y%m%d%H')

    # Rename the current wrfinput file to save for later
    os.system('cp wrfinput_d01 wrfinput_d01_%s' % datem)
    os.system('mv wrfinput_d01 wrfinput_d01_orig')

    # Copy in the met_em file for the right end time
    #if not os.path.exists('%s/met_em.d01.%s.nc' % (dir_wrf_dom,end.strftime('%Y-%m-%d_%H:%M:%S'))):
    #    error_handler('Unable to find met_em.d01 file for end time %s' % datef, \
    #                  'perturb_bcs') 
    #os.system('cp %s/met_em.d01.%s.nc .' % (dir_wrf_dom, end.strftime('%Y-%m-%d_%H:%M:%S')))

    # Write the namelist and use real.exe to convert met_em file to wrfinput file
    if not os.path.exists('write_namelists.py'):
       os.system('ln -sf %s/write_namelists.py .' % dir_wrf_dom)
    os.system('./write_namelists.py -d %s -l 0 -p -f' % datef)
    #print "Running real.exe for BC end time"
    #os.system('ln -sf %s/real.exe real.exe' % dir_src_wrf)
    #os.system('%s %s ./real.exe' % (mpi_run_command, mpi_numprocs_flag))
    #os.system('rm real.exe')
    # Check to be sure this worked
    #if not os.path.exists('wrfinput_d01'):
    #    error_handler('Program real.exe failed to produce wrfinput for BC pert.',\
    #                  'perturb_bcs')

    # Perturb this new wrfinput using wrfda
    #os.system('mv wrfinput_d01 fg')
    # ADDITIONS HERE FOR PRE-EXISTING WRFINPUT and WRFBDY
    os.system('cp %s/archive_bdy/%s_wrfinput_d01 fg' % (dir_wrf_dom,datef))
    os.system('cp %s/archive_bdy/%s_wrfbdy_d01 wrfbdy_d01' % (dir_wrf_dom,datem))

    bcint_hours = int(dlbc/60)
    #os.system('./write_namelists.py -d %s -e %d -l %d -f -p' % (datef,mem,bcint_hours))
    os.system('./write_namelists.py -d %s -e %d -l %d -f -p' % (datef,mem,0))
    os.system('ln -sf %s/build/da_wrfvar.exe da_wrfvar.exe' % dir_src_wrfvar)
    os.system('ln -sf %s/run/be.dat.cv3 be.dat' % dir_src_wrfvar)
    print "Running da_wrfvar.exe"
    os.system('%s %s ./da_wrfvar.exe' % (mpi_run_command, mpi_numprocs_flag))
    if not os.path.exists('wrfvar_output'):
        error_handler('Could not find wrfvar_output file.  da_wrfvar.exe failed.',
                      'perturb_bcs')
    os.system('mv wrfvar_output wrfinput_d01_%s' % datef)
    #os.system('ncdiff wrfinput_d01_%s fg diff.nc' % datef)
    os.system('rm da_wrfvar.exe be.dat') 

    # Copy in input.nml and rename files to run pert_wrf_bc
    os.system('cp %s/wrfdart/input.nml .' % dir_wrf_dom)
    os.system('mv wrfbdy_d01 wrfbdy_this')
    os.system('ln -sf wrfinput_d01_%s wrfinput_this' % datem)
    os.system('ln -sf wrfinput_d01_%s wrfinput_next' % datef)

    # Run pert_wrf_bc
    print "Running pert_wrf_bc"
    os.system('ln -sf %s/pert_wrf_bc pert_wrf_bc' % dir_src_dart)
    os.system('./pert_wrf_bc')
    os.system('rm pert_wrf_bc')

    # Now cleanup
    os.system('./write_namelists.py -d %s -l %d -p' % (datem,fct_len))
    os.system('mv wrfbdy_this wrfbdy_d01')
    os.system('mv wrfinput_d01_orig wrfinput_d01')
 
    os.system('rm wrfinput_d01_%s wrfinput_d01_%s' % (datem, datef))
    os.system('rm wrfinput_this wrfinput_next')
    os.system('rm fg buddy_check cost_fn grad_fn jo')
    os.system('rm dart_log.*')
    os.system('rm met_em.*')
    os.system('rm rsl.*')
    os.system('rm gts_omb_*')
    os.system('rm rej_obs_co*')
    os.system('rm unpert_obs*')

    # END of perurb_bcs sequence




def update_bcs(mem,start):
    """ Function to run update_wrf_bcs to match the initial BC field for this time
        to the updated initial conditions for this time"""
    os.system('ln -sf %s/update_wrf_bc .' % dir_src_dart)
    nml_good = write_dart_namelist(mem,start)
    if not nml_good:
        error_handler('Trouble writing input.nml for member %d' % mem, 'update_bcs')

    # Since we already passed wrf_prep, we know the namelist.input is of the right time.
    # All we have to do is run the utility.
    print "Updating WRF BCS"
    os.system('./update_wrf_bc >> uwb.out')

    # Check for errors in the update_wrf_bc sequence
    darterror = False
    logfile = open('uwb.out','r')
    for line in logfile:
        if re.search('error',str(line)):
            darterror = True

    # Error out if errors found
    if darterror:
        error_handler('Error running update_wrf_bc for member %d' % (mem), \
                      'update_bcs')
    else:
        # If no errors, clean up
        os.system('rm -f update_wrf_bc')
        #os.system('rm -f uwb.out')
    # End of routine


def wrf_prep(mem,start,end):
    """ Function to write new namelist.input and check to be
        sure files are in order 

    """

    for dn in range(max_dom+1)[1:]:
        # Use ncdump to get the time from the wrfinput_d?? file
        ncout = os.popen('ncdump -v Times wrfinput_d%02d' % dn)
        timeline = ncout.readlines()[-2]
        wrfin_timestring = re.search('"(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})"',timeline).groups()[0]
        wrfin_time = datetime.strptime(wrfin_timestring,'%Y-%m-%d_%H:%M:%S')
        #print "timeline:", timeline
        #print wrfin_time
        #print start
        # Make sure the wrfinput file is at the right time
        if wrfin_time != start:
            error_handler('wrfinput_d%02d time (%s) does not match start time (%s) for mem %d' \
                      % (dn,wrfin_time.strftime('%Y%m%d%H'), start.strftime('%Y%m%d%H'), mem), 'wrf_prep')

    # Get the files needed to set up namelists
    if os.path.exists('namelist.input'):
        os.system('rm namelist.input')
    if not os.path.exists('write_namelists.py'):
       os.system('ln -sf %s/write_namelists.py .' % dir_wrf_dom)
    if not os.path.exists('WRF_dart_param.py'):
        os.system('ln -sf %s/wrfdart/WRF_dart_param.py' % (dir_wrf_dom))

    # Find out how long the run should be 
    dtime = end - start
    dtime_hours = int(dtime.days * 24 + dtime.seconds / 3600)



    # Write the namelist
    os.system('./write_namelists.py -d %s -l %d -p -e %d' % (start.strftime('%Y%m%d%H'),dtime_hours, mem))
    if not os.path.exists('namelist.input'):
        error_handler('Unable to find namelist.input','wrf_prep')

    # End of wrf_prep



def run_wrf(mem):
    """ Function that actually runs WRF
        REQURES:
            From WRF_dart_param:
                mp_numprocs_member """

    # All we do here is run WRF
    mp_numprocs_member = 16 
    if not os.path.exists('wrf.exe'):
        os.system('ln -sf %s/wrf.exe' % WRFRUNDIR)
    #os.system('mpirun -np %d wrf.exe' % mp_numprocs_member)
    os.system('%s %s ./wrf.exe' % (mpi_run_command,mpi_numprocs_flag))



def post_wrf_cleanup(mem,start,end):
    """ Function to move wrfout files to appropriate places, handle
    calculating tendency if desired and process auxilliary outputs """
    print "Beginning post-wrf cleanup"

    # First, verify that WRF finished correctly
    for dn in range(int(max_dom)+1)[1:]:
        if not os.path.exists('wrfout_d%02d_%s' % (dn,end.strftime('%Y-%m-%d_%H:%M:%S'))):
            error_handler('File wrfout_d%02d_%s not found for member %d' \
                         % (dn,end.strftime('%Y-%m-%d_%H:%M:%S'),mem), 'post_wrf_cleanup')

    print "Wrfout file found at endtime.  Success!"    
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
            error_handler('Trouble writing input.nml for member %d' % mem, 'post_wrf_cleanup')

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


    # Now check if there are any tslists
    for file in os.listdir('.'):
        if file.endswith('.TS'):
            os.system('mv %s %s_%s' % (file,file,start.strftime('%Y%m%d%H')))
    
    # Look for Ryan Torn's pressure tendency output and process it
    for dn in range(int(max_dom)+1)[1:]:
        if os.path.exists('psfc_data_d%02d.nc' % (dn)):
           print "Processing psfc_data for domain", dn
           if not os.path.exists('input.nml'):
               nml_good = write_dart_namelist(mem,end)
               if not nml_good:
                   error_handler('Trouble writing input.nml for member %d' % mem, 'post_wrf_cleanup')
           os.system('mv psfc_data_d%02d.nc wrfout.nc' % dn)
           os.system('%s/domain_psfc_tend' % dir_utils)
           os.system('mv psfc_stats.nc psfc_stat_d%02d_%s_m%d.nc' % \
               (dn,start.strftime('%Y%m%d%H'),mem))
           os.system('rm wrfout.nc')

    # Actual file turnaround goes here
    print "Copying over wrfout to wrfinput for restart."
    for dn in range(int(max_dom)+1)[1:]:
        # Save the previous wrfinput file in case something goes wrong
        os.system('cp wrfinput_d%02d prev_wrfinput_d%02d' % (dn,dn))
        os.system('cp wrfout_d%02d_%s wrfinput_d%02d' % (dn,end.strftime('%Y-%m-%d_%H:%M:%S'),dn))   


    # Check to see if we keep the wrfout files
    wrfout_fmt = end.strftime('%Y-%m-%d_%H:%M:%S')
    if not flag_keep_wrfouts:
        print "Removing old wrfout files."
        filelist = os.listdir('%s/m%d' % (dir_members,mem))
        # Get wrfout files
        wrfouts = [f for f in filelist if f.startswith('wrfout_d')]
        # Don't delete the most recent wrfouts
        for dom in range(int(max_dom)+1)[1:]:
            wrfouts.remove('wrfout_d0%(dom)d_%(wrfout_fmt)s' % locals())
            wrfouts.remove('wrfout_d%02d_%s' % (dom,start.strftime('%Y-%m-%d_%H:%M:%S')))
        # Now remove what's left
        #print "To remove...:"
        #print wrfouts
        #raw_input()
        for file in wrfouts:
            os.system('rm %s/m%d/%s' % (dir_members,mem,file))
        # For now...
        #print "Remaining contents of directory:"
        #dirlist = os.listdir('%s/m%d' % (dir_members,mem))
        #print dirlist
        #raw_input()




    # End of post_wrf_cleanup


def run_wrf_to_dart(mem,end):
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
    if not os.path.exists('WRF_dart_param.py'):
        os.system('ln -sf %s/wrfdart/WRF_dart_param.py' % (dir_wrf_dom))
    if not os.path.exists('WRF_dart_obtypes.py'):
        os.system('ln -sf %s/wrfdart/WRF_dart_obtypes.py' % (dir_wrf_dom))
    
    # Link in make_namelist_dart
    if not os.path.exists('make_namelist_dart.py'):
        os.system('ln -sf %s/wrfdart/make_namelist_dart.py .' % (dir_wrf_dom))

    # Run make_namelist_dart
    os.system('./make_namelist_dart.py -d %s' % dtime.strftime('%Y%m%d%H'))

    # Check to be sure it was successful
    if not os.path.exists('input.nml'):
        return False
    else:
        return True


if __name__ == '__main__':
    main()

