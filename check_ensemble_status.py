#!/usr/bin/env python

# Script that checks the status of the entire model run
from __future__ import print_function, division
import os, sys, datetime, getopt, re
from netCDF4 import Dataset
from ens_dart_param import *
resub = False 
silent = False 

(opts, args) = getopt.getopt(sys.argv[1:],'d:s')
# Command line options to input the
# date we are checking and to run silently
# (for unsupervised runs)
for o,a in opts:
   if o == '-d':
      indate = a
   if o == '-s':
      silent = True

# Convert indate to seconds
indate = int(indate)
# Stat with several functions


def chkstat(chkdate):
    import os, sys, datetime, re
    # This is equivalent to sourcing the
    # WRF_dart_param file for just the
    # variables we need
    from ens_dart_param import Ne, dir_members, dir_dom

    # Convert the date we are checking to a datetime
    # and then re-write it in WRF file format   
    #chkdatedt = datetime.datetime.strptime(chkdate, '%Y%m%d%H')
    #wrftime = chkdatedt.strftime('%Y-%m-%d_%H:00:00')
    chktime = chkdate 

    # Loop through the ensemble directories and see if any
    # rst files for the given input time are missing

    # Four categories possible
    memsnotstarted = []
    memsdone = []
    memsnotdone = []
    memserror = []
   
    # Parse the queue output to see what's running
    if os.uname()[0] == 'AIX' or os.uname()[1].startswith('yslog'):
       qstat_out = os.popen('bjobs').readlines()
       #print qstat_out
       #raw_input()
       queue_members = {} 
       for line in qstat_out:
           if re.search('member_(\d{1,3})',line):
               linesp = line.split()
               #print line,linesp[6]
               # Job name could be in position 5 or 6
               try:
                   memnum = int(re.search('member_(\d{1,3})',linesp[6]).groups()[0])
               except:
                   memnum = int(re.search('member_(\d{1,3})',linesp[5]).groups()[0])
               runstat = linesp[2]
               # Convert the job status to Sun Grid Manager format
               if runstat.strip() == 'RUN':
                   runstat = 'r'
               elif runstat.strip() == 'PEND':
                   runstat = 'qw'
               jobid = int(linesp[0])
               queue_members[memnum] = (runstat, jobid)

    else:
        # On Local Queue
        qstat_out = os.popen('qstat').readlines()
        #print qstat_out
        queue_members = {} 
        for line in qstat_out:
            if re.search('m(\d{1,3})_run_',line):
                linesp = line.split()
                memnum = int(re.search('m(\d{1,3})_run_',linesp[2]).groups()[0])
                runstat = linesp[4]
                jobid = int(linesp[0])
                queue_members[memnum] = (runstat, jobid)


    # Loop through all ensemble members
    for mem in range(1,Ne+1):
        # Loop through all domainis
        #for dom in range(max_dom+1)[1:]:
        # NEW LOGIC HERE FOR DETERMINING STATUS
        # First check to see if this member exists in the queue
        if mem in queue_members.keys():
            if queue_members[mem][0].strip() == 'r':
                memsnotdone.append(mem)
            elif queue_members[mem][0].strip().lower() in ['t','qw']:
                memsnotstarted.append(mem)
            else:
                # Any other status is probably not good
                memserror.append(mem)
        else:
            # If we're here, then the member is not in the queue.  Check for wrfout
            rstfile = '{:s}/m{:d}/cm1out_rst_000001.nc'.format(dir_members,mem)
            if not os.path.exists(rstfile) and mem not in memserror:
                print("Restart file not found")
                # Member is not in the queue, so there must be a failure.
                error_found = check_logfile(mem)
                if not error_found:
                    # Not sure why this member is missing, but add it to error anyway
                    print("Member {:d} Model crashed.  Unknown reason.".format(mem))
                    memserror.append(mem)

            # If we're here, then restart file exists. Now the true test...is it at the
            # right time?
            with Dataset(rstfile, 'r') as rstnc:
                rst_time = int(rstnc.variables['time'][0])
                # If the restart time doesnt match and we are not in the queue or running,
                # the member most have crashed
                if (rst_time != chktime):
                    # Because we're in this part of the if statement, we know
                    # the member is not running or waiting
                    print("Restart time:", rst_time, "Does not match check time:", chktime)
                    error_found = check_logfile(mem)
                    if not error_found:
                        # Not sure why this member is missing, but add it to error anyway
                        print("Member {:d} Model crashed.  Unknown reason.".format(mem))
                        memserror.append(mem)

            # Final check if we're not doing the netcdf IO
            if not os.path.exists('{:s}/assimilation/filter_ic_old.{:04d}'.format(dir_dom,mem))\
                and (mem not in memserror) and not flag_direct_netcdf_io:
                # Something went wrong on WTD.  Only can proceed if this file exists
                print("Member {:d} Error in POST-MODEL or MODEL-TO-DART".format(mem))
                memserror.append(mem)
            # If we're here, the member is no longer in the queue, and restart
            # (and filter_ic_old) files exists for the correct time.  This member is done. 
            if mem not in memserror:
                memsdone.append(mem)


    # The list(set()).sort() syntax removes duplicates if there are multiple domains
    # Otherwise, return a list of which members are in which category
    return sorted(list(set(memsdone))), sorted(list(set(memsnotdone))), sorted(list(set(memsnotstarted))), sorted(list(set(memserror)))


def check_logfile(memnum):
    error = 0
    # If we found a logfile of some kind (rsl.error or log.out)
    try:
        logfile = open('%s/m%d/rsl.error.0000' % (dir_members,memnum), 'r') 
    except:
        # Could not find a logfile, so return error = 0
        return error
    # Look at the last 20 lines to see if there is an error message
    for line in logfile.readlines()[-20:]:
        # Below are various common error messages
        # and what they mean.  If an error message
        # is found, add the model to the errored members list
        if re.search('used in new version', line) or re.search('cfl', line):
            print("Member {:d}: Probable CFL Error".format(ie))
            print("     ", line)
            error = 1
        if re.search('forrtl: error', line):
            print("Member {:d}: Model crashed, unknown reason.".format(ie))
            print("     ", line)
            error = 1
        if re.search('recursive I/O operation', line):
            print("Member {:d}: WRF I/O crashed, unknown reason.".format(ie))
            print("     ",line) 
            error = 1
        if re.search('Segmentation Fault', line) or re.search('segmentation fault', line):
            print("Member {:d}: Memory error (seg fault)".format(ie))
            print("     ", line)
            error = 1
        if re.search('WOULD GO OFF TOP', line):
            print("Member {:d}: CFL error with convective scheme".format(ie))
            print("     ", line)
            error = 1
    return error






def check_complete_cm1(indate):
    # Checks to see if any of the wrfout files
    # for time indate have NAN errors in them
    # Equivalent to sourcing ens_dart_param
    from ens_dart_param import Ne,dir_members
    try:
        from netCDF4 import Dataset
        nonan = False     
    except:
        # May not have this library on some machines
        print("netCDF4 library not found. skipping Nancheck")
        nonan = True
    try:
        import numpy as np
    except:
        print("Numpy not found.  Skipping.")


    # Loop through each ensemble member and see if the T2 field
    # is viable

    # Open the wrfout file for this time
    # Here, write the time in WRF format

    nan_members = []
    if not nonan:
       for mem in range(Ne+1)[1:]:
           # Loop through all members
           try:
               # Open the dataset if it exists
               ncfile = Dataset('{:s}/m{:d}/cm1out_rst_000001.nc'.format(dir_members,mem))
           except:
               print("rst file at time {:d} min. doesn't seem to exist for member {:d}!".format(indate, mem))
               exit(1)
           T2 = ncfile.variables['t2'][:,:]
           # Search for nan-ed members
           if True in np.isnan(T2):
               print("Found NAN in member {:d}".format(mem))
               nan_members.append(mem)
    return nan_members




def resubmit(merror):
   # Quick function to resubmit all crashed members
   for mem in merror:
      print("   Resubmitting member {:d}".format(mem))
      os.chdir('{:s}/m{:d}'.format(dir_members,mem))
      os.system('rm -rf rsl.*')
      # May need to reset environment variables here
      # to have proper start time for resubmission.  (TODO)
      #if os.uname()[1] in ['enkf']:
      #    os.system('qsub -V -pe ompi %d -q %s -wd %s/m%d m%d_run_member.py'\
      #              % (mpi_numprocs_member, queue_members, dir_member, mem, mem))
      os.chdir(dir_dom)

if silent:
    # SILENT CONTROL MODE
    import time
    # Relics from the old script--make sure we're starting from scratch
    if os.path.exists('ensemble_done_{:d}'.format(indate)):
      os.system('rm -rf ensemble_done_{:d}'.format(indate))
    if os.path.exists('master_ens_log'):
      os.system('rm -rf master_ens_log')

    # Initialize some variables and the log file
    memsdone = 0
    logfile = open('master_ens_log','w')
    resub = 0 
    timecheck = 0
    while memsdone < Ne:
        # Check every ten seconds to see if all members are done
        # if resub is True (1), then resubmit members as they crash
        time.sleep(10)
        # Master control lock -- check to see if this file exists                              
        # If it doesn't exist, exit the program                                                
        if not os.path.exists('{:s}/AUTO_RUN_IN_PROGRESS'.format(dir_dom)):                        
            print("File 'AUTO_RUN_IN_PROGRESS' is not present in main dir. Exiting.")         
            exit(0) 

        # Check the status for the date specified
        mdone, mnotdone, mnotstart, merror = chkstat(indate)
        if len(merror)>0 and resub:
            # If some members have crashed (more than zero members in merror), resubmit
            # (if flag is set)
            logfile.write("")
            logfile.write("Resubmitting crashed members:", merror)
            logfile.write("")
            # Call the resubmit function
            resubmit(merror)

        if timecheck == 180:
            # Periodically write to the log file
            nowtm = datetime.datetime.now()
            logfile.write("")
            logfile.write("***  Status as of {:%m/%d  %I:%M:%S %p}  ***".format(nowtm)) 
            logfile.write("-----------------------------------------")
            logfile.write("   {:02d} Members done: ".format(len(mdone)))
            logfile.write("   {:02d} Members in progress: ".format(len(mnotdone)))
            logfile.write("   {:02d} Members not started: ".format(len(mnotstart)))
            logfile.write("   {:02d} Members crashed: ".format(len(merror)))
            logfile.write("")
            timecheck = 0

        timecheck = timecheck + 1
        # Find how many members are done now for the next run through
        # the while loop
        memsdone = len(mdone)

    # Once the silent mode while loop exits, check again to be sure
    # all ensemble members are done
    mdone, mnotdone, mnotstart, merror = chkstat(indate)
    if len(mdone) == Ne:
        # If all files are done...
        nowtm = datetime.datetime.now()
        logfile.write("")
        logfile.write("***  Status as of {:%m/%d  %I:%M:%S %p}  ***".format(nowtm)) 
        logfile.write("******* ALL FILES PRESENT *******")
        print("******* ALL FILES PRESENT *******")
        # Now check for any NANed members
        # Using the check_complete_cm1 function
        nanmems = check_complete_cm1(indate) 
        if len(nanmems) > 0:
            #resub = raw_input("Resubmit NAN members (0 or 1)?  ")
            #if int(resub) == 1:
            # Loop through the members, copy in old restart file and resubmit
            #for mem in nanmems:
            #    os.system('cp {:s}/m{:d}/prev_cm1out_rst_000001.nc {:s}/m{:d}/cm1out_rst_000001.nc'.format(dir_members,mem,dir_member,mem))
            logfile.write("Resubmitting NANed members:", nanmems)
            resubmit(nanmems)
            # Wait for the nan-ed members to finish
            # LEM (TODO) NEED CODE HERE FOR MONITORING RESTART
            os.system('touch ensemble_done_{:d}'.format(int(indate/60)))
            logfile.close()
            exit(0)

        else:
            # If all is good, we are done here
            os.system('touch ensemble_done_{:d}'.format(int(indate/60)))
            logfile.close()
            exit(0)

    else:
        # We only go here if the while loop somehow exited but
        # not all ensemble members were actually done.
        nowtm = datetime.datetime.now()
        logfile.write("")
        logfile.write("***  Status as of {:%m/%d  %I:%M:%S %p}  ***".format(nowtm)) 
        logfile.write("Problem occurred--not all files present.  Ending.")
        os.system('touch ensemble_error_{:d}'.format(indate/60))
        exit(1)
 

else:
    # INTERACTIVE MODE
    # This just displays a text output of the current
    # Status of the ensemble 
    mdone, mnotdone, mnotstart, merror = chkstat(indate)
    nowtm = datetime.datetime.now()
    print("")
    print("***  Status as of {:%m/%d  %I:%M:%S %p}  ***".format(nowtm))
    print("-----------------------------------------")
    print("   {:02d} Members done: ".format(len(mdone)), mdone)
    print("   {:02d} Members in progress: ".format(len(mnotdone)), mnotdone)
    print("   {:02d} Members not started: ".format(len(mnotstart)), mnotstart)
    print("   {:02d} Members crashed: ".format(len(merror)), merror)
    print("")

    if len(mdone) >= int(Ne):
        # If all members happen to be done, check for NAN-ed members
        print("******* ALL FILES PRESENT *******")
        print("")
        # Now check for any NANed members
        nanmems = check_complete_cm1(indate) 
        print(nanmems)
        if len(nanmems) > 0:
            resub = raw_input("Resubmit NAN members (0 or 1)?  ")
            if int(resub) == 1:
                # Loop through the members, copy in the old restart file and resubmit
                #for mem in nanmems:
                #    os.system('cp {:s}/m{:d}/prev_cm1out_rst_000001.nc {:s}/m{:d}/cm1out_rst_000001.nc'.format(dir_members,mem,dir_member,mem))
                resubmit(nanmems)
            else:
                pass
    if len(merror) > 0:
        # If there are any crashed members, prompt to resubmit
        resub = raw_input("Resubmit crashed members (0 or 1)?  ")
    else:
        pass

    if int(resub) == 1:
        # Loop through all crashed members and resubmit
        # If requesetd
        resubmit(merror)






