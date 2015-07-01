#!/usr/bin/env python

# Script that checks the status of the entire model run

import os, sys, datetime, getopt, re
sys.path.append('./wrfdart')

from WRF_dart_param import *
resub = False 
silent = False 

(opts, args) = getopt.getopt(sys.argv[1:],'d:s')
# Command line options to input the
# date we are checking and to run silently
# (for unsupervised runs)
for o,a in opts:
   if o == '-d':
      datestr = a
   if o == '-s':
      silent = True


# Stat with several functions


def chkstat(chkdate):
   import os, sys, datetime, re
   # This is equivalent to sourcing the
   # WRF_dart_param file for just the
   # variables we need
   from WRF_dart_param import Ne, max_dom, dir_members, dir_wrf_dom

   # Convert the date we are checking to a datetime
   # and then re-write it in WRF file format   
   chkdatedt = datetime.datetime.strptime(chkdate, '%Y%m%d%H')
   wrftime = chkdatedt.strftime('%Y-%m-%d_%H:00:00')

   # Loop through the ensemble directories and see if any
   # wrfout files for the given input time are missing

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
       # On enkf
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
      for dom in range(max_dom+1)[1:]:
         # Check if the wrfout file for that time exists
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
             if not os.path.exists('%s/m%d/wrfout_d%02d_%s' % (dir_members,mem,dom,wrftime))\
                    and mem not in memserror:
                 # Member is not in the queue, so there must be a failure.
                 error_found = check_logfile(mem)
                 if not error_found:
                     # Not sure why this member is missing, but add it to error anyway
                     print "Member %d Dom %d Model crashed.  Unknown reason." % (mem,dom)
                     memserror.append(mem)      

             # If we're here, then wrfout file exists at the valid time. Check for the
             # filter_ic_old file for this member
             if not os.path.exists('%s/wrfdart/filter_ic_old.%04d' % (dir_wrf_dom,mem))\
                    and mem not in memserror:
                 # Something went wrong on WTD.  Only can proceed if this file exists
                 print "Member %d Dom %d Error in POST-WRF or WRF-TO-DART" % (mem,dom)
                 memserror.append(mem)
             # If we're here, the member is no longer in the queue, and both wrfout
             # and filter_ic_old files exists for the correct time.  This member is done. 
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
            print "Member %d Dom %d : Probable CFL Error" % (ie,dom)
            print "     %s" % line
            error = 1
        if re.search('forrtl: error', line):
            print "Member %d Dom %d: Model crashed, unknown reason." % (ie, dom)
            print "     %s" % line
            error = 1
        if re.search('recursive I/O operation', line):
            print "Member %d Dom %d: WRF I/O crashed, unknown reason." % (ie,dom)
            print "     %s" % line  
            error = 1
        if re.search('Segmentation Fault', line) or re.search('segmentation fault', line):
            print "Member %d Dom %d: Memory error (seg fault)" % (ie,dom)
            print "     %s" % line
            error = 1
        if re.search('WOULD GO OFF TOP', line):
            print "Member %d Dom %d: CFL error with convective scheme" % (ie,dom)
            print "     %s" % line
            error = 1
    return error






def check_complete_wrf(datestr):
    # Checks to see if any of the wrfout files
    # for time datestr have NAN errors in them
    import sys
    sys.path.append('./wrfdart')
    # Equivalent to sourcing WRF_dart_param
    from WRF_dart_param import Ne,dir_members
    try:
        from netCDF4 import Dataset
        nonan = False     
    except:
        # May not have this library on some machines
        print "netCDF4 library not found. skipping Nancheck"
        nonan = True
    from datetime import datetime
    try:
        import numpy as np
    except:
        print "Numpy not found.  Skipping."
    # Set the date time
    try:
        indate = datetime.strptime(datestr,'%Y%m%d%H')
    except:
        print "Unrecognized date/time: ", datestr
        exit(1)


    # Loop through each ensemble member and see if the T2 field
    # is viable

    # Open the wrfout file for this time
    # Here, write the time in WRF format
    wrftime = indate.strftime('%Y-%m-%d_%H:00:00')

    nan_members = []
    if not nonan:
       for mem in range(Ne+1)[1:]:
           # Loop through all members
           try:
               # Open the dataset if it exists
               ncfile = Dataset('%s/m%d/wrfout_d01_%s' % (dir_members,mem,wrftime))
           except:
               print "Wrfout file at time %s doesn't seem to exist for member %d!" % (wrftime, mem)
               exit(1)
           T2 = ncfile.variables['T2'][:,:]
           # Search for nan-ed members
           if True in np.isnan(T2):
               print "Found NAN in member %d" % mem
               nan_members.append(mem)
    return nan_members




def resubmit(merror):
   # Quick function to resubmit all crashed members
   for mem in merror:
      print "   Resubmitting member %d" % mem
      os.chdir('%s/m%d' % (dir_members,mem))
      os.system('rm -rf rsl.*')
      # May need to reset environment variables here
      # to have proper start time for resubmission.  (TODO)
      #if os.uname()[1] in ['enkf']:
      #    os.system('qsub -V -pe ompi %d -q %s -wd %s/m%d m%d_run_member.py'\
      #              % (mpi_numprocs_member, queue_members, dir_member, mem, mem))
      os.chdir(dir_wrf_dom)

if silent:
   # SILENT CONTROL MODE
   import time
   # Relics from the old script--make sure we're starting from scratch
   if os.path.exists('ensemble_done_%s' % datestr):
      os.system('rm -rf ensemble_done_%s' % datestr)
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
      if not os.path.exists('%s/AUTO_RUN_IN_PROGRESS' % dir_wrf_dom):                        
          print "File 'AUTO_RUN_IN_PROGRESS' is not present in main dir.  Exiting."          
          exit(0) 

      # Check the status for the date specified
      mdone, mnotdone, mnotstart, merror = chkstat(datestr)
      if len(merror)>0 and resub:
         # If some members have crashed (more than zero members in merror), resubmit
         # (if flag is set)
         print >>logfile, ""
         print >>logfile, "Resubmitting crashed members:", merror
         print >>logfile, ""
         # Call the resubmit function
         resubmit(merror)

      if timecheck == 180:
         # Periodically write to the log file
         nowtm = datetime.datetime.now()
         print >>logfile, ""
         print >>logfile, "***  Status as of %s  ***" % nowtm.strftime('%m/%d  %I:%M:%S %p') 
         print >>logfile, "-----------------------------------------"
         print >>logfile, "   %02d Members done: " % len(mdone), mdone
         print >>logfile, "   %02d Members in progress: " % len(mnotdone), mnotdone
         print >>logfile, "   %02d Members not started: " % len(mnotstart), mnotstart
         print >>logfile, "   %02d Members crashed: " % len(merror), merror
         print >>logfile, ""
         timecheck = 0

      timecheck = timecheck + 1
      # Find how many members are done now for the next run through
      # the while loop
      memsdone = len(mdone)

   # Once the silent mode while loop exits, check again to be sure
   # all ensemble members are done
   mdone, mnotdone, mnotstart, merror = chkstat(datestr)
   if len(mdone) == Ne:
      # If all files are done...
      nowtm = datetime.datetime.now()
      print >>logfile, ""
      print >>logfile, "***  Status as of %s  ***" % nowtm.strftime('%m/%d  %I:%M:%S %p') 
      print >>logfile, "******* ALL FILES PRESENT *******"
      print "******* ALL FILES PRESENT *******"
      # Now check for any NANed members
      # Using the check_complete_wrf function
      nanmems = check_complete_wrf(datestr) 
      if len(nanmems) > 0:
          #resub = raw_input("Resubmit NAN members (0 or 1)?  ")
          #if int(resub) == 1:
          # Loop through the members, copy in a new wrfbdy file and resubmit
          for mem in nanmems:
              os.system('cp %s/m%d/wrfbdy_d01 %s/m%d/wrfbdy_d01' % (dir_members,(mem-1),dir_members,mem))
              datedt = datetime.datetime.strptime(datestr,'%Y%m%d%H')
              os.system('rm %s/m%d/wrfout_d01_%s' % (dir_members,mem,datedt.strftime('%Y-%m-%d_%H:00:00')))
          print >>logfile, "Resubmitting NANed members:", nanmems
          resubmit(nanmems)
          # Wait for the nan-ed members to finish
          for mem in nanmems:
              while not os.path.exists('%s/m%d/wrfout_d01_%s' % (dir_members,mem,datedt.strftime('%Y-%m-%d_%H:00:00'))):
                  time.sleep(10)
          os.system('touch ensemble_done_%s' % datestr)
          logfile.close()
          exit(0)

      else:
         # If all is good, we are done here
         os.system('touch ensemble_done_%s' % datestr)
         logfile.close()
         exit(0)

   else:
      # We only go here if the while loop somehow exited but
      # not all ensemble members were actually done.
      nowtm = datetime.datetime.now()
      print >>logfile, ""
      print >>logfile, "***  Status as of %s  ***" % nowtm.strftime('%m/%d  %I:%M:%S %p') 
      print >>logfile, "Problem occurred--not all files present.  Ending."
      os.system('touch ensemble_error_%s' % datestr)
      exit(1)
 

else:
   # INTERACTIVE MODE
   # This just displays a text output of the current
   # Status of the ensemble 
   mdone, mnotdone, mnotstart, merror = chkstat(datestr)
   nowtm = datetime.datetime.now()
   print ""
   print "***  Status as of %s  ***" % nowtm.strftime('%m/%d  %I:%M:%S %p') 
   print "-----------------------------------------"
   print "   %02d Members done: " % len(mdone), mdone
   print "   %02d Members in progress: " % len(mnotdone), mnotdone
   print "   %02d Members not started: " % len(mnotstart), mnotstart
   print "   %02d Members crashed: " % len(merror), merror
   print ""

   if len(mdone) >= int(Ne):
      # If all members happen to be done, check for NAN-ed members
      print "******* ALL FILES PRESENT *******"
      print ""
      # Now check for any NANed members
      nanmems = check_complete_wrf(datestr) 
      print nanmems
      if len(nanmems) > 0:
          resub = raw_input("Resubmit NAN members (0 or 1)?  ")
          if int(resub) == 1:
             # Loop through the members, copy in a new wrfbdy file and resubmit
             for mem in nanmems:
                 os.system('cp %s/m%d/wrfbdy_d01 %s/m%d/wrfbdy_d01' % (dir_members,(mem-1),dir_members,mem))
                 datedt = datetime.datetime.strptime(datestr,'%Y%m%d%H')
                 os.system('rm %s/m%d/wrfout_d01_%s' % (dir_members,mem,datedt.strftime('%Y-%m-%d_%H:00:00')))
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






