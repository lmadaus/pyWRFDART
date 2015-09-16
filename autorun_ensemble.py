#!/usr/bin/env python

import os, sys, getopt
from ens_dart_param import *
from datetime import datetime,timedelta
import time

os.chdir(dir_dom)

#gfs_dir = '/glade/p/work/lmadaus/rtfiles'

def main():

    # Make a list of all model times
    starttime = 0
    endtime = exp_length
    deltatime = int(cycle_len)

    opts,args =getopt.getopt(sys.argv[1:],'d:')
    for o,a in opts:
        if o == '-d':
            curtime = int(a)



    all_timelist = []
    run_timelist = []
    while curtime < endtime:
        curtime += deltatime
        run_timelist.append(curtime)

    curtime = starttime
    while curtime <= endtime:
        all_timelist.append(curtime)
        curtime += deltatime

    #print timelist

    # Find out which are assimilation times
    if assim_start < 0:
        assim_times = []
    else:
        assim_times = all_timelist[assim_start:N_assim+1:assim_interval]

    #print "Assim times correct?", assim_times
    #raw_input()
    #print "Run_datelist", run_datelist

    # Turn on the control lock
    if not os.path.exists('%s/AUTO_RUN_IN_PROGRESS' % dir_wrf_dom):                           
        os.system('touch %s/AUTO_RUN_IN_PROGRESS' % dir_wrf_dom)                                            
    else:                                                                                     
        print "Found file 'AUTO_RUN_IN_PROGRESS' in main dir."                                
        print "Don't want to start a duplicate auto run."                                     
        print "Exiting."                                                                      
        exit(1)             
    for time in run_timelist:

        print("Auto-running in silent mode building to time", time, "min")
        os.system('./check_ensemble_status.py -d {:d} -s'.format(time))
        # Master control lock -- check to see if this file exists                              
        # If it doesn't exist, exit the program                                                
        if not os.path.exists('{:s}/AUTO_RUN_IN_PROGRESS'.format(dir_wrf_dom)):                        
            print("File 'AUTO_RUN_IN_PROGRESS' is not present in main dir.  Exiting.")      
            exit(0)         
        if time in assim_times:
            #os.chdir('%s/wrfdart' % dir_wrf_dom)
            os.system('./submit_filter.py -d {:d} -m {:d}'.format(time,mpi_numprocs_filter))
            #os.chdir(dir_wrf_dom)

        # Run WPS to get new BCs
        #run_wps_sequence(date,deltatime)

        # Submit the new ensemble
        os.system('./submit_all.py -d {:d}'.format(time))
if __name__ == '__main__':
    main()
