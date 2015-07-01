#!/usr/bin/env python

import os, sys, getopt
from WRF_dart_param import *
from datetime import datetime,timedelta
import time

os.chdir(dir_wrf_dom)

gfs_dir = '/glade/p/work/lmadaus/rtfiles'

def run_wps_sequence(ldate,fcst_len):
    # Clean up the old
    indate = datetime.strptime(ldate,'%Y%m%d%H')
    os.chdir(dir_wrf_dom)
    os.system('rm -f FILE:*')
    os.system('rm -f met_em*')
    os.system('rm -f GRIBFILE*')
    os.system('rm -f wrfbdy_d01')
    os.system('rm -f gfs*')

    os.system('./write_namelists.py -d %s -l %d -f -p' % (indate.strftime('%Y%m%d%H'),int(fcst_len.seconds/3600)))
    # Find the correct GFS forecasts
    if indate.hour in [0,6,12,18]:
        gfs_date = indate - timedelta(hours=6)
        os.system('cp %s/%s/gfs_f06.grb .' % (gfs_dir,gfs_date.strftime('%Y%m%d%H')))
        os.system('cp %s/%s/gfs_f09.grb .' % (gfs_dir,gfs_date.strftime('%Y%m%d%H')))
    else:
        gfs_date = indate - timedelta(hours=3)
        os.system('cp %s/%s/gfs_f03.grb .' % (gfs_dir,gfs_date.strftime('%Y%m%d%H')))
        os.system('cp %s/%s/gfs_f06.grb .' % (gfs_dir,gfs_date.strftime('%Y%m%d%H')))

    # Run linkgrib
    os.system('./link_grib.csh gfs_')
    print "Ungribbing..."
    # Run ungrib
    os.system('./ungrib.exe')
    # Run metgrid
    os.system('./metgrid.exe')
    # Run real
    os.system('./run_real.py')

    # Wait for wrfbdy to exist
    print "Waiting for real.exe to finish."
    while not os.path.exists('wrfbdy_d01'):
        time.sleep(5)

    # Wait a second to be sure the file is done
    time.sleep(5)
    # Copy the bc file to each directory
    print "Copying BC files..."
    for mem in range(1,Ne+1):
        os.system('cp wrfbdy_d01 mems/m%d/wrfbdy_d01' % mem)



# Make a list of all model dates
startdate = datetime.strptime(date_start,"%Y%m%d%H")
enddate = datetime.strptime(date_end,"%Y%m%d%H")
deltatime = timedelta(minutes=int(fct_len))

opts,args =getopt.getopt(sys.argv[1:],'d:')
for o,a in opts:
    if o == '-d':
        curdate = datetime.strptime(a,"%Y%m%d%H") 
#startdate = datetime(2011,4,17,3)



all_datelist = []
run_datelist = []
while curdate < enddate:
    curdate = curdate + deltatime
    run_datelist.append(datetime.strftime(curdate,"%Y%m%d%H"))

curdate = startdate
while curdate <= enddate:
    all_datelist.append(datetime.strftime(curdate,'%Y%m%d%H'))
    curdate = curdate + deltatime

#print datelist

# Find out which are assimilation times
assim_times = all_datelist[assim_start:N_assim]

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
for date in run_datelist:

    print "Auto-running in silent mode building to date", date
    os.system('./check_ensemble_status.py -d %s -s' % date)
    # Master control lock -- check to see if this file exists                              
    # If it doesn't exist, exit the program                                                
    if not os.path.exists('%s/AUTO_RUN_IN_PROGRESS' % dir_wrf_dom):                        
        print "File 'AUTO_RUN_IN_PROGRESS' is not present in main dir.  Exiting."          
        exit(0)         
    if date in assim_times:
        #os.chdir('%s/wrfdart' % dir_wrf_dom)
        os.system('./submit_filter.py -d %s -m %d' % (date,mpi_numprocs_filter))
        #os.chdir(dir_wrf_dom)

    # Run WPS to get new BCs
    #run_wps_sequence(date,deltatime)

    # Submit the new ensemble
    os.system('./submit_all.py -d %s' % date)
