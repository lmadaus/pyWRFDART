#!/usr/bin/env python

import os, sys, getopt
from datetime import datetime, timedelta
from WRF_dart_param import *
from namelist_utils import read_namelist, write_namelist, update_time_wrf, update_time_wps
# Quick grab from WRF_dart_param

time = int(date_start)
timedt = datetime.strptime(date_start, '%Y%m%d%H')
"""
This is a quickly-reformatted version of the make_ensemble script
from version 1 of this code that makes it compatible with the new
namelist version.  This script will make an ensemble (in serial) with
a number of members as specified in the namelist variable "Ne" or manually
specified from the command line.  The script will copy in all of the necessary
files and perturb the initial conditions for each member using WRFDA.

"""


(opts,args) = getopt.getopt(sys.argv[1:],'n:d:')
for o,a in opts:
   if o == '-n':
        Ne = int(a)
   if o == '-d':
        time = int(a)
        timedt = datetime.strptime(a, '%Y%m%d%H')

# And the end time
enddt = datetime.strptime(date_end, '%Y%m%d%H')


def main():
    # WRFRUNDIR and WRFVARDIR are now in WRF_dart_param
    DOMDIR = dir_wrf_dom
    os.chdir(DOMDIR)



    print("############### MAKING ENSEMBLE ###############")
    print("############## NUMBER OF MEMS: %d #############" % Ne)
    # Now make the ensemble directory
    if not os.path.exists('mems'):
        print("Making ensemble directory (mems)")
        os.system('mkdir mems')  
    os.chdir(DOMDIR+'/mems')


    # Make a mean directory and copy the input file there
    os.system('mkdir logs')
    os.system('mkdir mean')
    os.chdir('mean')
    # The wrfinput files for each domain and a wrfbdy_d01 file
    # must be present in dir_wrf_dom to be copied to each member
    for dom in range(int(max_dom)+1)[1:]:
        if os.path.exists('%s/wrfinput_d0%d' % (DOMDIR,dom)) and os.path.exists('%s/wrfbdy_d01' % (DOMDIR)): 
            os.system('cp %s/wrfinput_d0%d .' % (DOMDIR,dom))
        else:
            print("Error!")
            print("Unable to find wrfinput or wrfbdy file in main directory!")
            exit(1)
    os.system('cp %s/wrfbdy_d01 .' % DOMDIR)


    # Loop through each member
    print("######################################")
    print("Making ensemble member sub-directories")
    print("and perturbing initial conditions")
    print("######################################")
    print("")
    # Read in the namelist data
    nmlinput = read_namelist(os.path.join(DOMDIR,'namelist.input'))
    nmlwps = read_namelist(os.path.join(DOMDIR, 'namelist.wps'))
    # Update the time in each
    nmlinput = update_time_wrf(nmlinput, timedt, enddt)
    nmlwps = update_time_wps(nmlwps, timedt, enddt)

    
    # Get all WRF tables to link in
    tablefiles = ['TBL','DBL','DATA','formatted']
    tables = [f for f in os.listdir(WRFRUNDIR) if any(substr in f for substr in tablefiles)]
    #print(tables)
    for k in range(Ne):
    #for k in [0]:
        k = k+1 # Don't start at 0
        print("***** MEMBER %d *****" %k)
        os.chdir('%s' % DOMDIR)
        # hold onto that while making a directory for it
        os.chdir('%s/mems' % DOMDIR)
        os.system('mkdir %s/mems/m%d' % (DOMDIR,k))
        os.chdir('%s/mems/m%d' % (DOMDIR,k))
        # Run script to generate the namelist
        #os.system('./write_namelists.py -d %d -e %d -f' % (time,k))
        write_namelist(nmlinput, 'namelist.input')
        write_namelist(nmlwps, 'namelist.wps')
        # That will be written as namelist.input
        # Copy in the namelist
        #os.system('cp %s/namelist.input namelist.input' % DOMDIR)

        # Copy the input file as fg
        print("Copying in wrfinput and wrfbdy files...")
        for dom in range(int(max_dom)+1)[1:]:
            os.system('cp %s/mems/mean/wrfinput_d0%d .' % (DOMDIR,dom))
            #os.system('cp /glade/scratch/schwartz/stella/ensf/2017031400/wrf_rundir/ens_{:d}/wrfinput_d{:02d} .'.format(k,dom))
        os.system('cp /glade/scratch/schwartz/stella/ensf/2017031400/wrf_rundir/ens_{:d}/wrfout_d02_2017-03-14_00:00:00.gz .'.format(k))
        os.system('gunzip wrfout_d02_2017-03-14_00:00:00.gz')
        os.system('mv wrfout_d02_2017-03-14_00:00:00 wrfinput_d02')
        #os.system('mv wrfinput_d01 fg')
        # Copy in the wrfbdy file
        #os.system('cp /glade/scratch/schwartz/stella/ensf/2017031400/wrf_rundir/ens_{:d}/wrfbdy_d01 .'.format(k))
        os.system('cp {:s}/wrfbdy_d01 .'.format(DOMDIR))
        # Link in executables
        for f in tables:
            os.system('ln -sf {:s} .'.format(os.path.join(WRFRUNDIR, f)))
        #os.system('ln -sf %s/var/run/be.dat.cv3 be.dat' % WRFVARDIR)
        #os.system('ln -sf %s/var/build/da_wrfvar.exe .' % WRFVARDIR)
        os.system('ln -sf %s/wrf.exe .' % dir_src_wrf)
        """
        if int(max_dom) > 1:
            os.system('ln -sf %s/ndown.exe .' % dir_src_wrf)

        # Now execute WRFVAR for first domain
        print("Running WRFVAR to perturb domain 1")
        os.system('./da_wrfvar.exe')
        os.system('%s ./da_wrfvar.exe' % mpi_run_command)
        # See if output file exists
        if os.path.exists('wrfvar_output'):
            os.system('mv wrfvar_output wrfinput_d01')
        else:
            print("WRFVAR execution failed for ensmem %d." %k)
            exit(1)
   

        # If we have other domains, run ndown to re-perturb inner model
        if int(max_dom) > 1:
            # Re-run script to generate the namelists with proper number of domains
            os.chdir(DOMDIR)
            os.system('./write_namelists.py -d %d -e %d -l 0' % (time,k))
            #os.system('./write_wrfvar.py -d %d -e %d' % (time,k))
            os.system('cp namelist.input %s/mems/m%d/namelist.input' % (DOMDIR,k))
            os.chdir('%s/mems/m%d' % (DOMDIR,k))
            os.system('cp wrfinput_d01 wrfout_d01')
            for dom in range(int(max_dom)+1)[2:]:
                os.system('mv wrfinput_d0%d wrfndi_d02' % dom)
                print("Running NDOWN to interp. perturbed dom 1 to dom %d" % dom)
                os.system('ulimit -s unlimited; ulimit -d unlimited; ./ndown.exe >! ndown.out')
                if os.path.exists('wrfinput_d02'):
                    os.system('mv wrfinput_d02 wrfinput_d0%d' % dom)
                    os.system('rm wrfndi_d02')
                else:
                    print("Error. NDOWN failed for domain %d in member %d." % (dom,k))
                    exit(1)

        """

        # Save the initial wrfin files in case of restart
        #for dom in range(int(max_dom)+1)[1:]:
        #    os.system('cp wrfinput_d%02d wrfinput_d%02d_orig' % (dom,dom))


        # Remove old links
        os.system('rm be.dat da_wrfvar.exe')
        # Go back up a directory to start the next one
        os.chdir('%s/mems' % DOMDIR)

    print("### DONE MAKING DIRECTORIES w/ ICs ###")
    print("")

    os.chdir(dir_wrf_dom)



if __name__ == '__main__':
    main()
