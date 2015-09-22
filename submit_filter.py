#!/usr/bin/env python
from __future__ import print_function, division
desc = """
 SUBMIT_FILTER.PY
 Version 2 of Python WRF/DART System
 Written by Luke Madaus, July 2012
###########################################\n
 Function main() represents the
 sequence to execute the filter
 This script is designed to be used in conjuction
 with the run_member.py script
 As such, all wrf_to_dart and dart_to_wrf details
 are handled within run_member
"""

import os
import time
from datetime import datetime, timedelta
from WRF_dart_param import *
from optparse import OptionParser

parser = OptionParser(description=desc)

parser.add_option('-d','--datein',dest='datein',action='store',type='string',default=date_start,\
                  help='Date of assimilation cycle (YYYYMMDDHH)')
parser.add_option('-m','--mpi_procs',dest='mpi_procs',action='store',type='string',\
                  default=mpi_numprocs_filter, help='Number of processors to use for MPI run')

(opts,args) = parser.parse_args()
#datein = datetime.strptime(opts.datein,'%Y%m%d%H')
#prevdate = datein - timedelta(minutes=int(fct_len))
datein = int(opts.datein)
prevdate = datein - cycle_len
datein *= 60
prevdate *= 60
mpi_numprocs_filter = int(opts.mpi_procs)

# Change these values to only run certain sections of the code
PRE_CLEAN       = True
PRE_CHECK       = True
RUN_FILTER      = True
ARCHIVE_FILES   = True
POST_CLEAN      = True


def main():
 
    # First step--write the namelist
    # Need to know if we are trying to apply adaptive
    # inflation from the previous time or not.  If we are using it,
    # Copy in and unzip the previous inflation
    make_namelist_and_inflate(datein,prevdate)
    
    # Put us in the assimilation directory
    os.chdir(dir_assim)

    if PRE_CLEAN:
        # Start by cleaning up the DART directoy,
        # Only removing old filter_ic_new files
        clean_dart(True, False)
 

    if PRE_CHECK and not flag_direct_netcdf_io:
        # Now, check to be sure all filter_ic_old.#### files are in place
        for r in range(1,Ne+1):
            if not os.path.exists('filter_ic_old.%04d' % r):
                error_handler('Could not find filter_ic_old.%04d. Exiting.' % r, 'submit_filter')

    if RUN_FILTER:
        # We've got the files, so now copy in the observation file
        obfile = os.path.join(dir_obs, '{:d}_obs_seq.prior'.format(datein))
        if os.path.exists(obfile):   
            os.system('cp {:s} obs_seq.prior'.format(obfile))
        else:
            error_handler('Could not find {:s}'.format(obfile),\
                          'submit_filter')

    
        # Make sure latest filter is linked in
        os.system('rm filter')
        os.system('ln -sf {:s} filter'.format(os.path.join(dir_src_dart,'filter')))
        
        # Make sure wrfinput is linked in
        #for dom in range(1,max_dom+1):
        #    if not os.path.exists('./wrfinput_d{:02d}'.format(dom)):
        #        os.system('ln -sf {:s}/wrfinput_d{:02d} .'.format(dir_wrf_dom, dom))

        # Now write the submission script
        qsub_cmd, scriptname = write_filter_submit(datein)

        # Using the script and command provided, submit the filter
        if os.path.exists('filter_done'):
            os.system('rm filter_done')

        # Setting times will let us see how long it is taking
        t_0 = time.time()
        os.system('{:s} {:s}'.format(qsub_cmd,scriptname))
    
        # Now sleep while waiting for the filter to finish
        while not os.path.exists('filter_done'):
            time.sleep(5)   

        # Only continue once filter_done is found
        print("Filter finished!")
        os.system('rm filter_done')
        t_1 = time.time()
        print("Filter execution time:", t_1 - t_0)

        # Check to see if obs_seq.posterior has been created
        if not os.path.exists('obs_seq.posterior'):
            error_handler('Could not find obs_seq.posterior.  Problem.','submit_filter')
    

    if ARCHIVE_FILES:
        # Conclude by archiving files and cleaning, removing just old filter_ic files
        archive_files(datein)
    if POST_CLEAN:
        clean_dart(False, True)


def make_namelist_and_inflate(datein, prevdate):
    # Determine if we need to copy in the previous inflation
    # values.  If so, do it and unzip them.  If not, just
    # write the namelist.

    # Use inflate_start to find the times we don't need
    # Make a list of all assim times
    #start = datetime.strptime(date_start,'%Y%m%d%H')
    #assim_dt = timedelta(minutes=int(fct_len))
    start = 0
    assim_dt = int(cycle_len)


    # Start on the first actual assimilation time
    curdate = start + (assim_dt * (assim_start-1))
    no_adaptive_inf_dates = []
    step = 1
    # WRF dart param tells us which assimilation step to start using inflation

    while step < inflate_start:
        no_adaptive_inf_dates.append(curdate)
        curdate = curdate + assim_dt
        step = step + 1 

    # Now check to see if we are in the list of no inflation dates
    if datein not in no_adaptive_inf_dates: 
        print("Using adaptive inflation values from previous time")
        os.system('./make_namelist_dart.py -d {:d} -i'.format(datein/60))
        prior_inf_file = os.path.join(dir_longsave, '{:d}_prior_inf_ic.gz'.format(prevdate))
        if os.path.exists(prior_inf_file):
            os.system('cp {:s} {:s}/prior_inf_ic_old.gz'.format(prior_inf_file, dir_assim))
            os.system('gunzip prior_inf_ic_old.gz')
        else:
            error_handler('Could not find {:s}'.format(prior_inf_file),'submit_filter') 
    else:
        print("Using initial values for inflation mean and std")
        os.system('./make_namelist_dart.py -d {:d}'.format(datein))
    os.system('cp input.nml {:s}/input.nml'.format(dir_assim))


def clean_dart(new_flag, old_flag):
    os.system('rm Posterior_Diag.nc')
    os.system('rm Prior_Diag.nc')
    os.system('rm PriorDiag*')
    os.system('rm mean_d01.nc')
    os.system('rm sd_d01.nc')
    os.system('rm *_forward_op_errors')
    os.system('rm assim_model*')
    os.system('rm *.out')
    os.system('rm prior_inf_ic_old')
    # If new flag is given, remove filter_ic_new.*
    if new_flag:
        os.system('rm filter_ic_new.*')
    if old_flag:
        os.system('rm filter_ic_old.*')




def archive_files(datem):
    print("#################### ARCHIVING FILES #######################")
    # Run diagnostics on the posterior file
    os.system('ln -sf obs_seq.posterior obs_seq.diag')
    if not os.path.exists('obs_diag'):
        os.system('ln -sf {:s}/obs_diag .'.format(dir_src_dart))
    os.system('./obs_diag')
    if not os.path.exists('obs_diag_output.nc'):
        print("Failure to produce obs_diag_output.nc!")
        pass
    else:
        os.system('mv obs_diag_output.nc {:s}/{:06d}_obs_diag_output.nc'.format(dir_longsave,datem))
    os.system('unlink obs_seq.diag')

    # Convert the posterior obs sequence file to netcdf format for easier diagnosis later
    if not os.path.exists('obs_seq_to_netcdf'):
        os.system('ln -sf {:s}/obs_seq_to_netcdf .'.format(dir_src_dart))
    os.system('./obs_seq_to_netcdf')
    if not os.path.exists('obs_epoch_001.nc'):
        print("Failure to produce obs_epoch_001.nc!")
        pass
    else:
        os.system('mv obs_epoch_001.nc {:s}/{:06d}_obs_sequence.nc'.format(dir_longsave,datem))



    # Move obs_seq.prior and obs_seq.posterior file to longsave
    os.system('mv -f obs_seq.prior  {:s}/{:06d}_obs_seq.prior'.format(dir_longsave,datem))
    os.system('mv -f obs_seq.posterior {:s}/{:06d}_obs_seq.posterior'.format(dir_longsave,datem))

    # Move the prior and posterior inflation files to longsave
    """
    if os.path.exists('prior_inf_ic_new'):
        os.system('cp -f prior_inf_ic_new prior_inf_ic_old')
        os.system('gzip -f prior_inf_ic_new')
        os.system('mv -f prior_inf_ic_new.gz   %s/%s_prior_inf_ic.gz' % (dir_longsave, datem))
    if os.path.exists('post_inf_ic_new'):
        os.system('cp -f post_inf_ic_new   post_inf_ic_old')
        os.system('gzip -f post_inf_ic_new')
        os.system('mv -f post_inf_ic_new.gz   %s/%s_post_inf_ic.gz' % (dir_longsave, datem))

    # Move the prior and posterior inflation diag files to longsave
    if os.path.exists('prior_inf_diag'):
        os.system('gzip -f prior_inf_diag')
        os.system('mv -f prior_inf_diag.gz   %s/%s_prior_inf_diag.gz' % (dir_longsave, datem))
    if os.path.exists('post_inf_diag'):
        os.system('gzip -f post_inf_diag')
        os.system('mv -f post_inf_diag.gz   %s/%s_post_inf_diag.gz' % (dir_longsave, datem))


    # Move Prior and Posterior Diag files to longsave
    os.system('mv -f Prior_Diag.nc  %s/%s_Prior_Diag.nc' % (dir_longsave, datem))
    os.system('mv -f Posterior_Diag.nc %s/%s_Posterior_Diag.nc' % (dir_longsave, datem))
    os.system('
    """
    # LEM -- Revisions here for new diag format, with just the mean and sd
    if os.path.exists('prior_inf_ic_new_mean_d01'):
        os.system('mv -f prior_inf_ic_new_mean_d01 {:s}/{:06d}_prior_inf_ic_mean_d01'.format(dir_longsave,datem))
        os.system('mv -f prior_inf_ic_new_sd_d01 {:s}/{:06d}_prior_inf_ic_sd_d01'.format(dir_longsave,datem))
    if os.path.exists('mean_d01.nc') and os.path.exists('sd_d01.nc'):
        os.system('ncdiff mean_d01.nc PriorDiag_mean_d01.nc mean_increment.nc')
        os.system('ncdiff sd_d01.nc PriorDiag_sd_d01.nc sd_increment.nc')
        os.system('mv -f mean_increment.nc {:s}/{:06d}_mean_increment.nc'.format(dir_longsave,datem))
        os.system('mv -f sd_increment.nc {:s}/{:06d}_sd_increment.nc'.format(dir_longsave,datem))


    # Check to see if we are compressing the Diag files
    if flag_compress_diag:
        # Zip up the files
        curdir = os.getcwd()
        os.chdir(dir_longsave)
        #os.system('gzip -f %s_Prior_Diag.nc' % datem)
        os.system('gzip -f {:06d}_Posterior_Diag.nc'.format(datem))
        if os.path.exists('{:s}/{:06d}_Prior_Diag_int.nc'.format(dir_longsave, datem)):
            os.system('gzip -f {:06d}_Prior_Diag_int.nc'.format(datem))
            os.system('gzip -f {:06d}_Posterior_Diag_int.nc'.format(datem))
        os.chdir(curdir)







def write_filter_submit():
    # Function to determine which system we are on and write accordingly
    # Three possibilities now -- enkf,student cluster or bluefire

    if mpi_numprocs_filter == 1:
         # We're not requesting an mpirun, return just the filter command 
         return ('./filter','')
    
    node_name = os.uname()[1]


    if node_name.startswith('be') or node_name.startswith('ys'):
        # We're on bluefire or yellowstone
        print("Submitting on YELLOWSTONE")                                                               
        # Import special variables                                                                      
        from WRF_dart_param import NCAR_GAU_ACCOUNT, ADVANCE_TIME_FILTER, ADVANCE_QUEUE_FILTER, ADVANCE_CORES_FILTER, NCAR_ADVANCE_PTILE
        if os.path.exists('run_filter_mpi.csh'):
            os.system('rm run_filter_mpi.csh')
        # Write a new run_filter_mpi.csh
        with open('run_filter_mpi.csh','w') as outfile:
            outfile.write("#!/bin/csh")
            outfile.write("#==================================================================")
            outfile.write("#BSUB -J run_filter")
            outfile.write("#BSUB -o submit_filter.%J.log")
            outfile.write("#BSUB -e submit_filter.%J.err")
            outfile.write("#BSUB -P {:s}".format(NCAR_GAU_ACCOUNT))
            outfile.write("#BSUB -W {:s}".format(ADVANCE_TIME_FILTER))
            outfile.write("#BSUB -q {:s}".format(ADVANCE_QUEUE_FILTER))
            outfile.write("#BSUB -n {:d}".format(ADVANCE_CORES_FILTER))
            outfile.write("#BSUB -x")
            outfile.write('#BSUB -R "span[ptile={:s}]"'.format(NCAR_ADVANCE_PTILE))
            outfile.write("#==================================================================")
            outfile.write("limit stacksize unlimited")
            outfile.write("setenv OMP_STACKSIZE 200000000000")
            outfile.write("setenv MP_STACK_SIZE 200000000000")
            outfile.write("set start_time = `date +%s`")
            outfile.write('echo "host is " `hostname`')
            outfile.write("")
            outfile.write("cd {:s}".format(dir_dom))
            outfile.write("echo $start_time >& {:s}/filter_started".format(dir_dom))
            outfile.write("")
            outfile.write("#  run data assimilation system")
            outfile.write("setenv TARGET_CPU_LIST -1")
            outfile.write("mpirun.lsf job_memusage.exe ./filter")
            outfile.write("")
            outfile.write("touch {:s}/filter_done".format(dir_dom)
            outfile.write("set end_time = `date  +%s`")
            outfile.write("@ length_time = $end_time - $start_time")
            outfile.write('echo "duration = $length_time"')
        return ('bsub < ','run_filter_mpi.csh')

    else:
        # We're on a UW system
        if os.path.exists('run_filter_mpi.py'):
            os.system('rm run_filter_mpi.py')
        with open('run_filter_mpi.py','w') as outfile:
            outfile.write("#!/usr/bin/env python")
            outfile.write("")
            outfile.write("import os")
            outfile.write("# Change to directory" )
            curdir = os.getcwd()
            outfile.write("os.chdir('{:s}')".format(curdir))
            outfile.write("os.system('{:s} -np {:d} filter >> filter.out')".format(mpi_run_command,mpi_numprocs_filter))
            outfile.write("os.system('touch filter_done')")
        return ('qsub -pe ompi {:d} -V -q {:s} -e {:s} -o {:s}'.format(mpi_numprocs_filter,queue_filter,dir_assim,dir_assim),\
                'run_filter_mpi.py')

    #else:
    #    error_handler('Unable to determine which system we are on. Only enkf,student cluster and bluefire are supported now.',\
    #                   'write_filter_submit')




def error_handler(msg,routine):
    print('!!!!!! Error in routine {:s} !!!!!!'.format(routine))
    print(msg)
    exit(1)



if __name__ == '__main__':
    main()
