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
from ens_dart_param import *
from optparse import OptionParser
from namelist_utils import read_namelist, write_dart_namelist

parser = OptionParser(description=desc)

parser.add_option('-d','--datein',dest='datein',action='store',type='string',default=cycle_len,\
                  help='Date of assimilation cycle (YYYYMMDDHH)')
parser.add_option('-m','--mpi_procs',dest='mpi_procs',action='store',type='string',\
                  default=mpi_numprocs_filter, help='Number of processors to use for MPI run')

(opts,args) = parser.parse_args()
#datein = datetime.strptime(opts.datein,'%Y%m%d%H')
#prevdate = datein - timedelta(minutes=int(fct_len))
datein = int(opts.datein)
prevdate = datein - cycle_len
mpi_numprocs_filter = int(opts.mpi_procs)
# Get the namelist values
nmld = read_namelist('input.nml')
cm1nmld = read_namelist('namelist.input')

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
        # We've got the files, so now copy in the observation file
        obfile = os.path.join(dir_obs, '{:d}_obs_seq.prior'.format(datein))
        if os.path.exists(obfile):   
            os.system('cp {:s} obs_seq.prior'.format(obfile))
        else:
            error_handler('Could not find {:s}'.format(obfile),\
                          'submit_filter')

    
        # Make sure latest filter is linked in
        #os.system('rm -f filter')
        #os.system('ln -sf {:s} filter'.format(os.path.join(dir_src_dart,'filter')))
        
        # Make sure template is linked in
        if not os.path.exists('./cm1out_rst_000001.nc'):
            os.system('cp {:s}/cm1out_rst_000001.nc .'.format(dir_dom))
        if not os.path.exists('./namelist.input'):
            os.system('cp {:s}/namelist.input .'.format(dir_dom))

        # Now write the submission script
        qsub_cmd, scriptname = write_filter_submit()

        # Using the script and command provided, submit the filter
        if os.path.exists('filter_done'):
            os.system('rm filter_done')
 

    if PRE_CHECK and not flag_direct_netcdf_io:
        # Now, check to be sure all filter_ic_old.#### files are in place
        for r in range(1,Ne+1):
            if not os.path.exists('filter_ic_old.%04d' % r):
                error_handler('Could not find filter_ic_old.%04d. Exiting.' % r, 'submit_filter')
    else:
        with open('input_filelist.txt','w') as filelist:
            for m in xrange(1,Ne+1):
                filelist.write('{:s}/m{:d}/cm1out_rst_000001.nc\n'.format(dir_members, m))


    if RUN_FILTER:
        os.system('touch dart_log.out')
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
    curdate = start + (assim_dt * (assim_start))
    no_adaptive_inf_dates = []
    step = 1
    # WRF dart param tells us which assimilation step to start using inflation

    while step < inflate_start:
        no_adaptive_inf_dates.append(curdate)
        curdate = curdate + assim_dt
        step = step + 1
    print(no_adaptive_inf_dates)
    # Now check to see if we are in the list of no inflation dates
    if datein not in no_adaptive_inf_dates: 
        print("Using adaptive inflation values from previous time")
        nmld['filter_nml']['inf_initial_from_restart'] = [True, True]
        nmld['filter_nml']['inf_sd_initial_from_restart'] = [True, True]
        prior_inf_mean = os.path.join(dir_longsave, '{:06d}_inf_ic_mean.nc'.format(prevdate))
        prior_inf_sd = os.path.join(dir_longsave, '{:06d}_inf_ic_sd.nc'.format(prevdate))
        if os.path.exists(prior_inf_mean):
            os.system('cp {:s} {:s}/prior_inf_ic_old_mean.nc'.format(prior_inf_mean, dir_assim))
            os.system('cp {:s} {:s}/prior_inf_ic_old_sd.nc'.format(prior_inf_sd, dir_assim))
        else:
            error_handler('Could not find {:s}'.format(prior_inf_file),'submit_filter') 
    else:
        nmld['filter_nml']['inf_initial_from_restart'] = [False, False]
        nmld['filter_nml']['inf_sd_initial_from_restart'] = [False, False]
        print("Using initial values for inflation mean and std")
    
    # now figure out the "date" we are at based on the CM1 namelist
    cm1date = cm1nmld['param11']
    start_date = datetime(cm1date['year'], cm1date['month'], cm1date['day'], cm1date['hour'],\
                          cm1date['minute'], cm1date['second'])
   
    # Make sure num members is write
    nmld['filter_nml']['ens_size'] = Ne
    write_dart_namelist(nmld, date=start_date + timedelta(seconds=datein))
    os.system('cp input.nml {:s}/input.nml'.format(dir_assim))


def clean_dart(new_flag, old_flag):
    os.system('rm -f Posterior_Diag.nc')
    os.system('rm -f Prior_Diag.nc')
    os.system('rm -f PriorDiag*')
    os.system('rm -f mean_d01.nc')
    os.system('rm -f sd_d01.nc')
    os.system('rm -f *_forward_op_errors')
    os.system('rm -f assim_model*')
    os.system('rm -f *.out')
    os.system('rm -f prior_inf_ic_old')
    os.system('rm -f prior_member*.nc')
    # If new flag is given, remove filter_ic_new.*
    if new_flag:
        os.system('rm -f filter_ic_new.*')
    if old_flag:
        os.system('rm -f filter_ic_old.*')




def archive_files(datem):
    print("#################### ARCHIVING FILES #######################")
    # Run diagnostics on the posterior file
    os.system('ln -sf obs_seq.posterior obs_seq.diag')
    #if not os.path.exists('obs_diag'):
    #    os.system('ln -sf {:s}/obs_diag .'.format(dir_src_dart))
    os.system('{:s}/obs_diag'.format(dir_src_dart))
    if not os.path.exists('obs_diag_output.nc'):
        print("Failure to produce obs_diag_output.nc!")
        pass
    else:
        os.system('mv obs_diag_output.nc {:s}/{:06d}_obs_diag_output.nc'.format(dir_longsave,datem))
    os.system('unlink obs_seq.diag')

    # Convert the posterior obs sequence file to netcdf format for easier diagnosis later
    #if not os.path.exists('obs_seq_to_netcdf'):
    #    os.system('ln -sf {:s}/obs_seq_to_netcdf .'.format(dir_src_dart))
    os.system('{:s}/obs_seq_to_netcdf'.format(dir_src_dart))
    if not os.path.exists('obs_epoch_001.nc'):
        print("Failure to produce obs_epoch_001.nc!")
        pass
    else:
        os.system('mv obs_epoch_001.nc {:s}/{:06d}_obs_sequence.nc'.format(dir_longsave,datem))



    # Move obs_seq.prior and obs_seq.posterior file to longsave
    os.system('mv -f obs_seq.prior  {:s}/{:06d}_obs_seq.prior'.format(dir_longsave,datem))
    os.system('mv -f obs_seq.posterior {:s}/{:06d}_obs_seq.posterior'.format(dir_longsave,datem))

    # LEM -- Revisions here for new diag format, with just the mean and sd
    if os.path.exists('prior_inf_ic_new_mean.nc'):
        os.system('mv -f prior_inf_ic_new_mean.nc {:s}/{:06d}_inf_ic_mean.nc'.format(dir_longsave,datem))
        os.system('mv -f prior_inf_ic_new_sd.nc {:s}/{:06d}_inf_ic_sd.nc'.format(dir_longsave,datem))
    if os.path.exists('mean.nc') and os.path.exists('sd.nc'):
        os.system('ncdiff mean.nc PriorDiag_mean.nc mean_increment.nc')
        os.system('ncdiff sd.nc PriorDiag_sd.nc sd_increment.nc')
        os.system('mv -f mean_increment.nc {:s}/{:06d}_mean_increment.nc'.format(dir_longsave,datem))
        os.system('mv -f sd_increment.nc {:s}/{:06d}_sd_increment.nc'.format(dir_longsave,datem))
        os.system('mv -f mean.nc {:s}/{:06d}_mean.nc'.format(dir_longsave,datem))
        os.system('mv -f sd.nc {:s}/{:06d}_sd.nc'.format(dir_longsave,datem))


    # Check to see if we are compressing the Diag files
    if flag_compress_diag:
        # Zip up the files
        curdir = os.getcwd()
        os.chdir(dir_longsave)
        #os.system('gzip -f %s_Prior_Diag.nc' % datem)
        os.system('gzip -f {:06d}_mean.nc'.format(datem))
        os.system('gzip -f {:06d}_sd.nc'.format(datem))
        os.system('gzip -f {:06d}_mean_increment.nc'.format(datem))
        os.system('gzip -f {:06d}_sd_increment.nc'.format(datem))
        os.chdir(curdir)







def write_filter_submit():
    # Function to determine which system we are on and write accordingly
    # Three possibilities now -- enkf,student cluster or bluefire

    if mpi_numprocs_filter == 1:
         # We're not requesting an mpirun, return just the filter command 
         return ('{:s}/filter'.format(dir_src_dart),'')
    
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
            outfile.write("#!/bin/csh\n")
            outfile.write("#==================================================================\n")
            outfile.write("#BSUB -J run_filter\n")
            outfile.write("#BSUB -o submit_filter.%J.log\n")
            outfile.write("#BSUB -e submit_filter.%J.err\n")
            outfile.write("#BSUB -P {:s}\n".format(NCAR_GAU_ACCOUNT))
            outfile.write("#BSUB -W {:s}\n".format(ADVANCE_TIME_FILTER))
            outfile.write("#BSUB -q {:s}\n".format(ADVANCE_QUEUE_FILTER))
            outfile.write("#BSUB -n {:d}\n".format(ADVANCE_CORES_FILTER))
            outfile.write("#BSUB -x\n")
            outfile.write('#BSUB -R "span[ptile={:s}]"\n'.format(NCAR_ADVANCE_PTILE))
            outfile.write("#==================================================================\n")
            outfile.write("limit stacksize unlimited\n")
            outfile.write("setenv OMP_STACKSIZE 200000000000\n")
            outfile.write("setenv MP_STACK_SIZE 200000000000\n")
            outfile.write("set start_time = `date +%s`\n")
            outfile.write('echo "host is " `hostname`\n')
            outfile.write("\n")
            outfile.write("cd {:s}\n".format(dir_dom))
            outfile.write("echo $start_time >& {:s}/filter_started\n".format(dir_dom))
            outfile.write("\n")
            outfile.write("#  run data assimilation system\n")
            outfile.write("setenv TARGET_CPU_LIST -1\n")
            outfile.write("mpirun.lsf job_memusage.exe ./filter\n")
            outfile.write("\n")
            outfile.write("touch {:s}/filter_done\n".format(dir_dom))
            outfile.write("set end_time = `date  +%s`\n")
            outfile.write("@ length_time = $end_time - $start_time\n")
            outfile.write('echo "duration = $length_time"\n')
        return ('bsub < ','run_filter_mpi.csh')

    else:
        # We're on a UW system
        if os.path.exists('run_filter_mpi.py'):
            os.system('rm run_filter_mpi.py')
        with open('run_filter_mpi.py','w') as outfile:
            outfile.write("#!/usr/bin/env python\n")
            outfile.write("\n")
            outfile.write("import os\n")
            outfile.write("# Change to directory\n" )
            curdir = os.getcwd()
            outfile.write("os.chdir('{:s}')\n".format(curdir))
            outfile.write("os.system('{:s} -np {:d} {:s}/filter >> filter.out')\n".format(mpi_run_command,mpi_numprocs_filter, dir_src_dart))
            outfile.write("os.system('touch filter_done')\n")
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
