#!/usr/bin/env python

import os, sys
from WRF_dart_param import Ne,fct_len,dir_members,dir_utils,dir_longsave,queue_members,mpi_numprocs_member,date_start
from datetime import datetime, timedelta
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-d','--d',action='store',dest='startdate',type=str,default=date_start, help='Start date of this cycle')

parser.add_argument('-m','--mems', nargs='+', dest='mems', help='Only submit these members', type=int)


args = parser.parse_args()
startdate = datetime.strptime(args.startdate,'%Y%m%d%H%M%S')
tdelt = timedelta(minutes=int(fct_len))
enddate = startdate + tdelt
print(startdate, enddate)
gzipped_flag = False

submit_mems = args.mems
if submit_mems in ([], None):
    submit_mems = list(range(1,Ne+1))

# Define the restart sequence here
def restart_mem(memnum):
    print("Running diag2wrfinput")
    # Link in the utility if it doesn't exists
    if not os.path.exists('diag2wrfinput'):
        os.system('ln -sf %s/diag2wrfinput' % dir_utils)
    # Link in the old Posterior diag file, if it exists
    if not os.path.exists('%s_Posterior_Diag.nc' % startdate.strftime('%Y%m%d%H')):
        if os.path.exists('%s/%s_Posterior_Diag.nc' % (dir_longsave,startdate.strftime('%Y%m%d%H'))):
            os.system('ln -sf %s/%s_Posterior_Diag.nc %s_Posterior_Diag.nc' % (dir_longsave,\
                       startdate.strftime('%Y%m%d%H'),startdate.strftime('%Y%m%d%H')))
        elif os.path.exists('%s/%s_Posterior_Diag.nc.gz' % (dir_longsave,startdate.strftime('%Y%m%d%H'))):
            print("Unzipping Posterior file")
            gzipped_flag = True
            os.system('gunzip %s/%s_Posterior_Diag.nc.gz' % (dir_longsave, startdate.strftime('%Y%m%d%H')))
            os.system('ln -sf %s/%s_Posterior_Diag.nc %s_Posterior_Diag.nc' % (dir_longsave,\
                     startdate.strftime('%Y%m%d%H'), startdate.strftime('%Y%m%d%H'))) 
        else:
            print("!!!!!! Error !!!!!!!")
            print("Could not find file %s_Posterior_Diag.nc in longsave directory." % startdate.strftime('%Y%m%d%H'))
            exit(1)

    # Check for input.nml
    if not os.path.exists('input.nml'):
        os.system('cp %s/wrfdart/input.nml .' % dir_wrf_dom)

    # Run the utility for this member
    os.system('./diag2wrfinput %s_Posterior_Diag.nc %s %d' % (startdate.strftime('%Y%m%d%H'),\
               startdate.strftime('%Y-%m-%d_%H:%M:%S'), memnum))

    # Copy the new wrfinput file into the member's directory 
    # And clean the directory of previous wrfout files
    os.system('rm %s/m%d/wrfout_d0*' % (dir_members,memnum))
    os.system('cp wrfinput_d01 %s/m%d/wrfinput_d01'%(dir_members, memnum))

    # Clean the current directory
    os.system('rm dart_log.nml') 



def lsf_submit(mem):
    # Write a wrapper c-shell script to 
    # set the environment variables and run on
    # bluefire
    from WRF_dart_param import NCAR_GAU_ACCOUNT, ADVANCE_TIME_MEMBER, ADVANCE_QUEUE_MEMBER, ADVANCE_CORES_MEMBER, NCAR_ADVANCE_PTILE
    with open('run_member_{:02d}.csh'.format(mem), 'w') as outfile:
        outfile.write("#!/bin/csh\n")
        outfile.write("#==================================================================\n")
        outfile.write("#BSUB -J run_member_{:d}\n".format(mem))
        outfile.write("#BSUB -o run_member_{:d}.log\n".format(mem))
        outfile.write("#BSUB -e run_member_{:d}.err\n".format(mem))
        outfile.write("#BSUB -P {:s}\n".format(NCAR_GAU_ACCOUNT))
        outfile.write("#BSUB -W {:s}\n".format(ADVANCE_TIME_MEMBER))
        outfile.write("#BSUB -q {:s}\n".format(ADVANCE_QUEUE_MEMBER))
        outfile.write("#BSUB -n {:d}\n".format(ADVANCE_CORES_MEMBER))
        outfile.write("#BSUB -x\n")
        outfile.write('#BSUB -R "span[ptile={:s}]"\n'.format(NCAR_ADVANCE_PTILE))
        outfile.write("#==================================================================\n")
        outfile.write("\n")
        outfile.write("# Change into the correct directory\n")
        outfile.write("cd {:s}/m{:d}\n".format(dir_members,mem))
        outfile.write("\n")
        outfile.write("# Set the start and end time environment variables\n")
        outfile.write("setenv STARTDATE {:%Y%m%d%H%M%S}\n".format(startdate))
        outfile.write("setenv ENDDATE {:%Y%m%d%H%M%S}\n".format(enddate))
        outfile.write("\n")
        outfile.write("# Run the python run_member script\n")
        outfile.write("{:s}/m{:d}/m{:d}_run_member.py\n".format(dir_members, mem, mem))
    os.system('mv run_member_{:02d}.csh {:s}/m{:d}'.format(mem,dir_members,mem))
    os.system('bsub < {:s}/m{:d}/run_member_{:02d}.csh'.format(dir_members,mem,mem))


def pbs_submit(mem):
    """
    Write a wrapper c-shell script to set environment variables and run
    a simulation with the PBS queue system
    """
    from getpass import getuser
    from WRF_dart_param import NCAR_GAU_ACCOUNT, ADVANCE_TIME_MEMBER, ADVANCE_QUEUE_MEMBER, ADVANCE_CORES_MEMBER, NCAR_ADVANCE_PTILE, nnodes_member, numprocs_per_node
    with open('run_member_{:02d}.csh'.format(mem), 'w') as outfile:
        outfile.write('#!/bin/csh\n')
        outfile.write('#PBS -N run_member_{:d}\n'.format(mem))
        outfile.write('#PBS -A {:s}\n'.format(NCAR_GAU_ACCOUNT))
        outfile.write('#PBS -l walltime={:s}\n'.format(ADVANCE_TIME_MEMBER))
        outfile.write('#PBS -q {:s}\n'.format(ADVANCE_QUEUE_MEMBER))
        #outfile.write('#PBS -j oe\n')
        #outfile.write('#PBS -m abe\n')
        #outfile.write('#PBS -M {:s}@ucar.edu\n'.format(getuser()))
        outfile.write('#PBS -l select={:d}:ncpus={:d}:mpiprocs={:d}\n'.format(nnodes_member, numprocs_per_node, numprocs_per_node))
        outfile.write("#==================================================================\n")
        outfile.write("\n")
        outfile.write("# Change into the correct directory\n")
        outfile.write("cd {:s}/m{:d}\n".format(dir_members,mem))
        outfile.write("\n")
        outfile.write("# Set the start and end time environment variables\n")
        outfile.write("setenv STARTDATE {:%Y%m%d%H%M%S}\n".format(startdate))
        outfile.write("setenv ENDDATE {:%Y%m%d%H%M%S}\n".format(enddate))
        outfile.write("\n")
        outfile.write("# Run the python run_member script\n")
        outfile.write("{:s}/m{:d}/m{:d}_run_member.py\n".format(dir_members, mem, mem))
    os.system('mv run_member_{:02d}.csh {:s}/m{:d}'.format(mem,dir_members,mem))
    os.system('qsub < {:s}/m{:d}/run_member_{:02d}.csh'.format(dir_members,mem,mem))


def submit_all():
    # FIRST--set the environment variables correctly
    os.putenv('STARTDATE',startdate.strftime('%Y%m%d%H'))
    os.putenv('ENDDATE',enddate.strftime('%Y%m%d%H'))
    # Clear out the log files
    os.system('rm -f *.log')
    os.system('rm -f *.err')

    # Loop through each member
    for mem in submit_mems:
    #for mem in range(1,int(Ne)+1):
    #for mem in [46]:
    #for mem in [6,16,27,28,39]:
        print("Member %d..." % (mem),)
        # Copy the run_member script as run_member_%d.py in each directory
        os.system('cp run_member.py %s/m%d/m%d_run_member.py' % (dir_members,mem,mem))


        # Submit this member to the queue
        # Check for which system we are on
        if os.uname()[0] == 'AIX' and os.uname()[1].startswith('be'):
            # We're on bluefire
            lsf_submit(mem)
        elif os.uname()[1].startswith('ys'):
            # We're on Yellowstone.  Use the Bluefire submission.
            lsf_submit(mem)
        else:
            pbs_submit(mem)
            #os.system('qsub -pe ompi %d -q %s -V -o %s/m%d -e %s/m%d -wd %s/m%d %s/m%d/m%d_run_member.py' \
            #       % (mpi_numprocs_member,queue_members,dir_members,mem,dir_members,mem,dir_members,mem,\
            #         dir_members,mem,mem))
        print("Done.")


if __name__ == '__main__':
    submit_all()
