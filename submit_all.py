#!/usr/bin/env python

import os, sys
from WRF_dart_param import Ne,fct_len,dir_members,dir_utils,dir_longsave,queue_members,mpi_numprocs_member,date_start
from datetime import datetime, timedelta
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-d','--d',action='store',dest='startdate',type='string',default=date_start,
                  help='Start date of this cycle')

parser.add_option('-r','--restart',action='store_true',dest='restart_flag',\
                  default=False, help='Submit this cycle as a restart')


(opts,args) = parser.parse_args()
startdate = datetime.strptime(opts.startdate,'%Y%m%d%H')
tdelt = timedelta(minutes=int(fct_len))
enddate = startdate + tdelt

restart_flag = opts.restart_flag
gzipped_flag = False


# Define the restart sequence here
def restart_mem(memnum):
    print "Running diag2wrfinput"
    # Link in the utility if it doesn't exists
    if not os.path.exists('diag2wrfinput'):
        os.system('ln -sf %s/diag2wrfinput' % dir_utils)
    # Link in the old Posterior diag file, if it exists
    if not os.path.exists('%s_Posterior_Diag.nc' % startdate.strftime('%Y%m%d%H')):
        if os.path.exists('%s/%s_Posterior_Diag.nc' % (dir_longsave,startdate.strftime('%Y%m%d%H'))):
            os.system('ln -sf %s/%s_Posterior_Diag.nc %s_Posterior_Diag.nc' % (dir_longsave,\
                       startdate.strftime('%Y%m%d%H'),startdate.strftime('%Y%m%d%H')))
        elif os.path.exists('%s/%s_Posterior_Diag.nc.gz' % (dir_longsave,startdate.strftime('%Y%m%d%H'))):
            print "Unzipping Posterior file"
            gzipped_flag = True
            os.system('gunzip %s/%s_Posterior_Diag.nc.gz' % (dir_longsave, startdate.strftime('%Y%m%d%H')))
            os.system('ln -sf %s/%s_Posterior_Diag.nc %s_Posterior_Diag.nc' % (dir_longsave,\
                     startdate.strftime('%Y%m%d%H'), startdate.strftime('%Y%m%d%H'))) 
        else:
            print "!!!!!! Error !!!!!!!"
            print "Could not find file %s_Posterior_Diag.nc in longsave directory."\
                  % startdate.strftime('%Y%m%d%H')
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



def bluefire_submit(mem):
    # Write a wrapper c-shell script to 
    # set the environment variables and run on
    # bluefire
    from WRF_dart_param import NCAR_GAU_ACCOUNT, ADVANCE_TIME_MEMBER, ADVANCE_QUEUE_MEMBER, ADVANCE_CORES_MEMBER, NCAR_ADVANCE_PTILE_MEM
    outfile = open('run_member_%02d.csh' % mem, 'w')
    print >>outfile, "#!/bin/csh"
    print >>outfile, "#==================================================================\\"         
    print >>outfile, "#BSUB -J run_member_%d" % mem                                                   
    print >>outfile, "#BSUB -o run_member_%d.log" % mem                                               
    print >>outfile, "#BSUB -e run_member_%d.err" % mem                                               
    print >>outfile, "#BSUB -P %(NCAR_GAU_ACCOUNT)s" % locals()                                      
    print >>outfile, "#BSUB -W %(ADVANCE_TIME_MEMBER)s" % locals()                                          
    print >>outfile, "#BSUB -q %(ADVANCE_QUEUE_MEMBER)s" % locals()                                         
    print >>outfile, "#BSUB -n %(ADVANCE_CORES_MEMBER)s" % locals()                                         
    print >>outfile, "#BSUB -x"                                                                      
    print >>outfile, '#BSUB -R "span[ptile=%(NCAR_ADVANCE_PTILE_MEM)s]"' % locals()                      
    print >>outfile, "#=================================================================="           
    print >>outfile, ""
    print >>outfile, "# Change into the correct directory"
    print >>outfile, "cd %s/m%d" % (dir_members,mem)
    print >>outfile, ""
    print >>outfile, "# Set the start and end time environment variables"
    print >>outfile, "setenv STARTDATE %s" % startdate.strftime('%Y%m%d%H')
    print >>outfile, "setenv ENDDATE %s" % enddate.strftime('%Y%m%d%H')
    print >>outfile, ""
    print >>outfile, "# Run the python run_member script"
    print >>outfile, "%s/m%d/m%d_run_member.py" % (dir_members, mem, mem)
    outfile.close()
    os.system('mv run_member_%02d.csh %s/m%d' % (mem,dir_members,mem))
    os.system('bsub < %s/m%d/run_member_%02d.csh' % (dir_members,mem,mem))


# FIRST--set the environment variables correctly
os.putenv('STARTDATE',startdate.strftime('%Y%m%d%H'))
os.putenv('ENDDATE',enddate.strftime('%Y%m%d%H'))
# Clear out the log files
os.system('rm -f *.log')
os.system('rm -f *.err')

# Loop through each member
for mem in range(1,int(Ne)+1):
#for mem in [46]:
#for mem in [6,16,27,28,39]:
    print "Member %d..." % (mem),
    # Check for restart
    if restart_flag:
        restart_mem(mem)
    # Copy the run_member script as run_member_%d.py in each directory
    os.system('cp run_member.py %s/m%d/m%d_run_member.py' % (dir_members,mem,mem))


    # Submit this member to the queue
    # Check for which system we are on
    if os.uname()[0] == 'AIX' and os.uname()[1].startswith('be'):
        # We're on bluefire
        bluefire_submit(mem)
    elif os.uname()[1].startswith('ys'):
        # We're on Yellowstone.  Use the Bluefire submission.
        bluefire_submit(mem)
    else:
        os.system('qsub -pe ompi %d -q %s -V -o %s/m%d -e %s/m%d -wd %s/m%d %s/m%d/m%d_run_member.py' \
               % (mpi_numprocs_member,queue_members,dir_members,mem,dir_members,mem,dir_members,mem,\
                 dir_members,mem,mem))
    print "Done."
