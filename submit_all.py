#!/usr/bin/env python

import os, sys
from ens_dart_param import Ne,fct_len,dir_members,dir_utils,dir_longsave,queue_members,mpi_numprocs_member
from datetime import datetime, timedelta
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('-d','--d',dest='starttime',type=int,default=0,
                  help='Start time of this cycle')
parser.add_argument('-l','--length',action='store',dest='length',type=int,default=int(fct_len),
                  help='Override fct_len: length of cycle in minutes')
parser.add_argument('-m','--mems',nargs='+',action='store',dest='sub_mems',default=[],
                  help='Only submit lited members')



opts = parser.parse_args()
starttime = opts.starttime
tdelt = opts.length
endtime = starttime + tdelt

sub_mems = opts.sub_mems
if len(sub_mems) == 0:
    # Submit all members
    memlist = range(1,int(Ne)+1)
else:
    memlist = [int(x) for x in sub_mems if (int(x) > 0) and (int(x)<=Ne)]


gzipped_flag = False



def bluefire_submit(mem):
    # Write a wrapper c-shell script to 
    # set the environment variables and run on
    # bluefire
    from ens_dart_param import NCAR_GAU_ACCOUNT, ADVANCE_TIME_MEMBER, ADVANCE_QUEUE_MEMBER, ADVANCE_CORES_MEMBER, NCAR_ADVANCE_PTILE_MEMBER
    outfile = open('run_member_{:02d}.csh'.format(mem), 'w')
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
    print >>outfile, '#BSUB -R "span[ptile=%(NCAR_ADVANCE_PTILE_MEMBER)s]"' % locals()                      
    print >>outfile, "#=================================================================="           
    print >>outfile, ""
    print >>outfile, "# Change into the correct directory"
    print >>outfile, "cd %s/m%d" % (dir_members,mem)
    print >>outfile, ""
    print >>outfile, "# Set the start and end time environment variables"
    print >>outfile, "setenv STARTTIME {:d}".format(starttime)
    print >>outfile, "setenv ENDTIME {:d}".format(endtime)
    print >>outfile, ""
    print >>outfile, "# Run the python run_member script"
    print >>outfile, "{:s}/m{:d}/m{:d}_run_member.py".format(dir_members, mem, mem)
    outfile.close()
    os.system('mv run_member_{:02d}.csh {:s}/m{:d}'.format(mem,dir_members,mem))
    os.system('bsub < {:s}/m{:d}/run_member_{:02d}.csh'.format(dir_members,mem,mem))


# FIRST--set the environment variables correctly
os.putenv('STARTTIME',str(starttime))
os.putenv('ENDTIME',str(endtime))
# Clear out the log files
os.system('rm -f *.log')
os.system('rm -f *.err')

# Loop through each member
#for mem in range(1,int(Ne)+1):
#for mem in [1]:
#for mem in [6,16,27,28,39]:
for mem in memlist:
    print "Member %d..." % (mem),
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
