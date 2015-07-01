#!/usr/bin/env python
from __future__ import print_function, division
from datetime import datetime, timedelta
import os

nMems = 50
#archive_dir = '/home/disk/trof/lmadaus/nobackup/july27_eval_only'
archive_dir = '/home/disk/trof/lmadaus/nobackup/july27_standard_only'
fdate = datetime(2014,7,26,12)

def run_cycle(fdate=fdate, archive=archive_dir):
    # Repopulate the ensemble
    populated = repopulate_ensemble(fdate,archive,posterior=True)
    if not populated:
        print("Re-population of ensemble members failed!")
        return None
    # Do the assimilation
    assimilation = run_assimilation_sequence(fdate)

def run_assimilation_sequence(fdate):
    """ Switch to the wrfdart directory and run the assimilation sequence for
    this time """
    #curdir = os.getcwd()
    #os.chdir('./wrfdart')
    os.system('./submit_filter.py -d {:%Y%m%d%H}'.format(fdate))
    #os.chdir(curdir)
    return True




def repopulate_ensemble(fdate, fdir, posterior=False):
    """ Grab Prior Diag files from an archive and repopulate the ensemble """
    # Look for a Diag file
    if posterior:
        master_file = 'Posterior_Diag.nc'
    else:
        master_file = 'Prior_Diag.nc'
    # Clear out an old file
    if os.path.exists(master_file):
        os.system('rm -f {:s}'.format(master_file))
    # Search for a new file
    search_path = '/'.join([fdir,'_'.join((fdate.strftime('%Y%m%d%H'),master_file))])

    # It's either zipped or not
    if os.path.exists(search_path):
        os.system('ln -sf {:s} {:s}'.format(search_path, master_file))
    elif os.path.exists('.'.join([search_path,'gz'])):
        # Have to unzip the file
        print("Must unzip archive file! Unzipping...")
        os.system('gunzip {:s}.gz'.format(search_path))
        os.system('ln -sf {:s} {:s}'.format(search_path, master_file))
    else:
        # File doesn't exists
        print("Could not find file at:", search_path)
        return False

    # Now attempt to run diag2wrfinput
    for mem in range(nMems):
        print("Member {:d}".format(mem+1))
        os.system('unlink wrfinput_d01')
        os.system('ln -sf mems/m{:d}/wrfinput_d01 wrfinput_d01'.format(mem+1))
        os.system('./diag2wrfinput {:s} {:s} {:d}'.format(master_file,\
                                                          fdate.strftime('%Y-%m-%d_%H:%M:%S'),\
                                                          mem+1))
        #os.system('cp wrfinput_d01 mems/m{:d}/wrfinput_d01'.format(mem+1))
    return True


if __name__ == '__main__':
    curdate = datetime(2014,7,26,12)
    enddate = datetime(2014,7,26,12)
    while curdate <= enddate:
        run_cycle(curdate)
        curdate += timedelta(hours=1)
