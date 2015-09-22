#!/usr/bin/env python

# set_all_scripts.py
# IMPORTANT TO RUN THIS SCRIPT BEORE STARTING
# This script ensures you have everything you need for the
# ensemble to run.  It also goes into all control scripts and
# changes the path to WRF_dart_param.py to the correct path.  This
# script must be run from the main dir_wrf_dom directory

import os
from ens_dart_param import dir_longsave, dir_obs, dir_assim


req_files = ['run_member.py','submit_all.py','autorun_ensemble.py','write_cm1_namelist.py',\
             'make_namelist_dart.py','submit_filter.py','make_ensemble.py',\
             'ens_dart_param.py','ens_dart_obtypes.py','make_osse_obs.py']

def main():
    check_files()
    change_sys_path()
    os.system('chmod +x *.py')
    make_sub_dirs()


def make_sub_dirs():
    """ Function to go through the directories we expect
    to find and create them if they don't exist """
    check_dirs = [dir_longsave, dir_obs, dir_assim]
    for d in check_dirs:
        if not os.path.exists(d):
            os.system('mkdir {:s}'.format(d))

def check_files():
    # Loop through the required files and 
    # make sure they exists
    curpath = os.getcwd()

    for file in req_files:
        filefound = does_it_exist(curpath+'/'+file)
        if not filefound:
            exit(1)

    print "Success! All required files found."

def change_sys_path():
    # Function that goes into each script that contains
    # a 'sys.path.append' line and replaces it with one that
    # reflects the current working directory structure
    print "Re-setting wrfdart directory in scripts."
    new_wrfdart_path = os.getcwd()
    for file in req_files:
        print file
        oldfile = open(file, 'r')
        newfile = open('new.temp','w')
        for line in oldfile:
            if not line.startswith('sys.path.append'):
                print >>newfile,line[0:-1]  # Skip the newline character at the end
            else:
                print >>newfile,"sys.path.append('%s')" % new_wrfdart_path
        newfile.close()
        oldfile.close()
        os.system('mv new.temp %s' % file)



def does_it_exist(filepath):
    # Function to determine if a certain file exists.
    # If not, return an error.
    if os.path.exists(filepath):
        print "Found file %s" % filepath
        return True
    else:
        print "Error--could not find file %s" % filepath
        return False 





if __name__ == '__main__':
    main()

