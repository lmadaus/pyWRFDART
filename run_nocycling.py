#!/usr/bin/env python
from __future__ import print_function, division
from datetime import datetime, timedelta
from netCDF4 import Dataset
from numpy import subtract, float32
import os

nMems = 50
archive_dir = '/home/disk/jabba/lmadaus/nobackup/july27/july27_eval_only'
#archive_dir = '/home/disk/trof/lmadaus/nobackup/july27_standard_only'
fdate = datetime(2014,7,26,12)

def run_cycle(fdate=fdate, archive=archive_dir):
    # Repopulate the ensemble
    #populated = repopulate_ensemble(fdate,archive,posterior=False)
    #if not populated:
    #    print("Re-population of ensemble members failed!")
    #    return None
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




def repopulate_ensemble(fdate, fdir, posterior=False, compute_tendency=True):
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

    os.system('cp old_input.nml input.nml')
    # Now attempt to run diag2wrfinput
    for mem in range(nMems):
        print("Member {:d}".format(mem+1))
        os.system('unlink wrfinput_d01')
        os.system('ln -sf mems/m{:d}/wrfinput_d01 wrfinput_d01'.format(mem+1))
        os.system('./diag2wrfinput {:s} {:s} {:d}'.format(master_file,\
                                                          fdate.strftime('%Y-%m-%d_%H:%M:%S'),\
                                                          mem+1))
        # Now compute tendency
        if compute_tendency:
            print("Computing altimeter tendency...")
            this_dset = Dataset('mems/m{:d}/wrfinput_d01'.format(mem+1),'r+')
            last_dset = Dataset('mems/m{:d}/prev_psfc.nc'.format(mem+1),'r')
            # Get the elevation
            elev = this_dset.variables['HGT'][0]
            this_psfc = this_dset.variables['PSFC'][0]
            last_psfc = last_dset.variables['PSFC'][0,0]
            # Compute altimeters
            this_alt = altimeter(this_psfc, elev)
            last_alt = altimeter(last_psfc, elev)
            # Compute tendency
            tend = subtract(this_alt, last_alt)
            # Make a new variable
            if 'ALT_TEND' not in this_dset.variables.keys():
                this_dset.createVariable('ALT_TEND','f', ('Time','south_north','west_east',))
            tendvar = this_dset.variables['ALT_TEND']
            tendvar.setncattr('FieldType',this_dset.variables['PSFC'].FieldType)
            tendvar.setncattr('MemoryOrder',"XY ")
            tendvar.setncattr('description', "ALTIMETER TENDENCY")
            tendvar.setncattr('units', "hPa")
            tendvar.setncattr('stagger', "") 
            tendvar.setncattr('coordinates', "XLONG XLAT")
            this_dset.variables['ALT_TEND'][0,:,:] = tend.astype(float32)
            this_dset.close()
            last_dset.close()
            os.system('rm -f mems/m{:d}/prev_psfc.nc'.format(mem+1))
            os.system('ncecat -v PSFC mems/m{:d}/wrfinput_d01 mems/m{:d}/prev_psfc.nc'.format(mem+1,mem+1))


        #os.system('cp wrfinput_d01 mems/m{:d}/wrfinput_d01'.format(mem+1))
    return True

def altimeter(psfc, elev):
    return ((psfc/100. - 0.3) ** 0.190284 + 8.4228807E-5 * elev) ** 1/0.190284


if __name__ == '__main__':
    curdate = datetime(2014,7,26,15)
    enddate = datetime(2014,7,26,15)
    while curdate <= enddate:
        run_cycle(curdate)
        curdate += timedelta(hours=1)
