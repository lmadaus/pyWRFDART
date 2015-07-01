#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, getopt
from datetime import datetime, timedelta
from numpy import arange
sys.path.append('/home/disk/pvort/nobackup/lmadaus/WRF/DOMAINS/july27')
from WRF_dart_param import *

# Set some defaults
calcend = 0
seed_num = 100000
add_nml_line = False

# Read in command line parameters
(opts,args) = getopt.getopt(sys.argv[1:],'i:d:l:e:f2p')
for o,a in opts:
   if o == '-i':
      today = datetime.datetime.today()
      date_start = today.strftime('%Y%m%d') + str(a)
   if o == '-d':
      date_start = a
   if o == '-l':
      calcend = 1
      lengthdt = datetime.timedelta(hours=int(a))
   if o == '-e':
      seed_num = int(a) * 100000
   if o == '-f':
      max_dom = 1
   if o == '-p':
      add_nml_line = True
   if o == '-2':
      # Rewrite this as if domain 2 was domain 1
      max_dom = '1'
      wps_namelist['geogrid']['e_we'] = [wps_namelist['geogrid']['e_we'][1]]
      wps_namelist['geogrid']['e_sn'] = [wps_namelist['geogrid']['e_sn'][1]]
      wrf_namelist['domains']['e_we'] = [wps_namelist['geogrid']['e_we'][1]]
      wrf_namelist['domains']['e_sn'] = [wps_namelist['geogrid']['e_sn'][1]]
     




def reformat_namelists(wps_namelist,wrf_namelist):
    """ Function that makes time-dependent adjustments to the 
    WRF and WPS namelist dictionaries depending on the parameters
    given on the command line """
    # Format times in the correct way
    startdatedt = datetime.strptime(date_start,'%Y%m%d%H')
    if calcend:
       enddatedt = startdatedt + lengthdt
    else:
       enddatedt = datetime.strptime(date_end,'%Y%m%d%H')

    totalrun = enddatedt - startdatedt
    run_days = int(totalrun.days)
    run_hours = int(totalrun.seconds/3600)
    run_minutes = int((totalrun.seconds-(run_hours*3600))/60)

    # Other variables to convert
    int_seconds = dlbc * 60

    # Compute the dfi times
    # Even if DFI is not requested
    bdfi = startdatedt - timedelta(minutes=dfi_bckstop_window)
    fdfi = startdatedt + timedelta(minutes=dfi_fwdstop_window)
    dfi_len = fdfi - bdfi
    dfi_cutoff_seconds = dfi_len.seconds

    #Seed num for skebs
    nens_skebs = str(int(seed_num/100000))+startdatedt.strftime('%m%d%H')


    # Change the WPS dates
    wps_namelist['share']['start_date'] = ['{:%Y-%m-%d_%H:%M:%S}'.format(startdatedt)] * max_dom
    wps_namelist['share']['end_date'] = ['{:%Y-%m-%d_%H:%M:%S}'.format(enddatedt)] * max_dom

    # Add the use_baseparam line if requested
    if add_nml_line:
        wrf_namelist['dynamics']['use_baseparam_fr_nml'] = True

    # Set the stochastic perturbation for skebs
    wrf_namelist['stoch']['nens'] =  nens_skebs

    # Set the wrf namelist times
    wrf_namelist['time_control']['run_days']         = [run_days] * max_dom
    wrf_namelist['time_control']['run_hours']        = [run_hours] * max_dom
    wrf_namelist['time_control']['run_minutes']      = [run_minutes] * max_dom
    wrf_namelist['time_control']['run_seconds']      =   [0] * max_dom
    wrf_namelist['time_control']['start_year']       =  [startdatedt.year] * max_dom
    wrf_namelist['time_control']['start_month']      =  [startdatedt.month] * max_dom
    wrf_namelist['time_control']['start_day']        =  [startdatedt.day] * max_dom
    wrf_namelist['time_control']['start_hour']       =  [startdatedt.hour] * max_dom
    wrf_namelist['time_control']['start_minute']     =  [startdatedt.minute] * max_dom
    wrf_namelist['time_control']['start_second']     =  [startdatedt.second] * max_dom
    wrf_namelist['time_control']['end_year']         =  [enddatedt.year] * max_dom
    wrf_namelist['time_control']['end_month']        =  [enddatedt.month] * max_dom
    wrf_namelist['time_control']['end_day']          =  [enddatedt.day] * max_dom
    wrf_namelist['time_control']['end_hour']         =  [enddatedt.hour] * max_dom
    wrf_namelist['time_control']['end_minute']       =  [enddatedt.minute] * max_dom
    wrf_namelist['time_control']['end_second']       =  [enddatedt.second] * max_dom
    wrf_namelist['time_control']['interval_seconds'] =  int_seconds 
    # Set the wrf DFI times
    """
    wrf_namelist['dfi_control']['dfi_cutoff_seconds'] = dfi_cutoff_seconds,
    wrf_namelist['dfi_control']['dfi_bckstop_year']   = [bdfi.year] * max_dom,
    wrf_namelist['dfi_control']['dfi_bckstop_month']  = [bdfi.month] * max_dom,
    wrf_namelist['dfi_control']['dfi_bckstop_day']    = [bdfi.day] * max_dom,
    wrf_namelist['dfi_control']['dfi_bckstop_hour']   = [bdfi.hour] * max_dom,
    wrf_namelist['dfi_control']['dfi_bckstop_minute'] = [bdfi.minute] * max_dom,
    wrf_namelist['dfi_control']['dfi_bckstop_second'] = [bdfi.second] * max_dom,
    wrf_namelist['dfi_control']['dfi_fwdstop_year']   = [fdfi.year] * max_dom,
    wrf_namelist['dfi_control']['dfi_fwdstop_month']  = [fdfi.month] * max_dom,
    wrf_namelist['dfi_control']['dfi_fwdstop_day']    = [fdfi.day] * max_dom,
    wrf_namelist['dfi_control']['dfi_fwdstop_hour']   = [fdfi.hour] * max_dom,
    wrf_namelist['dfi_control']['dfi_fwdstop_minute'] = [fdfi.minute] * max_dom,
    wrf_namelist['dfi_control']['dfi_fwdstop_second'] = [fdfi.second] * max_dom,
    """
    # And WRF-VAR seeds and starts for perturbations
    wrf_namelist['wrfvar11']['seed_array1']     =  int('{:%Y%m%d%H}'.format(startdatedt)),
    wrf_namelist['wrfvar11']['seed_array2']     = seed_num,
    wrf_namelist['wrfvar18']['analysis_date']   = '{:%Y-%m-%d_%H:%M:%S}'.format(startdatedt),
    return wps_namelist, wrf_namelist

# This function formats single values
def format_single(v):
    if isinstance(v, int):
        # Boolean is a subclass of int
        if isinstance(v, bool):
            if v:
                return ".true."
            else:
                return ".false."
        else:
            # Format as integer
            return "{:d}".format(v)
    elif isinstance(v, float):
        # Truncate at 4 decimals
        return "{:5.4f}".format(v)
    else:
        # Return the value with single quotes
        return "'{:s}'".format(v)

def var_format(invar):
    """ Decide how to format this namelist value depending on what it is
    """
    # If it's a list, format each item individually and
    # join with commas
    if isinstance(invar, list):
        if len(invar) == 1:
            return format_single(invar[0])
        values = [format_single(v) for v in invar]
        return ', '.join(values[:max_dom])
    else:
        # Return the single format
        return format_single(invar)

def write_namelist(nmld, outfname):
    """ Formats a namelist from a dictionary and writes it to file given by
    outfname
    """
    # Open the output file for writing
    outfile = open(outfname, 'wb')
    # Set the newline character and a default indentation
    nl = '\n'
    indent = '   '

    # Loop through all the sections of the namelist
    for section_name, section in nmld.items():
        # Write the header
        outfile.write(''.join(('&',section_name,nl)))
        # Now loop through all variables and write
        for varname, value in section.items():
            outfile.write(''.join((indent,varname.ljust(24), '= ', var_format(value), ',', nl)))
        # Close the section
        outfile.write(''.join(('/',nl)))
        outfile.write(nl)

    # Close the file
    outfile.close()




if __name__ == '__main__':
    # Update the namelist dictionaries with new parameters
    adjust_wps_namelist, adjust_wrf_namelist = reformat_namelists(wps_namelist, wrf_namelist)
    print("Writing namelists...")
    print(wps_namelist['share']['start_date'][0]," to ", wps_namelist['share']['end_date'][0])
    # Write each of these
    write_namelist(adjust_wps_namelist, 'namelist.wps')
    write_namelist(adjust_wrf_namelist, 'namelist.input')






