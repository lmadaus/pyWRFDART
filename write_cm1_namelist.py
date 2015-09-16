#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, getopt
from datetime import datetime, timedelta
from numpy import arange
from collections import OrderedDict
from netCDF4 import Dataset
sys.path.append('/glade/p/work/lmadaus/cm1/pyCM1DART')
from ens_dart_param import *

# Set some defaults
calcend = 0
seed_num = 100000
add_nml_line = False
irst = 0
rstnum = 1

# Set up default total_len
if fcst_len < 0:
    total_len = exp_length
else:
    total_len = cycle_len + fcst_len

# Read in command line parameters
(opts,args) = getopt.getopt(sys.argv[1:],'r:l:')
for o,a in opts:
    if o == '-l':
        # Overwrite total_len from ens_dart_param
        total_len = int(a)*60
        #exp_length = int(a)*60
    elif o == '-r':
        # Change the restart number to match this file name
        irst = 1
        if a.startswith('cm1out_rst'):
            # Strip the integer from the filename
            rstnum = int(a.split('_')[-1][:-3])
        else:
            # Just read as an integer
            rstnum = int(a)
# Figure out what to put the restart time as
if irst == 1:
    try:
        rstfile = Dataset('cm1out_rst_{:06d}.nc'.format(rstnum),'r')
        rsttime = rstfile.variables['time'][0]
        # New timax is added to this
        #exp_length = rsttime + cycle_len
    except:
        print("Unable to access file cm1out_rst_{:06d}.nc to find restart time!".format(rstnum))
        pass


     

def set_namelist_defaults():
    namelist = OrderedDict()

    namelist['param0'] = {
        'nx'           :  104,
        'ny'           :  104,
        'nz'           :  78,
        'nodex'        : 4,
        'nodey'        : 4,
        'ppnode'       : mpi_numprocs_member,
        'timeformat'   : 1,
        'timestats'    : 1,
        'terrain_flag' : False,
        'procfiles'    : False,
    }

    
    namelist['param1'] = {
        'dx'       : grid_resolutions,
        'dy'       : grid_resolutions,
        'dz'       : 290.0,
        'dtl'      : dt,
        'timax'    : float(exp_length)*60,
        'run_time' : float(total_len),
        'tapfrq'   : 300.0,
        'rstfrq'   : float(cycle_len),
        'statfrq'  : 60.*15,
        'prclfrq'  : 60.,
    }

    namelist['param2'] = {
     'adapt_dt'  :  1,
     'irst'      :  irst,
     'rstnum'    :  rstnum,
     'iconly'    :  0,
     'hadvordrs' :  5,
     'vadvordrs' :  5,
     'hadvordrv' :  5,
     'vadvordrv' :  5,
     'pdscheme'  :  1,
     'apmasscon' :  1,
     'advwenos'  :  2,
     'advwenov'  :  0,
     'idiff'     :  0,
     'mdiff'     :  0,
     'difforder' :  6,
     'imoist'    :  1,
     'iturb'     :  0,
     'tconfig'   :  2,
     'bcturbs'   :  1,
     'dns'       :  0,
     'irdamp'    :  1,
     'hrdamp'    :  0,
     'psolver'   :  3,
     'nsound'    :  6,
     'ptype'     :  2,
     'ihail'     :  0,
     'iautoc'    :  0,
     'icor'      :  0,
     'pertcor'   :  0,
     'eqtset'    :  2,
     'idiss'     :  1,
     'efall'     :  0,
     'rterm'     :  0,
     'wbc'       :  1,
     'ebc'       :  1,
     'sbc'       :  1,
     'nbc'       :  1,
     'bbc'       :  3,
     'tbc'       :  1,
     'irbc'      :  4,
     'roflux'    :  0,
     'isnd'      :  7,
     'iwnd'      :  1,
     'itern'     :  0,
     'iinit'     :  0,
     'irandp'    :  1,
     'ibalance'  :  0,
     'iorigin'   :  1,
     'axisymm'   :  0,
     'imove'     :  0,
     'iptra'     :  0,
     'npt'       :  1,
     'pdtra'     :  0,
     'iprcl'     :  0,
     'nparcels'  :  1,
    }

    namelist['param3'] = {
     'kdiff2'  :   75.0,
     'kdiff6'  :   0.040,
     'fcor'    : 0.00005,
     'kdiv'    : 0.10,
     'alph'    : 0.60,
     'rdalpha' : 3.3333333333e-3,
     'zd'      : 15000.0,
     'xhd'     : 100000.0,
     'umove'   : 12.5,
     'vmove'   :  3.0,
     'v_t'     :      7.0,
     'l_h'     :   1000.0,
     'l_inf'   :    100.0,
     'ndcnst'  :    250.0,
    }

    namelist['param11'] = {
     'radopt'  :        1,
     'dtrad'   :    300.0,
     'ctrlat'  :    33.1789,
     'ctrlon'  :   -86.7822,
     'year'    :     2014,
     'month'   :       8,
     'day'     :       8,
     'hour'    :       12,
     'minute'  :       00,
     'second'  :       00,
    }

    namelist['param12'] = {
     'isfcflx'    :      1,
     'sfcmodel'   :      2,
     'oceanmodel' :      1,
     'ipbl'       :      1,
     'initsfc'    :      1,
     'tsk0'       :  296.2,
     'tmn0'       :  296.2,
     'xland0'     :    1.0,
     'lu0'        :     14,
     'season'     :      1,
     'cecd'       :      3,
     'pertflx'    :      0,
     'cnstce'     :  0.001,
     'cnstcd'     :  0.001,
     'isftcflx'   :      0,
     'iz0tlnd'    :      0,
     'oml_hml0'   :   50.0,
     'oml_gamma'  :   0.14,
    }

    namelist['param4'] = {
     'stretch_x' :      0,
     'dx_inner'  :    1000.0,
     'dx_outer'  :    7000.0,
     'nos_x_len' :   40000.0,
     'tot_x_len' :  120000.0,
    }

    namelist['param5'] = {
     'stretch_y' :      0,
     'dy_inner'  :    1000.0,
     'dy_outer'  :    7000.0,
     'nos_y_len' :   40000.0,
     'tot_y_len' :  120000.0,
    }

    namelist['param6'] = {
     'stretch_z' :  1,
     'ztop'      :  18000.0,
     'str_bot'   :   3200.0,
     'str_top'   :   9000.0,
     'dz_bot'    :   80.0,
     'dz_top'    :   500.0,
    }

    namelist['param7'] = {
     'bc_temp'   : 1,
     'ptc_top'   : 250.0,
     'ptc_bot'   : 300.0,
     'viscosity' : 25.0,
     'pr_num'    : 0.72,
    }

    namelist['param8'] = {
     'var1'      :   0.0,
     'var2'      :   0.0,
     'var3'      :   0.0,
     'var4'      :   0.0,
     'var5'      :   0.0,
     'var6'      :   0.0,
     'var7'      :   0.0,
     'var8'      :   0.0,
     'var9'      :   0.0,
     'var10'     :   0.0,
    }

    namelist['param9'] = {
     'output_path'      : './',
     'output_basename'  : 'cm1out',
     'output_format'    : 2,
     'output_filetype'  : 1,
     'output_interp'    : 0,
     'restart_format'   : 2,
     'restart_filetype' : 1,
     'output_rain'      : 1,
     'output_sws'       : 0,
     'output_svs'       : 0,
     'output_sps'       : 0,
     'output_srs'       : 0,
     'output_sgs'       : 0,
     'output_sus'       : 1,
     'output_shs'       : 1,
     'output_coldpool'  : 0,
     'output_sfcflx'    : 1,
     'output_sfcparams' : 1,
     'output_sfcdiags'  : 1,
     'output_zs'        : 0,
     'output_zh'        : 0,
     'output_basestate' : 1,
     'output_th'        : 1,
     'output_thpert'    : 1,
     'output_prs'       : 1,
     'output_prspert'   : 0,
     'output_pi'        : 0,
     'output_pipert'    : 1,
     'output_rho'       : 0,
     'output_rhopert'   : 0,
     'output_tke'       : 1,
     'output_km'        : 1,
     'output_kh'        : 1,
     'output_qv'        : 1,
     'output_qvpert'    : 1,
     'output_q'         : 1,
     'output_dbz'       : 1,
     'output_buoyancy'  : 1,
     'output_u'         : 0,
     'output_upert'     : 1,
     'output_uinterp'   : 1,
     'output_v'         : 0,
     'output_vpert'     : 1,
     'output_vinterp'   : 1,
     'output_w'         : 1,
     'output_winterp'   : 1,
     'output_vort'      : 0,
     'output_pv'        : 0,
     'output_uh'        : 0,
     'output_pblten'    : 0,
     'output_dissten'   : 0,
     'output_dissheat'  : 0,
     'output_mptend'    : 0,
     'output_fallvel'   : 0,
     'output_nm'        : 0,
     'output_def'       : 1,
     'output_turbten'   : 0,
     'output_impdiften' : 0,
     'output_radten'    : 1,
     'output_psfc'      : 1,
    }
    namelist['param10'] = {
     'stat_w'        : 1,
     'stat_u'        : 1,
     'stat_v'        : 1,
     'stat_rmw'      : 1,
     'stat_pipert'   : 1,
     'stat_prspert'  : 1,
     'stat_thpert'   : 1,
     'stat_q'        : 1,
     'stat_tke'      : 1,
     'stat_km'       : 1,
     'stat_kh'       : 1,
     'stat_div'      : 1,
     'stat_rh'       : 1,
     'stat_rhi'      : 1,
     'stat_the'      : 1,
     'stat_cloud'    : 1,
     'stat_sfcprs'   : 1,
     'stat_wsp'      : 1,
     'stat_cfl'      : 1,
     'stat_vort'     : 1,
     'stat_tmass'    : 1,
     'stat_tmois'    : 1,
     'stat_qmass'    : 1,
     'stat_tenerg'   : 1,
     'stat_mo'       : 1,
     'stat_tmf'      : 1,
     'stat_pcn'      : 1,
     'stat_qsrc'     : 1,
    }

    namelist['param13'] = {
     'prcl_th'       : 1,
     'prcl_t'        : 1,
     'prcl_prs'      : 1,
     'prcl_ptra'     : 1,
     'prcl_q'        : 1,
     'prcl_nc'       : 1,
     'prcl_km'       : 1,
     'prcl_kh'       : 1,
     'prcl_tke'      : 1,
     'prcl_dbz'      : 1,
     'prcl_b'        : 1,
     'prcl_vpg'      : 1,
     'prcl_vort'     : 1,
     'prcl_rho'      : 1,
     'prcl_qsat'     : 1,
     'prcl_sfc'      : 1,
    }

    namelist['nssl2mom_params'] = {
       'alphah'  : 0,     # shape parameter of graupel
       'alphahl' : 0.5,   # shape parameter of hail
       'ccn'     : 2.0e9,  # base ccn concentration
       'cnor'    : 8.e6,  # for single moment only
       'cnoh'    : 4.e4,  # for single moment only
    }


    return namelist





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

def write_namelist(nmld, outfname='namelist.input'):
    """ Formats a namelist from a dictionary and writes it to file given by
    outfname
    """
    # Open th output file for writing
    # Set the newline character and a default indentation
    nl = '\n'
    indent = ' '
    with open(outfname, 'wb') as outfile:
        print('', file=outfile)
        # Loop through all the sections of the namelist
        #section_names = nmld.keys()
        #section_names.sort()
        #section_names.reverse()
        for section_name, section in nmld.items():#section_names:
            #section = nmld[section_name]
            # Write the header
            #outfile.write(''.join((' &',section_name,nl)))
            print(''.join((' &',section_name)), file=outfile)
            # Now loop through all variables and write
            for varname, value in section.items():
                print(''.join((indent,varname.ljust(15), '= ',
                               var_format(value), ',')), file=outfile)
            # Close the section
            print(' /', file=outfile)
            print('', file=outfile)
        outfile.close()





if __name__ == '__main__':
    # Update the namelist dictionaries with new parameters
    print("Writing namelist...")
    # Write based on the defaults
    namelist_values = set_namelist_defaults()
    write_namelist(namelist_values, 'namelist.input')






