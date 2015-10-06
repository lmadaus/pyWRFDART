#! /usr/bin/python

import os
#*******************************************************************************
#  
#  Use this file to modify all the parameters needed to run an experiment.  
#  Each experiment should have one of these files
#
#*******************************************************************************

exp_name='kdvn_ensemble'                # Name of the experiment
 
#**************************************************************
#
#   This set of flags controls some of the basic execution of 
#   the script for running a case study
#
#**************************************************************

flag_compute_tendency = False      # True to compute the altimeter tendencies for each out file
flag_make_interpol    = False      # True to perform the interpolation
flag_keep_raw         = True       # True to keep original files
flag_compress_diag    = False       # True to gzip the diag files in longsave
flag_keep_outs        = True      # True to preserve output files in each member's directory
flag_precip_diag      = False      # True te extract precip diag files (NOT ENABLED)
flag_keep_diags       = False       # True to keep the Prior and Posterior Diag files

flag_obs_diag         = False      # perform observation diagnostics
flag_inf_diag         = False      # perform inflation diagnostics

flag_direct_netcdf_io = True       # Use experimental DART IO to read/write
                                   # directly from/to restart netcdf files

#**************************************************************
#
#   This section sets the location of where files are to be 
#   created and where necessary files can be found
#
#   If basic directory structure is followed, the only two things
#   that need to be changed for all of these to work are the
#   exp_name parameter at the very top of this file and the 
#   dir parameter below.
#**************************************************************  

dir            = '/home/disk/pvort/nobackup/lmadaus/cm1' # MAIN DIRECTORY
dir_dom    = dir + '/DOMAINS/' + exp_name            # Directory where main experiment is run 
dir_longsave   = dir_dom + '/longsave'               # storage directory 
dir_obs        = dir_dom + '/obs'                    # repository of obs.
dir_utils      = dir + '/UTILS/bin'                      # location of utility codes
dir_members    = dir_dom + '/mems'                   # members directory
dir_assim      = dir_dom + '/assimilation'           # Directory where DART will be run

# Directories where CM1 and DART are found
dir_src_model      = dir + '/r18/cm1r18/run'        # Where the model executable is located
dir_src_dart      = dir + '/DART_CM1/models/cm1/work'   # Where DART executables are located

#**************************************************************
#
#  The following set of parameters change depending on the 
#  cluster that the experiment is being performed on.
#
#  Currently, the scripts will automatically determine the
#  But features like the queue name, the mpi_run_command, and the
#  number of processors come from here
#**************************************************************  

cluster_name        = 'enkf'                   # name of cluster nodes
mpi_run_command     = 'mpirun' # Command to run MPI
queue_members       = 'fast'                    # Queue to run members in
queue_filter        = 'fast'                    # Queue to run filter in          
mpi_numprocs_member = 16                       # Number of processors for member
mpi_numprocs_filter = 64                       # Number of processors for filter
#mpi_numprocs_flag   = ''
mpi_numprocs_flag   = '-np %d' % mpi_numprocs_member      # Flag for numprocs in code
                                                          # for member.  Bluefire does
                                                          # blank for bluefire

# Extra parameters for running on Bluefire
NCAR_GAU_ACCOUNT     = 'UWAS0031'                   # Account to charge to at NCAR
ADVANCE_TIME_FILTER  = '0:45'                # Estimate of time for filter to run
ADVANCE_TIME_MEMBER  = '0:20'                # Estimate of time for a single member to run 
ADVANCE_QUEUE_FILTER = queue_filter          # Name of NCAR queue to use for filter 
ADVANCE_QUEUE_MEMBER = queue_members         # Name of NCAR queue to use for members
ADVANCE_CORES_FILTER = mpi_numprocs_filter  # Number of cores to use (multiples of 32)
ADVANCE_CORES_MEMBER = mpi_numprocs_member  # Number of cores to use (multiples of 32)
NCAR_ADVANCE_PTILE_FILTER   = '16'                  # How many processes per core on Bluefire
NCAR_ADVANCE_PTILE_MEMBER   = '16'                  # How many processes per core on Bluefire



#**************************************************************
#
#   This set of parameters deals with parameters related to the 
#   ensemble itself
#
dt                 = 5.          # model time step (in sec)
grid_resolutions   = 1000.
out_int            = 10*60       # Interval to write out files (in seconds)
cycle_len          = 60*60          # Interval to write restart files (in seconds)
                                 # This is also the cycling frequency
fcst_len           = 300*60        # If 0 --> no forecasts will be produced beyone the cycling interval
                                 # If <0 --> At each cycle, after dropping a restart file at cycle_len
                                 #           minutes, the forecast will continue until reaching the
                                 #           equivalent of exp_len minutes
                                 # Else --> As above, but each forecast will be for an additional
                                 # fcst_len minutes beyond cycle_len

exp_length        = 10*3600        # TOTAL length of simulation (in seconds))
Ne                = 50           # Number of ensemble members

assim_start       = 1          	 # Cycle to start assimilation
assim_interval    = 1       	 # how often to assimilate obs (1=every cycle)
inflate_start     = 2            # Which assimlation step to start inflation



#**************************************************************
#
#   Additional EnKF namelist parameters
#   -- These will (in the future) control which variables are a 
#      part of the DART state vector
#
#**************************************************************
  
# assim_tools_nml section
# Horizontal localization and assimilation methods are here
assim_infl_meth   = 1          	 # inflation method (see README file for options)
assim_meth              = 1      # Type of filter to run 
                                 # (1=EAKF, 2=ENKF, 3=Kernel filter, 
                                 # 4=particle filter. 7=Boxcar kernel filter)
assim_loc_meth          = 1      # select localization method (1=Gaspari-Cohn, 2=boxcar,
                                 # 3=ramped boxcar)
adaptive_loc_threshold  = 1600   # Used to dynamically shrink localization in areas of
                                 # dense observations.  If the number of observations within
                                 # 2* the localization radius is greater than this value, it
                                 # shrinks the localization radius until there are only
                                 # this many observations within that radius.
assim_locrad            = 15.  # localization radius (km)
# cov_cutoff is the actual value put in as the localization radius. DART
# needs the radius specified in global radians.  Do not edit this
# calculation here -- edit the assim_locrad above
#cov_cutoff              = str(float(assim_locrad) / 2.0 / 6370.0)[0:5]
cov_cutoff = str(assim_locrad * 1000.)



# A few additional calcualtions
# Calculate the number of integration time steps
time_step=str((float(cycle_len)*60)/float(dt))
cycle_len_hrs=str(float(cycle_len)/60)

