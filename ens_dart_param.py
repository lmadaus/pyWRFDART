#! /usr/bin/python

import os
#*******************************************************************************
#  
#  Use this file to modify all the parameters needed to run an experiment.  
#  Each experiment should have one of these files
#
#*******************************************************************************

exp_name='pyCM1DART'                # Name of the experiment
 
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

dir            = '/glade/u/home/lmadaus' # MAIN DIRECTORY
dir_dom    = dir + '/DOMAINS/' + exp_name            # Directory where main experiment is run 
dir_longsave   = dir_dom + '/longsave'               # storage directory 
dir_obs        = dir_dom + '/obs'                    # repository of obs.
dir_utils      = dir + '/UTILS/bin'                      # location of utility codes
dir_members    = dir_dom + '/mems'                   # members directory
dir_assim      = dir_dom + '/assimilation'           # Directory where DART will be run

# Directories where CM1 and DART are found
dir_src_model      = dir + '/cm1/cm1r18/run'        # Where the model executable is located
dir_src_dart      = dir + '/DART/models/wrf/work'   # Where DART executables are located

#**************************************************************
#
#  The following set of parameters change depending on the 
#  cluster that the experiment is being performed on.
#
#  Currently, the scripts will automatically determine the
#  But features like the queue name, the mpi_run_command, and the
#  number of processors come from here
#**************************************************************  

cluster_name        = 'yellowstone'                   # name of cluster nodes
mpi_run_command     = 'mpirun.lsf' # Command to run MPI
queue_members       = 'economy'                    # Queue to run members in
queue_filter        = 'economy'                    # Queue to run filter in          
mpi_numprocs_member = 16                       # Number of processors for member
mpi_numprocs_filter = 128                       # Number of processors for filter
mpi_numprocs_flag   = ''
#mpi_numprocs_flag   = '-np %d' % mpi_numprocs_member      # Flag for numprocs in code
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
dt                 = 5.           # model time step (in sec)
grid_resolutions   = 1000.
out_int            = 20            # Interval to write out files (in MINUTES)
cycle_len          = 60.           # Interval to write restart files (in MINUTES)
                                     # This is also the cycling frequency
fcst_len           = -1          # If 0 --> no forecasts will be produced beyone the cycling interval
                                 # If <0 --> At each cycle, after dropping a restart file at cycle_len
                                 #           minutes, the forecast will continue until reaching the
                                 #           equivalent of exp_len minutes
                                 # Else --> As above, but each forecast will be for an additional
                                 # fcst_len minutes beyond cycle_len

exp_length        = 18*60        # TOTAL length of simulation (in MINUTES)
Ne                = 10           # Number of ensemble members
N_assim           = 500          # Number of assimilations (spinup=0)

assim_start       = 1          	 # Cycle to start assimilation
assim_interval    = 1       	 # how often to assimilate obs (1=every cycle)
inflate_start     = 3            # Which assimlation step to start inflation

## The only inflation option supported is assim_infl_meth=1
## To adjust the inflation mean and std, see the DART namelist
## paramters section below.
assim_infl_meth   = 1          	 # inflation method (see README file for options)
#assim_infl='1.00'            	 # inflation factor (see README file for options)
#assim_inf2='0.25'            	 # inflation factor II (see README file for options)
#assim_bmeth='3'         	 # boundary method (see README file for options) (NOT ENABLED)



#**************************************************************
#
#   Additional EnKF namelist parameters
#   -- These will (in the future) control which variables are a 
#      part of the DART state vector
#
#**************************************************************
# These features are not enabled, but would be easy to do this in the future
update_tsk     = False        # update skin temperature as state variable (NOT ENABLED)
update_lsm     = False        # update land surface state variables (NOT ENABLED)
qual_cont_std  = 3.0  	      # quality control standard deviation (NOT ENABLED)
  
#**************************************************************
# 
#   DART namelist parameters
#   THESE ARE IMPORTANT!!!
#**************************************************************
# filter_nml
async                = 2         # run the model in serial (async=2) or parallel (async=4)
                                 # DO NOT CHANGE async!  This system runs WRF separate from
                                 # DART.
input_qc_threshold   = 3.0       # quality control threshold for obs
outlier_threshold    = 3.0       # outlier threshold (in stds)
debug_filter         = '.true.'	 # trace_execution -->put lots of extra info in the log file
infl_flavor_prior    = 2         # Type of prior inflation (2=spatially-varying state-space infl)   
infl_mean_init_prior = 1.00      # Initial mean prior inflation (usually 1.0 or slightly higher)
infl_sd_init_prior   = 0.6       # Initial std of prior inflation (0.6 is standard)
infl_damp_prior      = 0.9       # Dampening of prior inflation (0.0 to 1.0) --> how quickly
                                 # the inflation can change from cycle to cycle (0.0 shuts off)
infl_lb_prior        = 1.0       # Lower bound of prior inflation mean (will not go below this)
infl_ub_prior        = 1000000.0 # Upper bound of prior inflation mean (will not go above this)
infl_sd_lb_prior     = 0.6       # Lower bound of prior inflation std
infl_flavor_post     = 0         # Type of posterior inflation (0 disables)
infl_mean_init_post  = 1.00      # Initial mean posterior inflation (usually 1.0 or slightly higher) 
infl_sd_init_post    = 0.5       # Initial std of posterior inflation
infl_damp_post       = 1.0       # Dampening of posterior inflation
infl_lb_post         = 1.0       # Lower bound of posterior inflation mean
infl_ub_post         = 1000000.0 # Upper bound of posterior inflation mean
infl_sd_lb_post      = 0.1       # Lower bound of posterior inflation std

# assim_tools_nml section
# Horizontal localization and assimilation methods are here
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
assim_locrad            = 1000.  # localization radius (km)
# cov_cutoff is the actual value put in as the localization radius. DART
# needs the radius specified in global radians.  Do not edit this
# calculation here -- edit the assim_locrad above
cov_cutoff              = str(float(assim_locrad) / 2.0 / 6370.0)[0:5]


# model_nml
# Describes the model state vector
minimal_vars        = True    # ONLY update PH aloft (focus on surface, NOCYCLING ONLY!)
update_surf         = True	  # update U10, V10, T2, TH2, Q2 as the state variable
update_psfc         = True	  # update PSFC as the state variable
num_moist_vars      = 6	          # Depending on the microphysics scheme used,
                                  # the number of moisture variables present
                                  # in the wrfout file [3,5,6]
vert_loc_coord      = 3	          # Gives the vertical coordinate to use when making
                                  # computations on the WRF state and doing vertical
                                  # localization
                                  # 1 = model level, 2 = pressure, 3 = height. Default is 3



# location_nml
# This deals with vertical localization
horizontal_localization_only    =  '.false.'   # True to not use vertical localization

# You only need to set the values for the assim_locrad_vert_****  parameter that
# corresponds to the vert_loc_coord you specified just above.  Once again, DART wants
# these coordinates in global radians, so don't change the calculation parameters

assim_locrad_vert_lev='5.0'     		# vertical localization radius (mod lev)
vert_norm_lev='%.3f' % (float(assim_locrad_vert_lev)/2.0/float(cov_cutoff))

assim_locrad_vert_pres='700.0'   		# vertical localization radius (hPa)
vert_norm_pres='%.3f' % (float(assim_locrad_vert_pres)*100.0/2.0 /float(cov_cutoff))

assim_locrad_vert_hght='16000.0'         	# vertical localization radius (m)
vert_norm_hght=str(int(float(assim_locrad_vert_hght) / 2.0 / float(cov_cutoff)))


#### THESE SECTIONS OF THE DART NAMELIST HAVE TO DO WITH OBSERVATION PROCESSING ####

# obs_sequence_tool
window_minutes    = 15             # Assimilation window--assimilate obs only within
                                   # this many MINUTES of the assimilation time 

# obs_diag_nml
save_obs_locs     = '.true.'       # create observaion_locations.xxx.dat for each bin
obs_diag_verbose  = '.true.'       # verbosity in obs_diag


# wrf_obs_preproc_nml
obs_bdy              =  5.0        # no. of gridpoints near the boundary to remove obs
inc_bdy_obs_err      = '.false.'   # increase obs. error of near-boundary obs
max_obs_fac          =  2.5        # max. factor by which to increase obs. error
obs_dist_bdy         =  10.0       # linearly ramp up the obs. error up to this grid point
thin_aircraft        = '.true.'    # True to thin aircraft to grid spacing 
thin_sat_winds       = '.true.'    # True to thin satellite derived winds
aircraft_horiz_dist  =  45.0       # super-ob horizontal distance for acars/aircraft (km)
aircraft_vert_dist   =  2500.0     # super-ob pressure interval for acars/aircraft (Pa)
sat_wind_horiz_dist  =  45.0       # super-ob horizontal distance for ctw/scat (km)
sat_wind_vert_dist   =  2500.0     # super-ob pressure interval for acars/aircraft (Pa)
soundings_sig_wndtmp = '.false.'   # include sig. levels in rawindsonde
sfc_elev_check       = '.true.'    # sfc elevation check
sfc_elev_tol         =  10000.0      # sfc elevation difference tolerance (meters)
obs_pres_top         =  10000.0    # Exclude obs above this pressure level (Pa)
obs_hght_top         =  20000.0    # Exclude obs above this height level (meters)
thin_sat_winds       = '.true.'    # True to thin SAT winds to grid spacing
thin_acars           = '.true.'    # True to thin ACARS to grid spacing
thin_ctw             = '.true.'    # True to thin cloud winds to 1 degree grid
thin_scat            = '.true.'    # True to thin scatterometer winds to grid spacing
thin_aircraft        = '.true.'    # True to thin Aircraft to grid spacing



# interpol_diag_nml
# This section describes what output is present in the interpolated
# files (if you enabled flag_make_interpol above) 
interpol_diag_include_slp                 = '.true.'
interpol_diag_include_wind_components     = '.true.'
interpol_diag_include_height_on_pres      = '.true.'
interpol_diag_include_temperature         = '.true.'
interpol_diag_include_rel_humidity        = '.true.'
interpol_diag_include_surface_fields      = '.true.'
interpol_diag_include_sat_ir_temp         = '.true.'
# Newer versions of interpol_diag utility want the pressure levels
# in hPa, but older vesions want them in Pa.  You'll have to experiment.
#interpol_diag_pres_levels='925., 850., 700., 500., 300., 250., 200., 150., 100., 50.'
interpol_diag_pres_levels='92500., 85000., 70000., 50000., 30000., 25000., 20000.'
interpol_diag_height_levels='10., 50., 100., 500., 1000.'

# covariance_relax_nml
# NOT ENABLED at this time
use_cov_relax         = False   # (NOT ENABLED)
cov_relax_prior_scale = 0.75    # (NOT ENABLED)

# inflate_ens_nml
# External inflation-- NOT ENABLED at this time
init_assim_scale      = 1.50           # factor to scale initial ensemble
inflate_lsm           = '.false.'      # inflate LSM during init_assim
inflate_tsk           = '.true.'      # inflate TSK during init_assim

#**************************************************************
# 
#   The parameters in this section deal with how to produce 
#   ensemble forecasts that can be integrated for arbitrary time
#
#   AUTOMATION OF THIS IS NOT YET ENABLED IN THE SYSTEM
#**************************************************************

flag_run_ens_fcst='ye'     # yes if want ensemble forecasts
fcst_proc_start='-1'     # processor to start running ens. forecasts (-1 none)
fcst_proc_end='-1'       # processor to end running ens. forecasts (-1 none)

run_free_wrf='no'            # yes if want to run WRF from global model conditions
flag_keep_raw_fcst='no'      # yes to keep raw WRF output (no thinned files)
ens_fcst_bdy_meth='3'    # boundary perturbation method (1=GEV, 2=CTS, 3=FCP)
ens_fcst_len='48'        # number of hours an ensemble forecast lasts
ens_fcst_int='12'        # number of hours between ensemble forecast starts
ens_fcst_start='2'       # assimilation time to start ensemble forecasts
ens_fcst_end='20'        # assimilation time to end ensemble forecasts
ens_fcst_out_int='6'    # number of hours between WRF model output in forecast files
ens_fcst_first_out='12'  # first time to place in thinned output files
ens_fcst_fcp_final='2.4'     # scaling factor for last time boundary condition
ens_fcst_parent_dgbc='180'	# ensemble forecast parent BCs interval (min)
flag_make_opcen_ens_fcst='no'		# yes to interpolate grib files to grid
ens_bc_pscale='0.65'         # perturbation scale for WRF-VAR 48 hr forecast






# A few additional calcualtions
# Calculate the number of integration time steps
time_step=str((float(cycle_len)*60)/float(dt))
assim_len=str(float(N_assim)*float(cycle_len))
cycle_len_hrs=str(float(cycle_len)/60)

