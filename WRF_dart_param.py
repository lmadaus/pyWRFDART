#! /usr/bin/python

import os
#*******************************************************************************
#  
#  Use this file to modify all the parameters needed to run an experiment.  
#  Each experiment should have one of these files
#
#*******************************************************************************

exp_name='ncar_ensemble'                # Name of the experiment
 
#**************************************************************
#
#   This set of flags controls some of the basic execution of 
#   the script for running a case study
#
#**************************************************************

flag_compile          = False      # True to compile new filter (NOT ENABLED)
flag_make_domain      = False      # True to make the domain (NOT ENABLED)
flag_create_ens       = False      # True to create a new experiment. (NOT ENABLED)
flag_make_icbc        = False      # True to create initial and boundary conditions (NOT ENABLED)
flag_make_obs_seq     = False      # True to create obs. seq. files (NOT ENABLED)
flag_make_opcen       = False      # True to interpolate grib files to grid (NOT ENABLED)

flag_pert_bcs         = False      # True to run pert_wrf_bc.  If false, will still run
                                   # update_wrf_bc

flag_compute_tendency = False      # True to compute the altimeter tendencies for each wrfout
                                   # and add altimeter tendency to WRF state vector in DART
flag_make_interpol    = False      # True to perform the interpolation
flag_keep_raw         = True       # True to keep original wrf files
flag_compress_diag    = False       # True to gzip the diag files in longsave
flag_keep_wrfouts     = False      # True to preserve wrfout files in each member's directory
flag_precip_diag      = False      # True to extract precip diag files
flag_keep_diags       = False       # True to keep the Prior and Posterior Diag files

flag_obs_diag         = False      # perform observation diagnostics (NOT ENABLED)
flag_inf_diag         = False      # perform inflation diagnostics (NOT ENABLED)

flag_direct_netcdf_io = True       # Use experimental DART IO to read/write
                                   # directly from/to wrfinput_d01 netcdf files


#**************************************************************
#
#   This set of parameters deals with parameters related to the 
#   model itself
#
#**************************************************************
max_dom            = 2         # maximum number of domains
dt                 = 20           # model time step (in sec)
grid_resolutions   = [15000,3000]
wrfout_int        = 60           # Interval to write wrfout files (in MINUTES)
dlbc              = 60        	 # Interval of global boundary conditions (in MINUTES)
#**************************************************************
#
#   This set of parameters deals with the EnKF parameters 
#   itself.  For more information on some of them, see the 
#   original papers or Hamill's review paper.
#
#**************************************************************
date_start        = '2017051612' # Start of the experiment (YYYYMMDDHH)
date_end          = '2017051600' # End of the experiment (YYYYMMDDHH)
Ne                = 30           # Number of ensemble members
fct_len           = 60*12          # Time between 2 assimilation (in MINUTES)
#fct_len           = 10          # Time between 2 assimilation (in MINUTES)
N_assim           = 500          # Number of assimilations (spinup=0)
N_assim_max       = 500          # Maximum number of assimilations 

assim_start       = 1          	 # Cycle time to start assimilation
assim_interval    = 1       	 # how often to assimilate obs (1=every cycle) (NOT ENABLED)
inflate_start     = 100            # Which assimlation step to start inflation

## The only inflation option supported is assim_infl_meth=1
## To adjust the inflation mean and std, see the DART namelist
## paramters section below.
assim_infl_meth   = 1          	 # inflation method (see README file for options)
#assim_infl='1.00'            	 # inflation factor (see README file for options)
#assim_inf2='0.25'            	 # inflation factor II (see README file for options)
#assim_bmeth='3'         	 # boundary method (see README file for options) (NOT ENABLED)

# digital filtering options
use_dfi = 0                             # DFI option to use (0 is none)
dfi_bckstop_window  = 30                # integrate backward in MINUTES
dfi_fwdstop_window  = 15                # integrate forward  in MINUTES



#**************************************************************
#
#   The parameters in this section deal with how the ensemble 
#   is initialized
#
#**************************************************************
# ONLY INITIALIZATION METHOD 2 is supported
# 3 can be used, but changing this value here doesn't change anything
init_method   = 2          # Method to initialize the ensemble 
                           # (1=cts, 2=fcp, 3=parent ensemble)



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
dir_wrf_dom    = dir + '/DOMAINS/' + exp_name            # Directory where main experiment is run 
dir_longsave   = dir_wrf_dom + '/longsave'               # storage directory 
dir_obs        = dir_wrf_dom + '/obs'                    # repository of obs.
dir_utils      = dir + '/UTILS/bin'                      # location of utility codes
dir_members    = dir_wrf_dom + '/mems'                   # members directory
dir_assim      = dir_wrf_dom + '/assimilation'           # Directory where DART will be run

# Directories where WRF, WRFDA and DART are found
WRFVARDIR         = dir + '/WRFDA'
WRFRUNDIR         = dir + '/WRFV3_ncar/run'
dir_src_wps       = dir + '/WPS'                    # WPS location
dir_src_wps_geog  = dir + '/DATA/geog'              # Location of geo data for WPS
dir_src_wrf       = dir + '/WRFV3_ncar/main'             # WRF location
dir_src_wrfvar    = dir + '/WRFDA/var'                  # WRF-VAR location
dir_src_dart      = dir + '/DART_Manhattan/models/wrf/work'   # Where DART executables are located


#**************************************************************
#
#  The following set of parameters change depending on the 
#  cluster that the experiment is being performed on.
#
#  Currently, the scripts will automatically determine the
#  But features like the queue name, the mpi_run_command, and the
#  number of processors come from here
#**************************************************************  

cluster_name        = 'cheyenne'                   # name of cluster nodes
#mpi_run_command     = 'mpirun.lsf' # Command to run MPI
mpi_run_command     = 'mpiexec_mpt' # Command to run MPI
queue_members       = 'regular'                    # Queue to run members in
queue_filter        = 'regular'                    # Queue to run filter in     
numprocs_per_node   = 36                   # Number of processors on each node of the system
mpi_numprocs_member = 512                       # Number of processors for member
mpi_numprocs_filter = 128                       # Number of processors for filter
#mpi_numprocs_flag   = '-np %d' % mpi_numprocs_member      # Flag for numprocs in code
                                                          # for member.  Bluefire does
                                                          # blank for bluefire
mpi_numprocs_flag = ''
# Calculation to figure out how many nodes to request
nnodes_member = mpi_numprocs_member // numprocs_per_node
if mpi_numprocs_member < numprocs_per_node:
    nnodes_member = 1
elif mpi_numprocs_member % numprocs_per_node !=0:
    nnodes_member += 1

nnodes_filter = mpi_numprocs_filter // numprocs_per_node
if mpi_numprocs_filter < numprocs_per_node:
    nnodes_filter = 1
elif mpi_numprocs_filter % numprocs_per_node !=0:
    nnodes_filter += 1

# Extra parameters for running on Bluefire
NCAR_GAU_ACCOUNT     = 'NASP0002'                   # Account to charge to at NCAR
ADVANCE_TIME_FILTER  = '0:45:00'                # Estimate of time for filter to run
ADVANCE_TIME_MEMBER  = '2:00:00'                # Estimate of time for a single member to run 
ADVANCE_QUEUE_FILTER = queue_filter          # Name of NCAR queue to use for filter 
ADVANCE_QUEUE_MEMBER = queue_members         # Name of NCAR queue to use for members
ADVANCE_CORES_FILTER = mpi_numprocs_filter  # Number of cores to use (multiples of 32)
ADVANCE_CORES_MEMBER = mpi_numprocs_member  # Number of cores to use (multiples of 32)
NCAR_ADVANCE_PTILE   = '16'                  # How many processes per core on Bluefire




# A few additional calcualtions
# Calculate the number of integration time steps
time_step=str((float(fct_len)*60)/float(dt))
assim_len=str(float(N_assim)*float(fct_len))
fct_len_hrs=str(float(fct_len)/60)
#model_grid_ratio=str(float(model_gridspx2)/float(model_gridspx1))
#model_timestep_ratio=model_grid_ratio
gmodnam = 'gfs'
gmodnamu=os.popen('echo %s | tr "[:lower:]" "[:upper:]"' % gmodnam).readlines()[0][0:-1]
username=os.popen('whoami').readlines()[0][0:-1]      # user running this experiment


# boundary zone width and other parameters
dlbc_hrs=str(float(dlbc)/60.0)


#**************************************************************
#
#  These parameters all deal with the WRF model, but 
#  are not changed very often.  See the WRF README file for 
#  explanations of what each does.
#
#  Some WRF namelist parameters, particularly related to the 
#  time control options, are set in variables above
#
#  These parameters are almost exclusively used by
#  write_namelists.py
#**************************************************************  

#### WPS NAMELIST ####
wps_namelist = {}
wps_namelist['share'] = {
    'wrf_core'         : 'ARW',
    'max_dom'          : max_dom,
    'interval_seconds' : dlbc*60,
    'io_form_geogrid'  : 2,
    }

wps_namelist['geogrid'] = {
    'parent_id'         :   [1] + list(range(1,max_dom)),
    'parent_grid_ratio' :   [int(grid_resolutions[0] / n) for n in grid_resolutions],
    'i_parent_start'    :   [0,  93],
    'j_parent_start'    :   [0,  50],
    'e_we'              :   [560, 271],
    'e_sn'              :   [420, 226],
    'geog_data_res'     :   ['30s']*max_dom,
    'dx'                :   grid_resolutions[0],
    'dy'                :   grid_resolutions[0],
    'map_proj'          :   'lambert',
    'ref_lat'           :   39.5,
    'ref_lon'           :   -80.25,
    'truelat1'          :   38.5,
    'truelat2'          :   38.5,
    'stand_lon'         :   -92.5,
    'geog_data_path'    : dir_src_wps_geog,
    }

wps_namelist['ungrib'] = {
    'out_format' : 'WPS',
    'prefix'     : 'FILE',
    }

wps_namelist['metgrid'] = {
    'fg_name'         : 'FILE',
    'io_form_metgrid' : 2,
    }

# WRF NAMELIST
wrf_namelist = {}
wrf_namelist['time_control'] = {
    'input_from_file'    :  [True] * max_dom,
    'history_interval'   :  [wrfout_int] * max_dom,
    'frames_per_outfile' :  [1] * max_dom,
    'restart'            : False,
    'restart_interval'   : 5000,
    'io_form_history'    : 2,
    'io_form_restart'    : 2,
    'io_form_input'      : 2,
    'io_form_boundary'   : 2,
    'io_form_auxinput2'  : 2,
    'debug_level'        : 0,
    }

wrf_namelist['domains'] = {
    'time_step'               : dt,
    'time_step_fract_num'     : 0,
    'time_step_fract_den'     : 1,
    'max_dom'                 : max_dom,
    's_we'                    : [1] * max_dom,
    'e_we'                    : wps_namelist['geogrid']['e_we'],
    's_sn'                    : [1] * max_dom,
    'e_sn'                    : wps_namelist['geogrid']['e_sn'],
    's_vert'                  : [1] * max_dom,
    'e_vert'                  : [51] * max_dom,
    'dx'                      : grid_resolutions,
    'dy'                      : grid_resolutions,
    'grid_id'                 : list(range(1,max_dom+1)),
    'parent_id'               : list(range(max_dom)),
    'i_parent_start'          : wps_namelist['geogrid']['i_parent_start'],
    'j_parent_start'          : wps_namelist['geogrid']['j_parent_start'],
    'parent_grid_ratio'       : wps_namelist['geogrid']['parent_grid_ratio'],
    'parent_time_step_ratio'  : wps_namelist['geogrid']['parent_grid_ratio'],
    'feedback'                : 0,
    'smooth_option'           : 0,
    'lagrange_order'          : 1,
    'interp_type'             : 1,
    'hypsometric_opt'         : 2,
    'lowest_lev_from_sfc'     : False,
    'force_sfc_in_vinterp'    : 1,
    'zap_close_levels'        : 500,
    'smooth_cg_topo'          : True,
    'sfcp_to_sfcp'            : True,
    'use_levels_below_ground' : True,
    'adjust_heights'          : True,
    'eta_levels'              : [1.0000, 0.9980, 0.9940, 0.9870, 0.9750, 0.9590, 0.9390, 0.9160, 0.8920, 0.8650, 0.8350, 0.8020, 0.7660, 0.7270, 0.6850, 0.6400, 0.5920, 0.5420, 0.4970, 0.4565, 0.4205, 0.3877, 0.3582, 0.3317, 0.3078, 0.2863, 0.2670, 0.2496, 0.2329, 0.2188, 0.2047, 0.1906, 0.1765, 0.1624, 0.1483, 0.1342, 0.1201, 0.1060, 0.0919, 0.0778, 0.0657, 0.0568, 0.0486, 0.0409, 0.0337, 0.0271, 0.0209, 0.0151, 0.0097, 0.0047, 0.0000],
    'p_top_requested'         : 5000.,
    'use_adaptive_time_step'  : True,
    'step_to_output_time'     : True,
    'target_cfl'              : 1.2,
    'max_step_increase_pct'   : 5,
    'starting_time_step'      : 15,
    'min_time_step'           : 5,
    'max_time_step'           : 30,
    }

wrf_namelist['physics'] = {
    'mp_physics'                          : [8] * max_dom,
    'do_radar_ref'                        : 1,
    'ra_lw_physics'                       : [4] * max_dom,
    'ra_sw_physics'                       : [4] * max_dom,
    'radt'                                : [30] * max_dom,
    'sf_sfclay_physics'                   : [5] * max_dom,
    'sf_surface_physics'                  : [3] * max_dom,
    'CO2TF'                               : 1,
    'bl_pbl_physics'                      : [5] * max_dom,
    'bldt'                                : [0] * max_dom,
    'cu_physics'                          : [0] * max_dom,
    'cudt'                                : [0] * max_dom,
    'grav_settling'                       : 0,
    'isfflx'                              : 1,
    'ifsnow'                              : 1,
    'icloud'                              : 1,
    'surface_input_source'                : 1,
    'num_soil_layers'                     : 9,
    'sf_urban_physics'                    : 0,
    'rdlai2d'                             : True,
    'usemonalb'                           : True,
    'seaice_threshold'                    : 271.4,
    'mp_zero_out'                         : 2,
    'mp_zero_out_thresh'                  : 1.0E-12,
    'maxiens'                             : 1,
    'maxens'                              : 3,
    'maxens2'                             : 3,
    'maxens3'                             : 16,
    'ensdim'                              : 144,
    'num_land_cat'                        : 21,
    'mosaic_lu'                           : 1,
    'mosaic_soil'                         : 1,
    'prec_acc_dt'                         : 60,
    }

wrf_namelist['fdda'] = {}

wrf_namelist['dynamics'] = {
    'rk_ord'                              : 3,
    'diff_6th_opt'                        : 2,
    'diff_6th_factor'                     : 0.25,
    'w_damping'                           : 1,
    'diff_opt'                            : 1,
    'km_opt'                              : 4,
    'damp_opt'                            : 3,
    'zdamp'                               : 5000.,
    'base_temp'                           : 300.,
    'dampcoef'                            : 0.2,
    'khdif'                               : [0] * max_dom,
    'kvdif'                               : [0] * max_dom,
    'smdiv'                               : [0.1] * max_dom,
    'emdiv'                               : [0.01] * max_dom,
    'epssm'                               : [0.1] * max_dom,
    'time_step_sound'                     : [4] * max_dom,
    'non_hydrostatic'                     : [True] * max_dom,
    }

wrf_namelist['stoch'] = {
    'stoch_force_opt'                     : [1] * max_dom,
    'stoch_vertstruc_opt'                 : [0] * max_dom,
    'perturb_bdy'                         : 1,
    'tot_backscat_psi'                    : 1.0E-04,
    'tot_backscat_t'                      : 5.0E-05,
    'ztau_psi'                            : 3600.0,
    'ztau_t'                              : 3600.0,
    'rexponent_psi'                       : -1.83,
    'rexponent_t'                         : -1.83,
    'kminforc'                            : 1,
    'lminforc'                            : 1,
    'kminforct'                           : 4,
    'lminforct'                           : 4,
    }

wrf_namelist['bdy_control'] = {
    'spec_bdy_width'                      : 5,
    'spec_zone'                           : 1,
    'relax_zone'                          : 4,
    'specified'                           : [True] + [False] * (max_dom-1),
    'nested'                              : [False] + [True] * (max_dom-1),
    }

wrf_namelist['grib2'] = {}

wrf_namelist['namelist_quilt'] = {
    'nio_tasks_per_group' : 0,
    'nio_groups' : 1,
    }

wrf_namelist['dfi_control'] = {
    'dfi_opt'                             : use_dfi,
    'dfi_nfilter'                         : 7,
    'dfi_write_filtered_input'            : True,
    'dfi_write_dfi_history'               : False,
    'dfi_time_dim'                        : 1000,
    }


wrf_namelist['diags'] = {
  'p_lev_diags'                         : 0,
  'num_press_levels'                    : 5,
  'press_levels'                        : [92500, 85000, 70000, 50000, 25000],
  'use_tot_or_hyd_p'                    : 2,
  }

# Now the WRF-VAR part of the namelist
# (Still in namelist.input)

wrf_namelist['wrfvar1'] = {
    'check_max_iv_print'     : False,
    'write_increments'       : False,
    }
wrf_namelist['wrfvar2'] = {}
wrf_namelist['wrfvar3'] = {}
wrf_namelist['wrfvar4'] = {
    'use_synopobs'    : False,
    'use_shipsobs'    : False,
    'use_metarobs'    : False,
    'use_soundobs'    : False,
    'use_pilotobs'    : False,
    'use_airepobs'    : False,
    'use_geoamvobs'   : False,
    'use_polaramvobs' : False,
    'use_bogusobs'    : False,
    'use_buoyobs'     : False,
    'use_profilerobs' : False,
    'use_satemobs'    : False,
    'use_gpspwobs'    : False,
    'use_gpsrefobs'   : False,
    'use_qscatobs'    : False,
    'use_radar_rv'    : False,
    'use_radar_rf'    : False,
    'use_airsretobs'  : False,
    }

wrf_namelist['wrfvar5'] = {
    'check_max_iv' : False,
    }

wrf_namelist['wrfvar6'] = {
    'max_ext_its' : 1,
    'ntmax'       : 200,
    'eps'         : 0.01,
    }

wrf_namelist['wrfvar7'] = {
    'cv_options'   : 3,
    'as1'          : [0.000001, 0.06, 0.15],
    'as2'          : [0.000001, 0.06, 0.15],
    'as3'          : [0.5, 0.10, 0.15],
    'as4'          : [0.1, 0.1, 0.3],
    'as5'          : [0.00000005, 0.10, 0.15],
    'rf_passes'    : 6,
    'var_scaling1' : 1.0,
    'var_scaling2' : 1.0,
    'var_scaling3' : 1.0,
    'var_scaling4' : 1.0,
    'var_scaling5' : 1.0,
    'len_scaling1' : 1.0,
    'len_scaling2' : 1.0,
    'len_scaling3' : 1.0,
    'len_scaling4' : 1.0,
    'len_scaling5' : 1.0,
    'je_factor'    : 1.0,
    }
wrf_namelist['wrfvar8'] = {}

wrf_namelist['wrfvar9'] = {
    'trace_csv'  : False,
    'use_html'   : False,
    }

wrf_namelist['wrfvar10'] = {}

wrf_namelist['wrfvar11'] = {
    'cv_options_hum'   : 1,
    'check_rh'         : 1,
    'set_omb_rand_fac' : 1.0,
    }

wrf_namelist['wrfvar12'] = {}


wrf_namelist['wrfvar13'] = {
    'vert_corr'     : 2,
    'vertical_ip'   : 0,
    'vert_evalue'   : 1,
    'max_vert_var1' : 99.0,
    'max_vert_var2' : 99.0,
    'max_vert_var3' : 99.0,
    'max_vert_var4' : 99.0,
    'max_vert_var5' : 0.0,
    }

wrf_namelist['wrfvar14'] = {}

wrf_namelist['wrfvar15'] = {
    'num_pseudo'    : 0,
    }

wrf_namelist['wrfvar16'] = {}

wrf_namelist['wrfvar17'] = {
    'analysis_type'  : "RANDOMCV",
    }

wrf_namelist['wrfvar18'] = {}
wrf_namelist['wrfvar19'] = {}
wrf_namelist['wrfvar20'] = {}
wrf_namelist['wrfvar21'] = {}
wrf_namelist['wrfvar22'] = {}
wrf_namelist['wrfvar23'] = {}


