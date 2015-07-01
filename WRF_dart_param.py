#! /usr/bin/python

import os
#*******************************************************************************
#  
#  Use this file to modify all the parameters needed to run an experiment.  
#  Each experiment should have one of these files
#
#*******************************************************************************

exp_name='july27'                # Name of the experiment
 
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

flag_pert_bcs         = True      # True to run pert_wrf_bc.  If false, will still run
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

model_Nx1       = 560           # number of x grid points
model_Ny1       = 420           # number of y grid points
model_Nz1       = 51            # number of grid points in z (specify levels below)
model_gridspx1  = 3000          # gridspacing in x (in meters!!!)
model_gridspy1  = 3000          # gridspacing in y (in meters!!!)
dt              = 20           # model time step (in sec)
	
model_proj      = 'lambert'     # model projection
model_centlat   = 39.5           # Center latitude of the model (don't forget (.))
model_centlon   = -80.25         # center longitude of model
model_stdlat1   = 38.5           # standard latitude 1
model_stdlat2   = 38.5          # standard latitude 2
model_stdlon    = -92.5        # standard longitude

#**************************************************************
#
#   This set of parameters deals with nested domains,
#   moving and fixed
#
#**************************************************************

max_dom            = 1         # maximum number of domains

model_Nx2          = 271        # number of x grid points
model_Ny2          = 226        # number of y grid points
model_Nz2          = model_Nz1  # number of y grid points
model_gridspx2     = 4000      # gridspacing in x (in meters!!!)
model_gridspy2     = 4000      # gridspacing in y (in meters!!!)
model_origin_ll_i2 = 93         # lower-left i of nested domain in parent
model_origin_ll_j2 = 50         # lower-left j of nested domain in parent
model_origin_ur_i2 = 61         # upper-right i of nested domain in parent
model_origin_ur_j2 = 48         # upper-right j of nested domain in parent

model_radt2        = 30         # radiation time step for the nest
model_cu_phys2     = 1          # cumulus scheme for the nest
	
move_nest_domain   = False      # True to use a moving nest domain (NOT ENABLED)
vortex_interval    = 15         # time to update vortex position (NOT ENABLED)
vortex_speed       = 20          # maximum speed of vortex (NOT ENABLED)

#**************************************************************
#
#   This set of parameters deals with the EnKF parameters 
#   itself.  For more information on some of them, see the 
#   original papers or Hamill's review paper.
#
#**************************************************************

date_start        = '2014072600' # Start of the experiment (YYYYMMDDHH)
date_end          = '2014072912' # End of the experiment (YYYYMMDDHH)
Ne                = 50           # Number of ensemble members
fct_len           = 60          # Time between 2 assimilation (in MINUTES)
N_assim           = 500          # Number of assimilations (spinup=0)
N_assim_max       = 500          # Maximum number of assimilations 
dlbc              = 60        	 # Interval of global boundary conditions (in MINUTES)
wrfout_int        = 60           # Interval to write wrfout files (in MINUTES)

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
#   Change the settings in this section if you are generating 
#   ensemble boundary conditions from WRF 3D-Var perturbations
#   assim_bmeth = 3
#
#**************************************************************
# assim_bmeth = 3 is the only type supported.  Can use BCs from a parent domain,
# and these parameters won't affect that.
# These paramters only affect the pert_wrf_bc sequence

autocorr_fcp  = 0.50           # autocorrelation coefficient between two boundary times
fcp_scale_bp  = 1.60           # factor to scale boundary perturbations
fcp_cov_scale = 2.0            # horizontal scale of perturbations (0.75 def)
fcp_npert     = str(4*int(Ne)) # number of perturbations to generate
bc_pscale     = 0.80           # perturbation scale for WRF-VAR 3 hr forecast (formerly bc_pert_scale)
bc_hscale     = 2.00           # horizontal scale for WRF-VAR 3 hr forecast
bc_vscale     = 1.50           # vertical scale for WRF-VAR 3 hr forecast
bc_pscale_vel   = 0.000001           # perturbation scale for WRF-VAR 3 hr forecast (formerly bc_pert_scale)
bc_hscale_vel   = 0.06           # horizontal scale for WRF-VAR 3 hr forecast
bc_vscale_vel   = 0.15           # vertical scale for WRF-VAR 3 hr forecast
bc_pscale_temp  = 0.5           # perturbation scale for WRF-VAR 3 hr forecast (formerly bc_pert_scale)
bc_hscale_temp  = 0.10           # horizontal scale for WRF-VAR 3 hr forecast
bc_vscale_temp  = 0.15           # vertical scale for WRF-VAR 3 hr forecast
bc_pscale_moist = 0.1           # perturbation scale for WRF-VAR 3 hr forecast (formerly bc_pert_scale)
bc_hscale_moist = 0.10           # horizontal scale for WRF-VAR 3 hr forecast
bc_vscale_moist = 0.3           # vertical scale for WRF-VAR 3 hr forecast
bc_pscale_pres  = 0.00000005           # perturbation scale for WRF-VAR 3 hr forecast (formerly bc_pert_scale)
bc_hscale_pres  = 0.10           # horizontal scale for WRF-VAR 3 hr forecast
bc_vscale_pres  = 0.15           # vertical scale for WRF-VAR 3 hr forecast






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

dir            = '/home/disk/pvort/lmadaus/nobackup/WRF' # MAIN DIRECTORY
dir_wrf_dom    = dir + '/DOMAINS/' + exp_name            # Directory where main experiment is run 
dir_longsave   = dir_wrf_dom + '/longsave'               # storage directory 
dir_obs        = dir_wrf_dom + '/obs'                    # repository of obs.
dir_utils      = dir + '/UTILS/bin'                      # location of utility codes
dir_members    = dir_wrf_dom + '/mems'                   # members directory
dir_assim      = dir_wrf_dom + '/assimilation'           # Directory where DART will be run

# Directories where WRF, WRFDA and DART are found
WRFVARDIR         = dir + '/WRFDA'
WRFRUNDIR         = dir + '/WRFV3/run'
dir_src_wps       = dir + '/WPS'                    # WPS location
dir_src_wps_geog  = dir + '/DATA/geog'              # Location of geo data for WPS
dir_src_wrf       = dir + '/WRFV3/main'             # WRF location
dir_src_wrfvar    = dir + '/WRFDA/var'                  # WRF-VAR location
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

cluster_name        = 'enkf'                   # name of cluster nodes
mpi_run_command     = '/usr/rels/openmpi/bin/mpirun' # Command to run MPI
queue_members       = 'reg'                    # Queue to run members in
queue_filter        = 'reg'                    # Queue to run filter in          
mpi_numprocs_member = 16                       # Number of processors for member
mpi_numprocs_filter = 128                       # Number of processors for filter
mpi_numprocs_flag   = '-np %d' % mpi_numprocs_member      # Flag for numprocs in code
                                                          # for member.  Bluefire does
                                                          # not specify numprocs in the
                                                          # MPI run command, so make this
                                                          # blank for bluefire

# Extra parameters for running on Bluefire
NCAR_GAU_ACCOUNT     = '0'                   # Account to charge to at NCAR
ADVANCE_TIME_FILTER  = '0:45'                # Estimate of time for filter to run
ADVANCE_TIME_MEMBER  = '0:20'                # Estimate of time for a single member to run 
ADVANCE_QUEUE_FILTER = queue_filter          # Name of NCAR queue to use for filter 
ADVANCE_QUEUE_MEMBER = queue_members         # Name of NCAR queue to use for members
ADVANCE_CORES_FILTER = mpi_numprocs_filter  # Number of cores to use (multiples of 32)
ADVANCE_CORES_MEMBER = mpi_numprocs_member  # Number of cores to use (multiples of 32)
NCAR_ADVANCE_PTILE   = '32'                  # How many processes per core on Bluefire



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
	
# physics options
model_mp_phys        = 8            	# microphysics option 
model_lw_phys        = 4		# model long wave scheme  (1=RRTM, 3=CAM) 
model_sw_phys        = 4        	# model short wave scheme (1=RRTM, 3=CAM)
model_radt           = 30		# radiation time step (in minutes)
model_sfclay_phys    = 5		# surface layer physics (if = 1, see model_use_surf_flux below)
model_surf_phys      = 3		# land surface model
model_pbl_phys       = 5		# pbl physics   
model_bldt           = 0		# boundary layer time steps (0 : each time steps, in min)
model_cu_phys        = 0     		# cumuls param
model_cudt           = 5          	# cumuls time step
model_use_surf_flux  = 1   	        # 1 is yes, 0 is no (only works when model_sfclay_phys=1)
model_use_snow       = 1		# 1 is yes, 0 is no (only works when model_sfclay_phys=1)
model_use_cloud      = 1		# for lw_phys = 1 and sw_phys = 1 only.
model_sf_urban_phys  = 0		# urban physics (0, 1, 2)
model_mp_zero_out    = 2		# for non-zero mp_physics options, DO NOT CHANGE
model_mp_zero_out_thresh= '1.e-12'
model_maxiens        = 1		# DO NOT CHANGE THESE
model_maxens         = 3		# DO NOT CHANGE THESE
model_maxens2        = 3		# DO NOT CHANGE THESE
model_maxens3        = 16		# DO NOT CHANGE THESE
model_ensdim         = 144		# DO NOT CHANGE THESE
if model_lw_phys == 3 and model_sw_phys == 3:
	model_cam_abs_freq_s  = 21600   # FOR CAM RADIATION SCHEME ONLY 
	model_levsiz          = 59      # FOR CAM RADIATION SCHEME ONLY
	model_paerlev         = 29      # FOR CAM RADIATION SCHEME ONLY 
	model_cam_abs_dim1    = 4       # FOR CAM RADIATION SCHEME ONLY 

# dynamics options
model_w_damping      = 1                # vertical velocity damping (1=ON, 0=OFF)
model_diff_opt       = 1                # turbulence and mixing option
model_diff_6thopt    = 2                # numerical hyperdiffusion option
model_diff_6thfact   = 0.25	        # numerical diffusion factor
model_damp_opt       = 3                # upper level damping option
model_dampcoef       = 0.2              # damping coefficient
model_zdamp          = 5000.            # damping depth from model top (in meters)
model_km_opt         = 4                # eddy coeff. option
model_khdif          = 0                # horiz. diffusion const.
model_kvdif          = 0                # vert. diffusion const.
model_smdiv          = 0.1              # divergence damping
model_emdiv          = 0.01             # external mode filter
model_epssm          = 0.1              # time-off centering for vert. sound waves
model_pd_moist       = 1                # pos. definite advection for moisture fields

# time control options
model_debug_level    = 0 		# WRF model debug_level 
model_num_in_output  = 1 

# domains options
# length of model_sigma must match model_Nz above
#model_sigma='1.000, 0.995, 0.990, 0.985, 0.980, 0.970, 0.960, 0.950, 0.940, 0.930, 0.920, 0.910, 0.900, 0.880, 0.860, 0.830, 0.800, 0.770, 0.740, 0.710, 0.680, 0.640, 0.600, 0.560, 0.520, 0.480, 0.440, 0.400, 0.360, 0.320, 0.280, 0.240, 0.200, 0.160, 0.120, 0.080, 0.040, 0.000' 
model_sigma='1.0000, 0.9980, 0.9940, 0.9870, 0.9750, 0.9590, 0.9390, 0.9160, 0.8920, 0.8650, 0.8350, 0.8020, 0.7660, 0.7270, 0.6850, 0.6400, 0.5920, 0.5420, 0.4970, 0.4565, 0.4205, 0.3877, 0.3582, 0.3317, 0.3078, 0.2863, 0.2670, 0.2496, 0.2329, 0.2188, 0.2047, 0.1906, 0.1765, 0.1624, 0.1483, 0.1342, 0.1201, 0.1060, 0.0919, 0.0778, 0.0657, 0.0568, 0.0486, 0.0409, 0.0337, 0.0271, 0.0209, 0.0151, 0.0097, 0.0047, 0.0000,'
model_pres_top =  5000.             # Fills in p_top_requested   
  
# boundary control options
model_spec_zone     = 1                 # boundary points specified by outer model
model_relax_zone    = 4                 # boundary points with mix of dynamics, outer model
model_bdy_width     = model_spec_zone + model_relax_zone

# digital filtering options
model_dfi_opt       = 0                 # DFI Scheme (0: NO DFI, 1-3)
model_dfi_nfilter   = 7                 # filter type (0-7)
model_dfi_write_filt_input = '.true.'   # write filtered wrfinput file at initialization time
dfi_bckstop_window  = 30                # integrate backward in MINUTES
dfi_fwdstop_window  = 15                # integrate forward  in MINUTES

#************************************
# DO NOT MODIFY BELOW
#************************************

# Some useful nco utilities
ncks_bin='`which ncks`'
ncra_bin='`which ncra`'
ncwa_bin='`which ncwa`'
ncflint_bin='`which ncflint`'
ncrcat_bin='`which ncrcat`'
ncrename_bin='`which ncrename`'

# Calculate the number of integration time steps
time_step=str((float(fct_len)*60)/float(dt))
assim_len=str(float(N_assim)*float(fct_len))
fct_len_hrs=str(float(fct_len)/60)
model_grid_ratio=str(float(model_gridspx2)/float(model_gridspx1))
model_timestep_ratio=model_grid_ratio
gmodnam = 'gfs'
gmodnamu=os.popen('echo %s | tr "[:lower:]" "[:upper:]"' % gmodnam).readlines()[0][0:-1]
username=os.popen('whoami').readlines()[0][0:-1]      # user running this experiment

# Model base temperature
model_tbase=300.

# Decide on how many soil layers to use depending on LSM model
if model_surf_phys == '1':
	model_soil_layers='5'
elif model_surf_phys == '2':
  model_soil_layers='4'
elif model_surf_phys == '3':
  model_soil_layers='6'
else:
  model_soil_layers='0'



# boundary zone width and other parameters
assim_bzw='model_spec_zone+model_relax_zone'  
assim_lgts='anal_len_bc'
dlbc_hrs=str(float(dlbc)/60.0)


# These variable are used for wrfenkf.nl, their role is same, but names have changed
bc_pert_scale=bc_pscale

# parallel_model flag
if async == 2: 
	parallel_model='false'
elif async == 4: 
	parallel_model='true'

