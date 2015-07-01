#!/usr/bin/env python

import os, sys, getopt,datetime
sys.path.append('/home/disk/pvort/nobackup/lmadaus/WRF/DOMAINS/ens_july27/wrfdart')
from WRF_dart_param import *


calcend = 0
seed_num = 100000
add_nml_line = False

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
      model_Nx1 = model_Nx2
      model_Ny1 = model_Ny2
      



# Format times in the correct way
startdatedt = datetime.datetime.strptime(date_start,'%Y%m%d%H')
if calcend:
   enddatedt = startdatedt + lengthdt
else:
   enddatedt = datetime.datetime.strptime(date_end,'%Y%m%d%H')

wps_startdate = startdatedt.strftime('%Y-%m-%d_%H:00:00')
wps_enddate = enddatedt.strftime('%Y-%m-%d_%H:00:00')
totalrun = enddatedt - startdatedt
run_days = int(totalrun.days)
run_hours = int(totalrun.seconds/3600)
run_minutes = int((totalrun.seconds-(run_hours*3600))/60)
fYYYY = startdatedt.strftime('%Y')
fMM = startdatedt.strftime('%m')
fDD = startdatedt.strftime('%d')
fHH = startdatedt.strftime('%H')
lYYYY = enddatedt.strftime('%Y')
lMM = enddatedt.strftime('%m')
lDD = enddatedt.strftime('%d')
lHH = enddatedt.strftime('%H')



# Other variables to convert
int_seconds = dlbc * 60
grid_rat = int(model_gridspx1/model_gridspx2)

# Compute the dfi times
# Even if DFI is not requested
bdfi = startdatedt - datetime.timedelta(minutes=dfi_bckstop_window)
fdfi = startdatedt + datetime.timedelta(minutes=dfi_fwdstop_window)
bdfi_YYYY = bdfi.strftime('%Y')
bdfi_MM   = bdfi.strftime('%m')
bdfi_DD   = bdfi.strftime('%d')
bdfi_HH   = bdfi.strftime('%H')
bdfi_mm   = bdfi.strftime('%M')
fdfi_YYYY = fdfi.strftime('%Y')
fdfi_MM   = fdfi.strftime('%m')
fdfi_DD   = fdfi.strftime('%d')
fdfi_HH   = fdfi.strftime('%H')
fdfi_mm   = fdfi.strftime('%M')
dfi_len = fdfi - bdfi
dfi_cutoff_seconds = dfi_len.seconds

#Seed num for skebs
nens_skebs = str(int(seed_num/100000))+startdatedt.strftime('%m%d%H')

print "Writing namelists..."
print wps_startdate
print wps_enddate

wpsfile = open('namelist.wps','w')

print >> wpsfile,"""
&share
 wrf_core = 'ARW',
 max_dom = %(max_dom)s,
 start_date = '%(wps_startdate)s', '%(wps_startdate)s',
 end_date   = '%(wps_enddate)s', '%(wps_enddate)s',
 interval_seconds = %(int_seconds)d,
 io_form_geogrid = 2,
/

&geogrid
 parent_id         =   1,   1,
 parent_grid_ratio =   1,   %(grid_rat)d,
 i_parent_start    =   0,  %(model_origin_ll_i2)d,
 j_parent_start    =   0,  %(model_origin_ll_j2)d,
 e_we              =  %(model_Nx1)d, %(model_Nx2)d,
 e_sn              =  %(model_Ny1)d, %(model_Ny2)d,
 geog_data_res     = '30s','30s',
 dx = %(model_gridspx1)f,
 dy = %(model_gridspy1)f,
 map_proj = '%(model_proj)s',
 ref_lat   =  %(model_centlat)f,
 ref_lon   = %(model_centlon)f,
 truelat1  =  %(model_stdlat1)f,
 truelat2  =  %(model_stdlat2)f,
 stand_lon = %(model_stdlon)f,
 geog_data_path = '%(dir_src_wps_geog)s', 
/

&ungrib
 out_format = 'WPS',
 prefix = 'FILE',
/

&metgrid
 fg_name = 'FILE'
 io_form_metgrid = 2, 
/
""" % locals()
wpsfile.close()

wrfinfile = open('namelist.input','w')
print >>wrfinfile, """&time_control
  run_days           =  %(run_days)d,
  run_hours          =  %(run_hours)d,
  run_minutes        =  %(run_minutes)d,
  run_seconds        = 0,
  start_year         = %(fYYYY)s, %(fYYYY)s,
  start_month        = %(fMM)s, %(fMM)s,
  start_day          = %(fDD)s, %(fDD)s,
  start_hour         = %(fHH)s, %(fHH)s,
  start_minute       = 00, 00,
  start_second       = 00, 00,
  end_year           = %(lYYYY)s, %(lYYYY)s,
  end_month          = %(lMM)s, %(lMM)s,
  end_day            = %(lDD)s, %(lDD)s,
  end_hour           = %(lHH)s, %(lHH)s,
  end_minute         = 00, 00,
  end_second         = 00, 00,
  interval_seconds   = %(int_seconds)d,
  input_from_file    = .true.,.true.,
  history_interval   = %(wrfout_int)d, %(wrfout_int)d,
  frames_per_outfile = %(model_num_in_output)d, %(model_num_in_output)d,
  restart            = .false.,
  restart_interval   = 5000,
  io_form_history    = 2,
  io_form_restart    = 2,
  io_form_input      = 2,
  io_form_boundary   = 2,
  io_form_auxinput2  = 2,
  debug_level        = %(model_debug_level)d,
  !auxhist23_outname   = "plevel_d<domain>.nc",
  !auxhist23_interval  = %(wrfout_int)d,
  !frames_per_auxhist23= 6*2000,
  !io_form_auxhist23    = 2,
  !ncd_nofill          = .true.
/

&domains
  time_step               = %(dt)d,
  time_step_fract_num     = 0,
  time_step_fract_den     = 1,
  max_dom                 = %(max_dom)d,
  s_we                    = 1, 1,
  e_we                    = %(model_Nx1)d, %(model_Nx2)d,
  s_sn                    = 1, 1,
  e_sn                    = %(model_Ny1)d, %(model_Ny2)d,
  s_vert                  = 1, 1,
  e_vert                  = %(model_Nz1)d, %(model_Nz2)d,
  dx                      = %(model_gridspx1)f, %(model_gridspx2)f,
  dy                      = %(model_gridspy1)f, %(model_gridspy2)f,
  grid_id                 = 1, 2,
  parent_id               = 0, 1,
  i_parent_start          = 1, %(model_origin_ll_i2)d,
  j_parent_start          = 1, %(model_origin_ll_j2)d,
  parent_grid_ratio       = 1, %(grid_rat)d,
  parent_time_step_ratio  = 1, %(grid_rat)d,
  feedback                = 0,
  smooth_option           = 0,
  lagrange_order          = 1,
  interp_type             = 1,
  hypsometric_opt         = 2,
  lowest_lev_from_sfc     = .false.,
  force_sfc_in_vinterp    = 1,
  zap_close_levels        = 500
  smooth_cg_topo          = .true.,
  sfcp_to_sfcp            = .true.,
  use_levels_below_ground = .true.,
  adjust_heights          = .false.
  eta_levels              = %(model_sigma)s,
  p_top_requested         = %(model_pres_top)f,
  use_adaptive_time_step  = .true.
  step_to_output_time     = .true.
  target_cfl              = 1.2
  max_step_increase_pct   = 5
  starting_time_step      = 15
  min_time_step           = 5
  max_time_step           = 30
/

&physics
  mp_physics                          = %(model_mp_phys)d, %(model_mp_phys)d,
  do_radar_ref                        = 1,
  ra_lw_physics                       = %(model_lw_phys)d, %(model_lw_phys)d,
  ra_sw_physics                       = %(model_sw_phys)d, %(model_sw_phys)d,
  radt                                = %(model_radt)d, %(model_radt)d,
  sf_sfclay_physics                   = %(model_sfclay_phys)d, %(model_sfclay_phys)d,
  sf_surface_physics                  = %(model_surf_phys)d, %(model_surf_phys)d,
  CO2TF                               = 1,
  bl_pbl_physics                      = %(model_pbl_phys)d, %(model_pbl_phys)d,
  bldt                                = %(model_bldt)d, %(model_bldt)d,
  cu_physics                          = %(model_cu_phys)d, %(model_cu_phys)d,
  cudt                                = %(model_cudt)d, %(model_cudt)d,
  grav_settling                       = 0,
  isfflx                              = %(model_use_surf_flux)d,
  ifsnow                              = %(model_use_snow)d,
  icloud                              = %(model_use_cloud)d,
  surface_input_source                = 1,
  num_soil_layers                     = 9,
  sf_urban_physics                    = 0,
  rdlai2d                             = .true.,
  usemonalb                           = .true.,
  seaice_threshold                    = 271.4,
  mp_zero_out                         = %(model_mp_zero_out)d,
  mp_zero_out_thresh                  = %(model_mp_zero_out_thresh)s,
  maxiens                             = %(model_maxiens)d,
  maxens                              = %(model_maxens)d,
  maxens2                             = %(model_maxens2)d,
  maxens3                             = %(model_maxens3)d,
  ensdim                              = %(model_ensdim)d,
  num_land_cat                        = 21,
  mosaic_lu                           = 1,
  mosaic_soil                         = 1,
  prec_acc_dt                         = 60,
/

&fdda
/

&dynamics
  rk_ord                              = 3,
  diff_6th_opt                        = %(model_diff_6thopt)d,
  diff_6th_factor                     = %(model_diff_6thfact)f,
  w_damping                           = %(model_w_damping)d,
  diff_opt                            = %(model_diff_opt)d,
  km_opt                              = %(model_km_opt)d,
  damp_opt                            = %(model_damp_opt)d,
  zdamp                               = %(model_zdamp)f,
  base_temp                           = %(model_tbase)f,
  dampcoef                            = %(model_dampcoef)f, %(model_dampcoef)f,
  khdif                               = %(model_khdif)d, %(model_khdif)d,
  kvdif                               = %(model_kvdif)d, %(model_kvdif)d,
  smdiv                               = %(model_smdiv)f, %(model_smdiv)f,
  emdiv                               = %(model_emdiv)f, %(model_emdiv)f,
  epssm                               = %(model_epssm)f, %(model_epssm)f,
  time_step_sound                     = 4, 4,
  non_hydrostatic                     = .true., .true.,
""" % locals()
if add_nml_line:
    print >> wrfinfile,"  use_baseparam_fr_nml                = .t.,"

print >> wrfinfile, """/

&stoch
  stoch_force_opt                     = 1, 1,
  stoch_vertstruc_opt                 = 0, 0,
  perturb_bdy                         = 1
  nens                                = %(nens_skebs)s,
  tot_backscat_psi                    = 1.0E-04
  tot_backscat_t                      = 5.0E-05
  ztau_psi                            = 3600.0
  ztau_t                              = 3600.0
  rexponent_psi                       = -1.83
  rexponent_t                         = -1.83
  kminforc                            = 1
  lminforc                            = 1
  kminforct                           = 4
  lminforct                           = 4
/

&bdy_control
  spec_bdy_width                      = %(model_bdy_width)d,
  spec_zone                           = %(model_spec_zone)d,
  relax_zone                          = %(model_relax_zone)d,
  specified                           = .true., .false.,
  nested                              = .false., .true.,
/

&grib2
/

&namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
/

&dfi_control
  dfi_opt                             = %(model_dfi_opt)d,
  dfi_nfilter                         = %(model_dfi_nfilter)d,
  dfi_write_filtered_input            = %(model_dfi_write_filt_input)s,
  dfi_write_dfi_history               = .false.,
  dfi_time_dim                        = 1000,
  dfi_cutoff_seconds                  = %(dfi_cutoff_seconds)d,
  dfi_bckstop_year                    = %(bdfi_YYYY)s,
  dfi_bckstop_month                   = %(bdfi_MM)s,
  dfi_bckstop_day                     = %(bdfi_DD)s,
  dfi_bckstop_hour                    = %(bdfi_HH)s,
  dfi_bckstop_minute                  = %(bdfi_mm)s,
  dfi_bckstop_second                  = 0,
  dfi_fwdstop_year                    = %(fdfi_YYYY)s,
  dfi_fwdstop_month                   = %(fdfi_MM)s,
  dfi_fwdstop_day                     = %(fdfi_DD)s,
  dfi_fwdstop_hour                    = %(fdfi_HH)s,
  dfi_fwdstop_minute                  = %(fdfi_mm)s,
  dfi_fwdstop_second                  = 0,
 
/


&diags
  p_lev_diags                         = 0 
  num_press_levels                    = 5
  press_levels                        = 92500, 85000, 70000, 50000, 25000
  use_tot_or_hyd_p                    = 2
/

  


! Begin WRF-VAR section

&wrfvar1
  check_max_iv_print     = .false.,
  write_increments       = .false.,
/
&wrfvar2
/
&wrfvar3
/
&wrfvar4
  use_synopobs    = .false.,
  use_shipsobs    = .false.,
  use_metarobs    = .false.,
  use_soundobs    = .false.,
  use_pilotobs    = .false,
  use_airepobs    = .false.,
  use_geoamvobs   = .false.,
  use_polaramvobs = .false.,
  use_bogusobs    = .false.,
  use_buoyobs     = .false.,
  use_profilerobs = .false.,
  use_satemobs    = .false.,
  use_gpspwobs    = .false.,
  use_gpsrefobs   = .false.,
  use_qscatobs    = .false.,
  use_radar_rv    = .false.,
  use_radar_rf    = .false.,
  use_airsretobs  = .false.,
/
&wrfvar5
  check_max_iv=.false.,
/
&wrfvar6
  max_ext_its = 1,
  ntmax       = 200,
  eps         = 0.01,
/

&wrfvar7
  cv_options   = 3,
  as1          = %(bc_pscale_vel)f, %(bc_hscale_vel)f, %(bc_vscale_vel)f,
  as2          = %(bc_pscale_vel)f, %(bc_hscale_vel)f, %(bc_vscale_vel)f,
  as3          = %(bc_pscale_temp)f, %(bc_hscale_temp)f, %(bc_vscale_temp)f,
  as4          = %(bc_pscale_moist)f, %(bc_hscale_moist)f, %(bc_vscale_moist)f,
  as5          = %(bc_pscale_pres)f, %(bc_hscale_pres)f, %(bc_vscale_pres)f,
  rf_passes    = 6,
  var_scaling1 = 1.0,
  var_scaling2 = 1.0,
  var_scaling3 = 1.0,
  var_scaling4 = 1.0,
  var_scaling5 = 1.0,
  len_scaling1 = 1.0,
  len_scaling2 = 1.0,
  len_scaling3 = 1.0,
  len_scaling4 = 1.0,
  len_scaling5 = 1.0,
  je_factor    = 1.0,
/
&wrfvar8
/
&wrfvar9
  trace_csv  = .false.,
  use_html   = .false.,
/
&wrfvar10
/
&wrfvar11
  cv_options_hum   = 1,
  check_rh         = 1,
  set_omb_rand_fac = 1.0,
  seed_array1     = %(date_start)s,
  seed_array2     = %(seed_num)s,
/
&wrfvar12
/
&wrfvar13
  vert_corr     = 2,
  vertical_ip   = 0,
  vert_evalue   = 1,
  max_vert_var1 = 99.0,
  max_vert_var2 = 99.0,
  max_vert_var3 = 99.0,
  max_vert_var4 = 99.0,
  max_vert_var5 = 0.0,
/
&wrfvar14
/
&wrfvar15
  num_pseudo    = 0,
/
&wrfvar16
/
&wrfvar17
  analysis_type  = "RANDOMCV",
/
&wrfvar18
  analysis_date = "%(wps_startdate)s",
/
&wrfvar19
/
&wrfvar20
/
&wrfvar21
/
&wrfvar22
/
&wrfvar23
/

""" % locals()
wrfinfile.close()
print "Done writing namelists."
