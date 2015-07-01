#! /usr/bin/env python
########################################################################
#
#   make_namelist_dart.scr - script that makes the DART namelist 
#
#   $1 - parameter file
#
#   created Feb. 2009 R. Mahajan, University of Washington
#   modified for Python June 2011 L. Madaus, University of Washington
#
########################################################################

from WRF_dart_param import *
from WRF_dart_obtypes import *
import sys,getopt,re
from datetime import datetime, timedelta
from os import popen

timestr = date_start
write_prior_inf = False
(opts,args) = getopt.getopt(sys.argv[1:],'d:i')

for o,a in opts:
   if o == '-d':
      timestr = a
   if o == '-i':
      write_prior_inf = True

# Get the correct time
if re.search('^\d{12}$',timestr):
   assim_time = datetime.strptime(timestr,'%Y%m%d%H%M')
elif re.search('^\d{10}$', timestr):
   assim_time = datetime.strptime(timestr,'%Y%m%d%H')
else:
   print("Unrecognized time string: ", timestr)
   exit(1)

# Use the advance_time utility from DART to get the day and second
# in DART format of the assimilation window
#first_ob = os.popen('echo {:s} -{:d}m -g | {:s}/advance_time'.format(timestr,window_minutes,dir_src_dart)).read()
#last_ob = os.popen('echo {:s} +{:d}m -g | {:s}/advance_time'.format(timestr,window_minutes,dir_src_dart)).read()
epoch = datetime(1601,1,1,0)
start_window = assim_time - timedelta(minutes=window_minutes)
end_window = assim_time + timedelta(minutes=window_minutes)
first_ob = start_window - epoch
last_ob = end_window - epoch

# Forecast lengths
fct_len_days = str(int(float(fct_len_hrs) / 24))
fct_len_seconds = str(int(float(fct_len_hrs) % 24) * 3600)
bin_int_seconds = int(window_minutes)*120

# Just let us know where we are
print(timestr)


def main():
    """ Main execution of the writing.  First, populate the
    namelist dictionary and then call the function to write it """
    namelist = set_namelist_sectors()
    write_namelist(namelist)

def write_namelist(nmlst,append_nml=True,tab_width=4):
    """ Loop through the namelist dictionary and write the sections """
    tabstr = ' ' * tab_width
    with open('input.nml', 'wt') as outfile:
        # Loop through all sections
        for nml_section, sectionvals in nmlst.items():
            if append_nml:
                outfile.write('&{:s}_nml\n'.format(nml_section))
            else:
                outfile.write('&{:s}\n'.format(nml_section))
            # Now loop through each variable in the section
            for varname, value in sectionvals.items():
                outfile.write('{0:s}{1:<30s} = {2:s},\n'.format(tabstr,varname,value))
            # After all variables, close the section
            outfile.write('/\n')
            outfile.write('\n')
        outfile.close()
    # File closes automatically
    return




def set_namelist_sectors():
    """ Here we define a dictionary of dictionaries for each block of the
    namelist """

    namelist = {}
    # Here we define each block of the input.nml file separately
    namelist['perfect_model_obs'] = {}

    namelist['filter'] = {
        'async' : '{:d}'.format(async),
        'adv_ens_command' : '"./advance_model.csh"',
        'ens_size' :  '{:d}'.format(Ne),
        'start_from_restart' :  '.true.',
        'output_restart' :  '.true.',
        'obs_sequence_in_name' : '"obs_seq.prior"',
        'obs_sequence_out_name': '"obs_seq.posterior"',
        'restart_in_file_name' : '"filter_ic_old"',
        'restart_out_file_name' : '"filter_ic_new"', 
        'init_time_days' :  '-1',
        'init_time_seconds' : '-1',
        'first_obs_days' : '-1',
        'first_obs_seconds': '-1',
        'last_obs_days' : '-1',
        'last_obs_seconds'         : '-1',
        'num_output_state_members' : '{:d}'.format(Ne),
        'num_output_obs_members'   : '{:d}'.format(Ne),
        'output_interval'          : '1',
        'num_groups'               : '1',
        'input_qc_threshold'       : '{:f}'.format(input_qc_threshold),
        'outlier_threshold'        : '{:f}'.format(outlier_threshold),
        'output_forward_op_errors' : '.false.',
        'output_timestamps'        : '.false.',
        'output_inflation'         : '.true.',
        'trace_execution'          : '.true.',
        'silence'                  : '.false.',
        'inf_flavor'               : '{0:d}, {1:d}'.format(infl_flavor_prior,infl_flavor_post),
        'inf_initial_from_restart' : '{0:s}, {0:s}'.format(format_true_false(write_prior_inf)),
        'inf_output_restart'          : '.true.,   .true.,',
        'inf_deterministic'           : '.true.,   .true.,',
        'inf_in_file_name'            :  '"prior_inf_ic_old", "post_inf_ic_old"',
        'inf_out_file_name'           : '"prior_inf_ic_new", "post_inf_ic_new"',
        'inf_diag_file_name'          : '"prior_inf_diag", "post_inf_diag"',
        'inf_initial'                 : '{0:f}, {1:f}'.format(infl_mean_init_prior, infl_mean_init_post),
        'inf_sd_initial'              : '{0:f}, {1:f}'.format(infl_sd_init_prior, infl_sd_init_post),
        'inf_damping'                 : '{0:f}, {1:f}'.format(infl_damp_prior, infl_damp_post), 
        'inf_lower_bound'             : '{0:f}, {1:f}'.format(infl_lb_prior, infl_lb_post),
        'inf_upper_bound'             : '{0:f}, {1:f}'.format(infl_ub_prior, infl_ub_post),
        'inf_sd_lower_bound'          : '{0:f}, {1:f}'.format(infl_sd_lb_prior, infl_sd_lb_post),
        'direct_netcdf_read'          : format_true_false(flag_direct_netcdf_io),
        'direct_netcdf_write'         : format_true_false(flag_direct_netcdf_io),
    }

    namelist['io_filenames'] = {
        'restart_in_stub'       : '"../mems/m"',
        'restart_out_stub'      : '"../mems/m"',
        'overwrite_input'       : '.true.'
    }

    namelist['state_vector_io'] = {
        'limit_mem' : '2147483640',
        'limit_procs' : '100000',
        'create_restarts' : '.false.',
        'time_unlimited' : '.true.',
    }

    namelist['smoother'] = {
        'num_lags'              : '0',
        'start_from_restart'    : '.false.',
        'output_restart'        : '.false.',
        'restart_in_file_name'  : '"smoother_ics"',
        'restart_out_file_name' : '"smoother_restart"',
    }

    namelist['assim_tools'] = {
        'filter_kind'                   : '{:d}'.format(assim_meth),
        'cutoff'                          : cov_cutoff,
        'sort_obs_inc'                    : '.false.',
        'spread_restoration'              : '.false.',
        'sampling_error_correction'       : '.false.',
        'adaptive_localization_threshold' : '1600',
        'print_every_nth_obs'             : '100',
    }

    namelist['cov_cutoff'] = {
        'select_localization' : '{:d}'.format(assim_loc_meth) 
    }


    namelist['obs_sequence'] = {
        'write_binary_obs_sequence' : '.false.'
    }

    namelist['preprocess'] = {
        'input_obs_kind_mod_file'  : "'../../../obs_kind/DEFAULT_obs_kind_mod.F90'",
        'output_obs_kind_mod_file' : "'../../../obs_kind/obs_kind_mod.f90'",
        'input_obs_def_mod_file'   : "'../../../obs_def/DEFAULT_obs_def_mod.F90'",
        'output_obs_def_mod_file'  : "'../../../obs_def/obs_def_mod.f90'",
        'input_files'              : """'../../../obs_def/obs_def_reanalysis_bufr_mod.f90',
                                 '../../../obs_def/obs_def_altimeter_mod.f90',
                                 '../../../obs_def/obs_def_radar_mod.f90',
                                 '../../../obs_def/obs_def_metar_mod.f90',
                                 '../../../obs_def/obs_def_dew_point_mod.f90',
                                 '../../../obs_def/obs_def_gps_mod.f90',
                                 '../../../obs_def/obs_def_gts_mod.f90',
                                 '../../../obs_def/obs_def_QuikSCAT_mod.f90',
                                 '../../../obs_def/obs_def_vortex_mod.f90',
                                 '../../../obs_def/obs_def_uw_supplemental_mod.f90',
                                 '../../../obs_def/obs_def_alt_tendency_mod.f90'"""
    }

    namelist['obs_kind'] = {
        'assimilate_these_obs_types' : make_obs_list(eval=False),
        'evaluate_these_obs_types'   : make_obs_list(eval=True),
    }



    # Notes for obs_def_radar_mod_nml:
    # (1) Reflectivity limit can be applied both to observations or state (forward operator).
    # (2) Default lowest_reflectivity values DART will use (if apply_reflectivity_limit = .true.)
    #     is missing_r8. If you want to use the default, delete the line of respective
    #     lowest_reflectivity.
    # (3) As it is not clear how to assimilate Z (for now), "convert_to_dbz" is reset to .true.
    #     even if you set it to .false. here in the namelist.

    # Proceed with printing
    namelist['obs_def_radar_mod'] = {
        'apply_ref_limit_to_obs'     : '.false.',
        'reflectivity_limit_obs'     : '0.0',
        'lowest_reflectivity_obs'    : '0.0',
        'apply_ref_limit_to_fwd_op'   : '.false.',
        'reflectivity_limit_fwd_op'   : '0.0',
        'lowest_reflectivity_fwd_op'  : '0.0',
    }

    namelist['assim_model'] = {
      'write_binary_restart_files' : '.true.',
      'netCDF_large_file_support'  : '.true.',
    }

    # Notes for model_nml:
    # (1) vert_localization_coord must be one of:
    #     1 = model level
    #     2 = pressure
    #     3 = height
    # (2) see below for explanations of polar, periodic_x,
    #     periodic_y, and scm

    namelist['model'] = {
        'default_state_variables'     : '.false.',
        'wrf_state_variables'         : get_wrf_state_vars()[0],
        'wrf_state_bounds'            : get_wrf_state_vars()[1],   
        'output_state_vector'         : '.false.',
        'num_domains'                 : '{:d}'.format(max_dom),
        'calendar_type'               : '3',
        'sfc_elev_max_diff'           : '{:f}'.format(sfc_elev_tol),
        'assimilation_period_seconds' : '{:d}'.format(fct_len*60),
        'allow_obs_below_vol'         : '.false.',
        'vert_localization_coord'     : '{:d}'.format(vert_loc_coord),
        'center_search_half_length'   : '200000.',
        'center_spline_grid_scale'    : '5',
        'circulation_pres_level'      : '80000.0',
        'circulation_radius'          : '200000.0',
        'polar'                       : '.false.',
        'periodic_x'                  : '.false.',
        'periodic_y'                  : '.false.',
        'scm'                         : '.false.',
    }

    namelist['dart_to_wrf'] = {
        'model_advance_file' : '.false.',
        'dart_restart_name'  : '"dart_wrf_vector"',
        'adv_mod_command'    : '"./wrf.exe"',
    }


    namelist['wrf_to_dart'] = {
        'dart_restart_name'  : '"dart_wrf_vector"',
        'print_data_ranges'  : '.true.',
        'debug'              : '.true.',
    }


    namelist['location'] = {
        'horiz_dist_only'             : '{:s}'.format(horizontal_localization_only),
        'vert_normalization_pressure' : '{:s}'.format(vert_norm_pres), 
        'vert_normalization_height'   : '{:s}'.format(vert_norm_hght),
        'vert_normalization_level'    : '{:s}'.format(vert_norm_lev),
        'approximate_distance'        : '.false.', 
        'nlon'                        : '71',
        'nlat'                        : '36',
        'output_box_info'             : '.false.',
    }

    namelist['utilities'] = {
        'TERMLEVEL'      : '1',
        'logfilename'    : '"dart_log.out"',
        'nmlfilename'    : '"dart_log.nml"',
        'module_details' : '.false.',
    }

    namelist['reg_factor'] = {
        'select_regression'    : '1',
        'input_reg_file'       : '"time_mean_reg"',
        'save_reg_diagnostics' : '.false.',
        'reg_diagnostics_file' : '"reg_diagnostics"',
    }

    namelist['ensemble_manager'] = {
        'single_restart_file_in'  : '.false.',
        'single_restart_file_out' : '.false.',
        'perturbation_amplitude'  : '0.0',
    }

    namelist['closest_member_tool'] = {
        'input_file_name'         : '"filter_ic_old"',
        'output_file_name'        : '"closest_member_restart"',
        'ens_size'                : '{:d}'.format(Ne),
        'single_restart_file_in'  : '.false.',
        'difference_method'       : '4',
        'use_only_kinds'          : '""',
    }

    namelist['merge_obs_seq'] = {
        'num_input_files' : '2',
        'filename_seq'    : '"obs_seq.ONE", "obs_seq.TWO"',
        'filename_out'    : '"obs_seq.merged"',
    }

    namelist['obs_diag'] = {
        'obs_sequence_name'     : '"obs_seq.diag"',
        'first_bin_center'      : '{:%Y, %m, %d, %H, %M, %S}'.format(assim_time),
        'last_bin_center'       : '{:%Y, %m, %d, %H, %M, %S}'.format(assim_time),
        'bin_separation'        : '0,  0,  0,  {:s},  0,  0'.format(fct_len_hrs),
        'bin_width'             : '0,  0,  0,  {:s},  0,  0'.format(fct_len_hrs),
        'time_to_skip'          : '0,  0,  0,  0,  0,  0',
        'max_num_bins'          : '10000',
        'Nregions'              : '1',
        'lonlim1'               : '0.0',
        'lonlim2'               : '360.0',
        'latlim1'               : '-89.0',
        'latlim2'               :  '89.0', 
        'reg_names'             : '"Full Domain"',
        'print_mismatched_locs' : '.false.',
        'verbose'               : '{:s}'.format(obs_diag_verbose),
    }

    namelist['schedule'] = {
        'calendar'              : '"Gregorian"',  
        'first_bin_start'       : '{:%Y, %m, %d, %H, %M, %S}'.format(start_window), 
        'first_bin_end'         : '{:%Y, %m, %d, %H, %M, %S}'.format(end_window), 
        'last_bin_end'          : '{:%Y, %m, %d, %H, %M, %S}'.format(end_window), 
        'bin_interval_days'     : '0',
        'bin_interval_seconds'  : '{:d}'.format(bin_int_seconds),   
        'max_num_bins'          : '1',
        'print_table'           : '.true.',
    }

    namelist['obs_seq_to_netcdf'] = {
        'obs_sequence_name' : '"obs_seq.posterior"',
        'lonlim1'           : '0.0',
        'lonlim2'           : '360.0',
        'latlim1'           : '-89.9',
        'latlim2'           :  '89.9',
        'verbose'           :  '.true.',
    }

    namelist['restart_file_tool'] = {
        'input_file_name'              : '"restart_file_input"',
        'output_file_name'             : '"restart_file_output"',
        'ens_size'                     : '{:d}'.format(Ne),
        'single_restart_file_in'       : '.true.',
        'single_restart_file_out'      : '.true.',
        'write_binary_restart_files'   : '.true.',
        'overwrite_data_time'          : '.false.',
        'new_data_days'                : '-1',
        'new_data_secs'                : '-1',
        'input_is_model_advance_file'  : '.false.',
        'output_is_model_advance_file' : '.true.',
        'overwrite_advance_time'       : '.true.',
        'new_advance_days'             : '0',
        'new_advance_secs'             : '0', 
    }



    namelist['obs_sequence_tool'] = {
        'num_input_files'   : '""', 
        'filename_seq'      : '"{:%Y%m%d%H}_obs_seq.prior"'.format(assim_time),
        'filename_out'      : '"obs_seq.processed"',
        'filename_seq_list' : '"obfiles.txt"',
        'first_obs_days'    : '{:d}'.format(first_ob.days), 
        'first_obs_seconds' : '{:d}'.format(first_ob.seconds), 
        'last_obs_days'     : '{:d}'.format(last_ob.days), 
        'last_obs_seconds'  : '{:d}'.format(last_ob.seconds), 
        'obs_types'         : '""', 
        'keep_types'        : '.false.', 
        'print_only'        : '.false.', 
        'min_lat'           : '-90.0', 
        'max_lat'           : '90.0', 
        'min_lon'           : '0.0', 
        'max_lon'           : '360.0',
        'gregorian_cal'     : '.true.',
        'synonymous_copy_list' : '"observation", "NCEP observation"',
    }


    namelist['replace_wrf_fields'] = {
        'debug'                 : '.false.',
        'fail_on_missing_field' : '.false.',
        'fieldnames'            : '''"SNOWC",
                                     "ALBBCK",
                                     "TMN",
                                     "TSK",
                                     "SH20",
                                     "SMOIS",
                                     "SEAICE",
                                     "TSLB",
                                     "SST",
                                     "SNOWH",
                                     "SNOW"''',
        'fieldlist_file'        : '""',
    }

    namelist['restart_file_utility'] = {
        'input_file_name'              : '"restart_file_input"',
        'output_file_name'             : '"restart_file_output"',
        'ens_size'                     : '{:d}'.format(Ne),
        'single_restart_file_in'       : '.true.',
        'single_restart_file_out'      : '.true.',
        'write_binary_restart_files'   : '.true.',
        'overwrite_data_time'          : '.false.',
        'new_data_days'                : '-1',
        'new_data_secs'                : '-1',
        'input_is_model_advance_file'  : '.false.',
        'output_is_model_advance_file' : '.true.',
        'overwrite_advance_time'       : '.true.',
        'new_advance_days'             : '0',
        'new_advance_secs'             : '0',
    }

    namelist['wrf_obs_preproc'] = {
        'file_name_input'     : 'obs_seq.processed',
        'file_name_output'    : 'obs_seq.new',
        'overwrite_obs_time'  : '.false.',
        'obs_boundary'        : '{:f}'.format(obs_bdy),
        'increase_bdy_error'  : format_true_false(inc_bdy_obs_err),
        'maxobsfac'           : '{:f}'.format(max_obs_fac),
        'obsdistbdy'          : '{:f}'.format(obs_dist_bdy),
        'superob_aircraft'    : format_true_false(thin_aircraft),
        'superob_sat_winds'   : format_true_false(thin_sat_winds),
        'aircraft_horiz_int'  : '{:f}'.format(aircraft_horiz_dist),
        'aircraft_pres_int'   : '{:f}'.format(aircraft_vert_dist),
        'sat_wind_horiz_int'  : '{:f}'.format(sat_wind_horiz_dist),
        'sat_wind_pres_int'   : '{:f}'.format(sat_wind_vert_dist),
        'include_sig_data'    : format_true_false(soundings_sig_wndtmp),
        'sfc_elevation_check' : format_true_false(sfc_elev_check),
        'sfc_elevation_tol'   : '{:f}'.format(sfc_elev_tol),
        'obs_pressure_top'    : '{:f}'.format(obs_pres_top),
        'obs_height_top'      : '{:f}'.format(obs_hght_top), 
    }

    namelist['interpol_diag'] = {
        'in_diag_file'            : '"Diag.nc"',
        'out_diag_file'           : '"Diag_int.nc"',
        'include_slp'             : format_true_false(interpol_diag_include_slp),
        'include_wind_components' : format_true_false(interpol_diag_include_wind_components),
        'include_height_on_pres'  : format_true_false(interpol_diag_include_height_on_pres),
        'include_temperature'     : format_true_false(interpol_diag_include_temperature),
        'include_rel_humidity'    : format_true_false(interpol_diag_include_rel_humidity),
        'include_surface_fields'  : format_true_false(interpol_diag_include_surface_fields),
        'include_sat_ir_temp'     : format_true_false(interpol_diag_include_sat_ir_temp),
        'pres_levels'             : '{:s}'.format(interpol_diag_pres_levels),
        'height_levels'           : '{:s}'.format(interpol_diag_height_levels),
    }

    namelist['covariance_relax'] = {
        'restart_in_file_name'   : '"filter_ic_old"',
        'restart_out_file_name'  : '"filter_ic_new"',
        'ens_size'               : '{:d}'.format(Ne),
        'read_mean_from_diag'    : '.true.',
        'prior_scale'            : '{:f}'.format(cov_relax_prior_scale),
        'read_restart_members'   : '.false.',
        'write_restart_members'  : '.false.',
    }


    namelist['inflate_ens'] = {
        'ens_size'    : '{:d}'.format(Ne),
        'infl_scale'  : '{:f}'.format(init_assim_scale),
        'num_domains' : '{:d}'.format(max_dom),
        'inflate_tsk' : format_true_false(inflate_tsk),
        'inflate_lsm' : format_true_false(inflate_lsm),
    }

    namelist['obs_selection'] = {
        'filename_seq'       : '""',
        'filename_seq_list'  : '""',
        'num_input_files'    : '2',
        'filename_out'       : '"obs_seq.processed"',
        'selections_file'    : '"obsdef_mask.txt"',
        'print_only'         : '.false.',
        'calendar'           : '"gregorian"',
    }

    namelist['obs_seq_coverage'] = {
        'obs_sequence_name'       : '"obs_seq.in"',
        'obs_sequence_list'       : '""',
        'obs_of_interest'         : '""',
        'textfile_out'            : '"obsdef_mask.txt"',
        'netcdf_out'              : '"obsdef_mask.nc"',
        'calendar'                : '"gregorian"',
        'first_analysis'          : '{:%Y, %m, %d, %H, %M, %S}'.format(assim_time),
        'last_analysis'           : '{:%Y, %m, %d, %H, %M, %S}'.format(assim_time),
        'forecast_length_days'    :'{:s}'.format(fct_len_days),
        'forecast_length_seconds' :  '{:s}'.format(fct_len_seconds),
        'verification_interval_seconds' : '10800',
        'temporal_coverage_percent' : '100.0',
        'lonlim1'                 : '0.0',
        'lonlim2'                 : '360.0',
        'latlim1'                 : '-90.0',
        'latlim2'                 : '90.0',
        'verbose'                 : '.false.',    
    }

    namelist['obs_seq_verify'] = {
      'obs_sequence_name'       : '"obs_seq.processed"',
      'obs_sequence_list'       : '"obs_coverage_list.txt"',
      'station_template'        : '"obsdef_mask.nc"',
      'netcdf_out'              : '"forecast.nc"',
      'calendar'                : '"gregorian"',
      'obtype_string'           : '""',
      'verbose'                 : '.false.',
    }
    return namelist  




def make_obs_list(eval=False):
    """ Loop through the settings in WRF_dart_obtypes and format a list of
    observation kinds requested """
    eval = int(eval)
    obs_list = []

    # Evaluate the characteristics
    if use_obs_soundings:
        if rawin_sfc_pres[eval]:
            obs_list.append("'RADIOSONDE_SURFACE_PRESSURE'")
        if rawin_sfc_altimeter[eval]:
            obs_list.append("'RADIOSONDE_SURFACE_ALTIMETER'")
        if rawin_wind[eval]:
            obs_list.append("'RADIOSONDE_U_WIND_COMPONENT'")
            obs_list.append("'RADIOSONDE_V_WIND_COMPONENT'")
        if rawin_temp[eval]:
            obs_list.append("'RADIOSONDE_TEMPERATURE'")
        if rawin_humidity[eval]:
            obs_list.append("'RADIOSONDE_SPECIFIC_HUMIDITY'")

    # Now on to surface obs
    if use_obs_surface:
        if sfc_pressure[eval]:
            obs_list.append("'METAR_SURFACE_PRESSURE'")
        if sfc_altimeter[eval]:
            obs_list.append("'METAR_ALTIMETER'")
            obs_list.append("'LAND_SFC_ALTIMETER'")
            obs_list.append("'MARINE_SFC_ALTIMETER'")
        if sfc_alt_tendency[eval]:
            obs_list.append("'METAR_ALTIMETER_TENDENCY'")
       #if sfc_verify_alt[eval]:
       #   obs_list.append("'VERIFY_ALTIMETER',"
        if phone_altimeter[eval]:
            obs_list.append("'PHONE_ALTIMETER'")
        if wunderground_alt[eval]:
            obs_list.append("'WUNDERGROUND_ALTIMETER'")

        if sfc_wind[eval]:
            obs_list.append("'METAR_U_10_METER_WIND'")
            obs_list.append("'LAND_SFC_U_WIND_COMPONENT'")
            obs_list.append("'MARINE_SFC_U_WIND_COMPONENT'")
            obs_list.append("'METAR_V_10_METER_WIND'")
            obs_list.append("'LAND_SFC_V_WIND_COMPONENT'")
            obs_list.append("'MARINE_SFC_V_WIND_COMPONENT'")
        if wunderground_wind[eval]:
            obs_list.append("'WUNDERGROUND_U_WIND_COMPONENT'")
            obs_list.append("'WUNDERGROUND_V_WIND_COMPONENT'")
        if sfc_temp[eval]:
            obs_list.append("'METAR_TEMPERATURE_2_METER'")
            obs_list.append("'LAND_SFC_TEMPERATURE'")
            obs_list.append("'MARINE_SFC_TEMPERATURE'")
        if wunderground_temp[eval]:
            obs_list.append("'WUNDERGROUND_TEMPERATURE'")
        if sfc_humidity[eval]:
            obs_list.append("'METAR_SPECIFIC_HUMIDITY_2_METER'")
            obs_list.append("'LAND_SFC_SPECIFIC_HUMIDITY'")
            obs_list.append("'MARINE_SFC_SPECIFIC_HUMIDITY'")
        if wunderground_dewp[eval]:
            obs_list.append("'WUNDERGROUND_DEWPOINT'")
    if use_obs_pacnw:
        if pacnw_wind[eval]:
            obs_list.append("'PACNW_U_WIND_COMPONENT'")
            obs_list.append("'PACNW_V_WIND_COMPONENT'")
        if pacnw_altimeter[eval]:
            obs_list.append("'PACNW_ALTIMETER'")
        if pacnw_temp[eval]:
            obs_list.append("'PACNW_TEMPERATURE'")
        if pacnw_humidity[eval]:
            obs_list.append("'PACNW_SPECIFIC_HUMIDITY'")

    # Now for Aircraft
    if use_obs_aircraft:
        if acars_wind[eval]:
            obs_list.append("'ACARS_U_WIND_COMPONENT'")
            obs_list.append("'ACARS_V_WIND_COMPONENT'")
        if acars_temp[eval]:
            obs_list.append("'ACARS_TEMPERATURE'")
        if acars_humidity[eval]:
            obs_list.append("'ACARS_SPECIFIC_HUMIDITY'")

        if tamdar_wind[eval]:
            obs_list.append("'TAMDAR_U_WIND_COMPONENT'")
            obs_list.append("'TAMDAR_V_WIND_COMPONENT'")
        if tamdar_temp[eval]:
            obs_list.append("'TAMDAR_TEMPERATURE'")
        if tamdar_humidity[eval]:
            obs_list.append("'TAMDAR_SPECIFIC_HUMIDITY'")
    # Now for satellite cloud-top winds (CTW)
    if use_obs_ctw:
        if sat_ctw[eval]:
            obs_list.append("'SAT_U_WIND_COMPONENT'")
            obs_list.append("'SAT_V_WIND_COMPONENT'")

    # GPS Refractivity
    if use_obs_gpsro:
        if gpsro[eval]:
            obs_list.append("'GPSRO_REFRACTIVITY'")

    # NEXRAD Radar data
    if use_nexrad_data:
        if radar_vad_winds[eval]:
            obs_list.append("'VAD_U_WIND_COMPONENT'")
            obs_list.append("'VAD_V_WIND_COMPONENT'")
            #obs_list.append("'VAD_W_WIND_COMPONENT',"
        if radar_reflectivity[eval]:
            obs_list.append("'RADAR_REFLECTIVITY'")
        if radar_radial_velocity[eval]:
            obs_list.append("'DOPPLER_RADIAL_VELOCITY'")
        if radar_clear_reflectivity[eval]:
            obs_list.append("'RADAR_CLEARAIR_REFLECTIVITY'")
        if radar_drop_fall_speed[eval]:
            obs_list.append("'PRECIPITATION_FALL_SPEED'")

    # Dropsondes
    if use_obs_dropsondes:
        if drop_pres[eval]:
            obs_list.append("'DROPSONDE_PRESSURE'")
        if drop_wind[eval]:
            obs_list.append("'DROPSONDE_U_WIND_COMPONENT'")
            obs_list.append("'DROPSONDE_V_WIND_COMPONENT'")
        if drop_temp[eval]:
            obs_list.append("'DROPSONDE_TEMPERATURE',")
        if drop_humidity[eval]:
            obs_list.append("'DROPSONDE_SPECIFIC_HUMIDITY'")

    # Driftsondes
    if use_obs_driftsondes:
        if drift_wind[eval]:
            obs_list.append("'DRIFTSONDE_U_WIND_COMPONENT'")
            obs_list.append("'DRIFTSONDE_V_WIND_COMPONENT'")
        if drift_temp[eval]:
            obs_list.append("'DRIFTSONDE_TEMPERATURE'")
        if drift_humidity[eval]:
            obs_list.append("'DRIFTSONDE_SPECIFIC_HUMIDITY'")
       

    # TC Info
    if use_obs_tcinfo:
        if vortex_loc[eval]:
            obs_list.append("'VORTEX_LAT'")
            obs_list.append("'VORTEX_LON'")
        if vortex_pmin[eval]:
            obs_list.append("'VORTEX_PMIN'")
        if vortex_wmax[eval]:
            obs_list.append("'VORTEX_WMAX'")

    # Format the obs list
    obs_text = ",\n          ".join(obs_list)
    return obs_text

def get_wrf_state_vars():
    """ Gets the requested state vars based on parameters in WRF_dart_param """
    if minimal_vars:
        default_vars = '''    'PH','KIND_GEOPOTENTIAL_HEIGHT','TYPE_GZ','UPDATE','999',
                              'QVAPOR','KIND_VAPOR_MIXING_RATIO','TYPE_QV','UPDATE','999',
                              'T','KIND_POTENTIAL_TEMPERATURE','TYPE_T','UPDATE','999',
                              'MU','KIND_PRESSURE','TYPE_MU','UPDATE','999',
                              '''
    else:
        default_vars = '''    'U','KIND_U_WIND_COMPONENT','TYPE_U','UPDATE','999',
                                'V','KIND_V_WIND_COMPONENT','TYPE_V','UPDATE','999',
                                'W','KIND_VERTICAL_VELOCITY','TYPE_W','UPDATE','999',
                                'T','KIND_POTENTIAL_TEMPERATURE','TYPE_T','UPDATE','999',
                                'PH','KIND_GEOPOTENTIAL_HEIGHT','TYPE_GZ','UPDATE','999',
                                'MU','KIND_PRESSURE','TYPE_MU','UPDATE','999',
                                '''


    # Decide what additional variables to update
    if update_surf:
        default_vars = default_vars + ''''U10','KIND_U_WIND_COMPONENT','TYPE_U10','UPDATE','999',
                                'V10','KIND_V_WIND_COMPONENT','TYPE_V10','UPDATE','999',
                                'T2','KIND_TEMPERATURE','TYPE_T2','UPDATE','999',
                                'TH2','KIND_POTENTIAL_TEMPERATURE','TYPE_TH2','UPDATE','999',
                                'Q2','KIND_SPECIFIC_HUMIDITY','TYPE_Q2','UPDATE','999',
                                '''
                   
    if update_psfc:
        default_vars += "'PSFC','KIND_PRESSURE','TYPE_PS','UPDATE','999',\n"
    if flag_compute_tendency:
        default_vars +=  "'ALT_TEND','KIND_ALTIMETER_TENDENCY','TYPE_ALTTEND','UPDATE','999',\n"
    if update_tsk:
        default_vars += "'TSK','KIND_SKIN_TEMPERATURE','TYPE_TSK','UPDATE','999',\n"
    if update_lsm:
        # Updates Soil Temperature (TSLB), Soil Moisture (SMOIS), Unfrozen soil
        # moisture content (SH2O)
        default_vars += "                                 'TSLB','KIND_SOIL_TEMPERATURE', 'TYPE_TSLB','UPDATE','999',\n"
        default_vars += "                                 'SMOIS','KIND_SOIL_MOISTURE', 'TYPE_SMOIS','UPDATE','999',\n"
        default_vars += "                                 'SH20','KIND_SOIL_LIQUID_WATER', 'TYPE_SH2O','UPDATE','999',\n"
    # Update appropriate number of moisture variables in state vector
    if num_moist_vars >= 3:
        if not minimal_vars:
            default_vars +=  """                                'QVAPOR','KIND_VAPOR_MIXING_RATIO','TYPE_QV','UPDATE','999',
                                'QCLOUD','KIND_CLOUD_LIQUID_WATER','TYPE_QC','UPDATE','999',
                                'QRAIN','KIND_RAINWATER_MIXING_RATIO','TYPE_QR','UPDATE','999'"""

        wrf_state_bounds = """ 'QVAPOR','0.0','NULL','CLAMP',
                             'QRAIN','0.0','NULL','CLAMP',
                             'QCLOUD','0.0','NULL','CLAMP'"""
        if num_moist_vars == 3:
            default_vars +="""\n"""
            wrf_state_bounds +="""\n"""
        else:
            default_vars +=""",\n"""
            wrf_state_bounds +=""",\n"""

    if num_moist_vars >= 5:
        if not minimal_vars:
            default_vars +=  """                                'QICE','KIND_CLOUD_ICE','TYPE_QI','UPDATE','999',
                                'QSNOW','KIND_SNOW_MIXING_RATIO','TYPE_QS','UPDATE','999'"""
        wrf_state_bounds +=  """ 'QICE','0.0','NULL','CLAMP',
                             'QSNOW','0.0','NULL','CLAMP'"""
        if num_moist_vars == 5: 
            default_vars +="""\n"""
            wrf_state_bounds +="""\n"""
        else:
            default_vars +=""",\n"""
            wrf_state_bounds +=""",\n"""



    if num_moist_vars >= 6:
        if not minimal_vars:
            default_vars +=  """                               'QGRAUP','KIND_GRAUPEL_MIXING_RATIO','TYPE_QG','UPDATE','999'"""
        wrf_state_bounds += """ 'QGRAUP','0.0','NULL','CLAMP'"""

    return default_vars, wrf_state_bounds



def format_true_false(input):
    """ Quick function to write a Fortran friendly ".true." or ".false."
    based on the Python True or False given """
    if input:
        return '.true.'
    else:
        return '.false.'


if __name__ == '__main__':
    main()
