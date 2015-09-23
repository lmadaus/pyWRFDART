#! /usr/bin/env python
from __future__ import print_function, division
########################################################################
#
#   make_namelist_dart.scr - script that makes the DART namelist 
#
#   $1 - parameter file
#
#   created Feb. 2009 R. Mahajan, University of Washington
#   modified for Python June 2011 L. Madaus, University of Washington
#   modified for CM1 Sept 2015 L. Madaus, University of Washington
#
########################################################################

from ens_dart_param import *
from ens_dart_obtypes import *
import sys,getopt,re
from datetime import datetime, timedelta
from os import popen

timestr = cycle_len 
write_prior_inf = False
(opts,args) = getopt.getopt(sys.argv[1:],'d:i')

for o,a in opts:
   if o == '-d':
      timestr = a
   if o == '-i':
      write_prior_inf = True

assim_time_sec = int(timestr) * 60
# Use the advance_time utility from DART to get the day and second
# in DART format of the assimilation window
#first_ob = os.popen('echo {:s} -{:d}m -g | {:s}/advance_time'.format(timestr,window_minutes,dir_src_dart)).read()
#last_ob = os.popen('echo {:s} +{:d}m -g | {:s}/advance_time'.format(timestr,window_minutes,dir_src_dart)).read()
epoch = datetime(1601,1,1,0)
assim_time_date = epoch + timedelta(seconds=assim_time_sec) 
start_window = assim_time_date - timedelta(minutes=window_minutes)
end_window = assim_time_date + timedelta(minutes=window_minutes)
first_ob = start_window - epoch
last_ob = end_window - epoch

# Forecast lengths
cycle_len_days = 0 
cycle_len_seconds = int(cycle_len) * 60
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
    namelist['perfect_model_obs'] = {
        'adv_ens_command' : '../shell_scripts/advance_model.csh',
        'async'                    : '0',
        'first_obs_days'           : '-1',
        'first_obs_seconds'        : '-1',
        'init_time_days'           : '0',
        'init_time_seconds'        : '0',
        'last_obs_days'            : '-1',
        'last_obs_seconds'         : '-1',
        'obs_seq_in_file_name'     : 'obs_seq.in',
        'obs_seq_out_file_name'    : 'obs_seq.out',
        'output_forward_op_errors' : '.false',
        'output_interval'          : '1',
        'output_restart'           : '1',
        'output_timestamps'        : '.false.',
        'print_every_nth_obs'      : '-1',
        'restart_in_file_name'     : 'perfect_ics',
        'restart_out_file_name'    : 'perfect_restart',
        'silence'                  : '.false.',
        'start_from_restart'       : '.true.',
        'trace_execution'          : '.false.',

    }

    namelist['filter'] = {
        'async' : '{:d}'.format(async),
        'adv_ens_command' : '"./advance_model.csh"',
        'diagnostic_files' : '.true.',
        'direct_netcdf_read'          : format_true_false(flag_direct_netcdf_io),
        'direct_netcdf_write'         : format_true_false(flag_direct_netcdf_io),
        'ens_size' :  '{:d}'.format(Ne),
        'first_obs_days' : '-1',
        'first_obs_seconds': '-1',
        'inf_damping'                 : '{0:f}, {1:f}'.format(infl_damp_prior, infl_damp_post), 
        'inf_deterministic'           : '.true.,   .true.,',
        'inf_diag_file_name'          : '"prior_inf_diag", "post_inf_diag"',
        'inf_flavor'               : '{0:d}, {1:d}'.format(infl_flavor_prior,infl_flavor_post),
        'inf_in_file_name'            :  '"prior_inf_ic_old", "post_inf_ic_old"',
        'inf_initial'                 : '{0:f}, {1:f}'.format(infl_mean_init_prior, infl_mean_init_post),
        'inf_initial_from_restart' : '{0:s}, {0:s}'.format(format_true_false(write_prior_inf)),
        'inf_lower_bound'             : '{0:f}, {1:f}'.format(infl_lb_prior, infl_lb_post),
        'inf_out_file_name'           : '"prior_inf_ic_new", "post_inf_ic_new"',
        'inf_output_restart'          : '.true.,   .true.,',
        'inf_sd_initial'              : '{0:f}, {1:f}'.format(infl_sd_init_prior, infl_sd_init_post),
        'inf_sd_initial_from_restart' : '{0:s}, {0:s}'.format(format_true_false(write_prior_inf)),
        'inf_sd_lower_bound'          : '{0:f}, {1:f}'.format(infl_sd_lb_prior, infl_sd_lb_post),
        'inf_upper_bound'             : '{0:f}, {1:f}'.format(infl_ub_prior, infl_ub_post),
        'init_time_days' :  '-1',
        'init_time_seconds' : '-1',
        'last_obs_days' : '-1',
        'last_obs_seconds'         : '-1',
        'num_groups'               : '1',
        'num_output_obs_members'   : '{:d}'.format(Ne),
        'num_output_state_members' : '{:d}'.format(Ne),
        'obs_sequence_in_name' : '"obs_seq.prior"',
        'obs_sequence_out_name': '"obs_seq.posterior"',
        'output_forward_op_errors' : '.false.',
        'output_interval'          : '1',
        'output_restart' :  '.true.',
        'output_timestamps'        : '.false.',
        'restart_in_file_name' : '"filter_ic_old"',
        'restart_out_file_name' : '"filter_ic_new"', 
        'silence'                  : '.false.',
        'start_from_restart' :  '.true.',
        'trace_execution'          : '.true.',
    }
    namelist['model_mod_check'] = {
        'run_test'    : '2',
        'num_ens'     : '1',
    }

    namelist['io_filenames'] = {
        'restart_in_stub'       : '"../mems/m"',
        'restart_out_stub'      : '"../mems/m"',
        'overwrite_input'       : '.true.'
    }

    namelist['state_vector_io'] = {
        'perturbation_amplitude'  : '0.2',
        'single_restart_file_out' : '.true.',
        'write_binary_restart_files' : '.false.',
    }


    namelist['smoother'] = {
        'num_lags'              : '0',
        'start_from_restart'    : '.false.',
        'output_restart'        : '.false.',
        'restart_in_file_name'  : '"smoother_ics"',
        'restart_out_file_name' : '"smoother_restart"',
    }

    namelist['assim_tools'] = {
        'adaptive_localization_threshold' : '-1',
        'cutoff'                          : cov_cutoff,
        'filter_kind'                   : '{:d}'.format(assim_meth),
        'gaussian_likelihood_tails'       : '.false.',
        'localization_diagnostics_file'   : 'localization_diagnostics',
        'output_localization_diagnostics' : '.false.',
        'print_every_nth_obs'             : '100',
        'rectangular_quadrature'          : '.true.',
        'sampling_error_correction'       : '.false.',
        'sort_obs_inc'                    : '.false.',
        'spread_restoration'              : '.false.',
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
                                 '../../../obs_def/obs_def_alt_tendency_mod.f90'""",
        'overwrite_output'       : '.true.',
    
    }

    namelist['quality_control'] = {
        'enable_special_outlier_code' : '.false.',
        'input_qc_threshold'          : '3.0',
        'outlier_threshold'           : '-1.0',
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
      #'write_binary_restart_files' : '.true.',
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
        'output_state_vector'         : '.false.',
        'assimilation_period_days'    : '0',
        'assimilation_period_seconds' : '{:d}'.format(cycle_len*60),
        'model_perturbation_amplitude':'0.2',
        'model_restart_dirname'       : '.',
        'grid_filename'               : 'cm1out_rst_000001.nc',
        'calendar'                    : '"Gregorian"',  
        'debug'                       : '1',

        'model_variables'         : get_model_state_vars()[0],
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
        'bin_width_days'        : '-1',
        'bin_width_seconds'     : '-1',
        'create_rank_histogram' : '.true.',
        'int_skip_days'         : '0',
        'int_skip_seconds'      : '0',
        'lonlim1'               : '0.0',
        'lonlim2'               : '360.0',
        'Nregions'              : '1',
        'obs_sequence_name'     : '"obs_seq.diag"',
        'outliers_in_histogram' : '.true.',
        'reg_names'             : 'whole',
        'trusted_obs'           : 'null',
        'use_zero_error_obs'    : '.false',
        'verbose'               : '{:s}'.format(obs_diag_verbose),
    }

    namelist['schedule'] = {
        'calendar'              : '"Gregorian"',  
        'first_bin_start'       : '{:d}, {:d}, {:d}, {:d}, {:d}, {:d}'.format(\
            start_window.year, start_window.month, start_window.day,\
            start_window.hour, start_window.minute, start_window.second), 
        'first_bin_end'       : '{:d}, {:d}, {:d}, {:d}, {:d}, {:d}'.format(\
            end_window.year, end_window.month, end_window.day,\
            end_window.hour, end_window.minute, end_window.second), 
        'last_bin_end'       : '{:d}, {:d}, {:d}, {:d}, {:d}, {:d}'.format(\
            end_window.year, end_window.month, end_window.day,\
            end_window.hour, end_window.minute, end_window.second), 
        'bin_interval_days'     : '0',
        'bin_interval_seconds'  : '{:d}'.format(bin_int_seconds),   
        'max_num_bins'          : '1',
        'print_table'           : '.true.',
    }

    namelist['obs_seq_to_netcdf'] = {
        'append_to_netcdf'  : '.false.',
        'obs_sequence_name' : '"obs_seq.posterior"',
        'lonlim1'           : '0.0',
        'lonlim2'           : '360.0',
        'verbose'           :  '.true.',
    }

    namelist['restart_file_tool'] = {
        'ens_size'                     : '{:d}'.format(Ne),
        'gregorian_calendar'           : '.false.',
        'input_file_name'              : '"restart_file_input"',
        'input_is_model_advance_file'  : '.false.',
        'new_advance_days'                : '-1',
        'new_advance_secs'                : '-1',
        'new_data_days'                : '-1',
        'new_data_secs'                : '-1',
        'output_file_name'             : '"restart_file_output"',
        'output_is_model_advance_file' : '.true.',
        'overwrite_advance_time'       : '.true.',
        'overwrite_data_time'          : '.false.',
        'single_restart_file_in'       : '.true.',
        'single_restart_file_out'      : '.true.',
        'write_binary_restart_files'   : '.true.',
    }



    namelist['obs_sequence_tool'] = {
        'filename_out'      : '"obs_seq.processed"',
        'filename_seq'      : '"{:d}_obs_seq.prior"'.format(assim_time_sec),
        'filename_seq_list' : '""',
        'first_obs_days'    : '{:d}'.format(first_ob.days), 
        'first_obs_seconds' : '{:d}'.format(first_ob.seconds), 
        'gregorian_cal'     : '.true.',
        'last_obs_days'     : '{:d}'.format(last_ob.days), 
        'last_obs_seconds'  : '{:d}'.format(last_ob.seconds), 
        'print_only'        : '.false.', 
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
        'num_domains' : '1',
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
        'first_analysis'       : '{:d}, {:d}, {:d}, {:d}, {:d}, {:d}'.format(\
            assim_time_date.year, assim_time_date.month, assim_time_date.day,\
            assim_time_date.hour, assim_time_date.minute, assim_time_date.second), 
        'last_analysis'       : '{:d}, {:d}, {:d}, {:d}, {:d}, {:d}'.format(\
            assim_time_date.year, assim_time_date.month, assim_time_date.day,\
            assim_time_date.hour, assim_time_date.minute, assim_time_date.second), 
        'forecast_length_days'    :'{:d}'.format(cycle_len_days),
        'forecast_length_seconds' :  '{:d}'.format(cycle_len_seconds),
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

def get_model_state_vars():
    """ Gets the requested state vars based on parameters in WRF_dart_param """
    if minimal_vars:
        default_vars = '''    'qv','KIND_VAPOR_MIXING_RATIO','UPDATE','NULL','NULL',
                              'theta','KIND_POTENTIAL_TEMPERATURE','UPDATE','NULL','NULL',
                              'ppi','KIND_PRESSURE','UPDATE','NULL','NULL',
                              '''
    else:
        default_vars = '''    'ua','KIND_U_WIND_COMPONENT','UPDATE','NULL','NULL',
                                'va','KIND_V_WIND_COMPONENT','UPDATE','NULL','NULL',
                                'wa','KIND_VERTICAL_VELOCITY','UPDATE','NULL','NULL',
                                'theta','KIND_POTENTIAL_TEMPERATURE','UPDATE','NULL','NULL',
                                'ppi','KIND_PRESSURE','UPDATE','NULL','NULL',
                                '''


    # Decide what additional variables to update
    if update_surf:
        default_vars = default_vars + ''' 'u10','KIND_U_WIND_COMPONENT','UPDATE','NULL','NULL',
                                'v10','KIND_V_WIND_COMPONENT','UPDATE','NULL','NULL',
                                't2','KIND_TEMPERATURE','UPDATE','NULL','NULL',
                                'th2','KIND_POTENTIAL_TEMPERATURE','UPDATE','NULL','NULL',
                                'q2','KIND_SPECIFIC_HUMIDITY','UPDATE','NULL','NULL',
                                '''
                   
    if update_psfc:
        default_vars += "'psfc','KIND_PRESSURE','UPDATE','NULL','NULL',\n"
    if flag_compute_tendency:
        default_vars +=  "'alt_tend','KIND_ALTIMETER_TENDENCY','UPDATE','NULL','NULL',\n"
    if update_tsk:
        default_vars += "'tsk','KIND_SKIN_TEMPERATURE','UPDATE','NULL','NULL',\n"
    if update_lsm:
        # Updates Soil Temperature (TSLB)
        default_vars += "                                 'tslb1','KIND_SOIL_TEMPERATURE', 'UPDATE','NULL','NULL',\n"
        default_vars += "                                 'tslb2','KIND_SOIL_TEMPERATURE', 'UPDATE','NULL','NULL',\n"
        default_vars += "                                 'tslb3','KIND_SOIL_TEMPERATURE', 'UPDATE','NULL','NULL',\n"
        default_vars += "                                 'tslb4','KIND_SOIL_TEMPERATURE', 'UPDATE','NULL','NULL',\n"
        default_vars += "                                 'tslb5','KIND_SOIL_TEMPERATURE', 'UPDATE','NULL','NULL',\n"
    # Update appropriate number of moisture variables in state vector
    if num_moist_vars >= 3:
        if not minimal_vars:
            default_vars +=  """                                'qv','KIND_VAPOR_MIXING_RATIO','UPDATE','NULL','NULL',
                                'qc','KIND_CLOUD_LIQUID_WATER','UPDATE','NULL','NULL',
                                'qr','KIND_RAINWATER_MIXING_RATIO','UPDATE','NULL','NULL'"""

        if num_moist_vars == 3:
            default_vars +="""\n"""
        else:
            default_vars +=""",\n"""

    if num_moist_vars >= 5:
        if not minimal_vars:
            default_vars +=  """                                'qi','KIND_CLOUD_ICE','UPDATE','NULL','NULL',
                                'qs','KIND_SNOW_MIXING_RATIO','UPDATE','NULL','NULL'"""
        if num_moist_vars == 5: 
            default_vars +="""\n"""
        else:
            default_vars +=""",\n"""



    if num_moist_vars >= 6:
        if not minimal_vars:
            default_vars +=  """                                'qg','KIND_GRAUPEL_MIXING_RATIO','UPDATE','NULL','NULL'"""

    return [default_vars]



def format_true_false(input):
    """ Quick function to write a Fortran friendly ".true." or ".false."
    based on the Python True or False given """
    if input:
        return '.true.'
    else:
        return '.false.'


if __name__ == '__main__':
    main()
