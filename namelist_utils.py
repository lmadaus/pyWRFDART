#!/usr/bin/env python
from __future__ import print_function, division
from collections import OrderedDict
from datetime import datetime, timedelta

def read_namelist(infname='namelist.input'):
    """
    This reads in the namelist file given by
    infname and returns a dictionary of
    the namelist contents organized by subsection
    and then variable name
    """
    nmld = OrderedDict()
    # Putting this in a "with" statement will close the
    # file when we are done
    with open(infname, 'r') as infile:
        for line in infile.readlines():
            cleanline = line.strip() # This is like Fortran 'trim'
            # Skip blank lines and section ends
            if cleanline in ['', '/']:
                continue
            # If we start with an ampersand, this is a section header
            if cleanline[0] == '&':
                sectionhead = cleanline[1:]
                # Make a section in the nmld
                nmld[sectionhead] = {}
                continue
            # If we're not blank, not a section end, and
            # we dont begin with '&', this is a variable
            # Split into to variables on equals sign
            try:
                varname, value = cleanline.split('=')
                varname = varname.strip()
                value = value.strip()
            except ValueError:
                # If the above split failed, we are still
                # in the previous variable
                # Hold over varname and append to current value
                if isinstance(nmld[sectionhead][varname], str):
                    oldvalue = nmld[sectionhead][varname]
                else:
                    oldvalue = ','.join([str(x) for x in nmld[sectionhead][varname]])
                value = ','.join([oldvalue,cleanline.strip()])
            # Remove comma at end, if it exists
            value = value.strip(',')
            # If there are other commas in this, it should
            # be stored as a list.  Trying to split
            # on commas will returna list
            valsplit = value.split(',')
            # Convert to correct types
            outvals = str_to_value(valsplit)
            # If this is length 1, don't store
            # the list
            if len(outvals) == 1:
                nmld[sectionhead][varname] = outvals[0]
            else:
                nmld[sectionhead][varname] = outvals
            # This will automatically continue to the next line here
    # Exit the with statement to close the file
    # And return the namelist dictionary
    return nmld
            
def str_to_value(instr):
    """
    Function to determine if the contents
    of a list or single string should stay a
    string, or be converted to a float or integer
    """
    if isinstance(instr, str):
        # Enclose this in a list for our usage
        instr = [strval]
    out = []
    for x in instr:
        try:
            # If this fails, this must be a string
            # and it will jump to the exception below
            float(x)
            # If there's a decimal point, store as float
            if '.' in x:
                out.append(float(x))
            else:
                out.append(int(x))

        except ValueError:
            # Remove any leading/trailing space and quotes
            strval = x.strip().strip("'").strip('"')
            # Check if it's a boolean
            if strval in ['.true.','T']:
                strval = True
            elif strval in ['.false.','T']:
                strval = False
            out.append(strval)
    return out
        

def format_single(v):
    """ 
    This function formats single values
    for the namelist writeout 
    """
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
        # Return the stringvalue with single quotes
        return "'{:s}'".format(v)

def var_format(invar, multiline=None):
    """ 
    Decide how to format this namelist value depending on what it is
    """
    # If it's a list, format each item individually and
    # join with commas
    if isinstance(invar, list):
        if len(invar) == 1:
            return format_single(invar[0])
        values = [format_single(v) for v in invar]
        # Put on multiple lines if longer than a certain amount
        if multiline is not None:
            total_rows = int(len(values) / multiline)
            allrows = []
            for r in xrange(total_rows):
                thisrow = values[r*multiline:(r+1)*multiline]
                rowjoin = ', '.join(thisrow)
                allrows.append(rowjoin)
            return '\n                           '.join(allrows)
        else:
            return ', '.join(values)
        
        return ', '.join(values)
    else:
        # Return the single format
        return format_single(invar)

def write_namelist(nmld, outfname='namelist.input'):
    """ 
    Formats a namelist from a dictionary (nmld)
    and writes it to file given by outfname
    """
    # Open the output file for writing
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
                # Special check here for dart nml
                if section_name == 'model_nml' and varname in ['model_variables']:
                    print(''.join((indent,varname.ljust(24), '= ',
                               var_format(value, multiline=5), ',')), file=outfile)
                else:
                    print(''.join((indent,varname.ljust(24), '= ',
                               var_format(value), ',')), file=outfile)
            # Close the section
            print(' /', file=outfile)
            print('', file=outfile)


def write_dart_namelist(nmld=None, mem=1, date=None, window_minutes=15, outname='input.nml'):
    """ Function that makes a new input.nml file in the current working directory 
        nmld is "required", but following can be specified:
            mem => The current member number, if in the future we want
                    member-specific changes (currently does nothing)
            date => The time to specify in the namelist.  *Most* time-related
                    variables will be overwritten with this time
            window_minutes => For binning observations, how many minutes (+/-date) to
                    consider the same bin
            outname => the name of the output file to write

        RETURNS:
            True if successful
            writes input.nml file
    """
    if nmld is None:
        if 'dir_dom' in locals().keys():
            # Read the master namelist
            nmld = read_namelist('{:s}/input.nml'.format(dir_dom))
        else:
            try:
                nmld = read_namelist('input.nml')
            except:
                print("Unable to read namelist: input.nml")
                return False


    # If these variables are from a parameter file, use them
    if all (x in locals().keys() for x in ['ens_size','cov_cutoff','assim_loc_meth']):
        # Set certain variables based on ens_dart_param.py
        nmld['filter_nml']['ens_size'] = ens_size
        nmld['filter_nml']['num_output_obs_members'] = ens_size
        nmld['filter_nml']['num_output_state_members'] = ens_size
        nmld['assim_tools_nml']['cov_cutoff'] = cov_cutoff
        nmld['cov_cutoff_nml']['select_localization'] = assim_loc_meth
        nmld['closest_member_tool_nml']['ens_size'] = ens_size
        nmld['restart_file_tool_nml']['ens_size'] = ens_size
        nmld['restart_file_utility_nml']['ens_size'] = ens_size
        nmld['covariance_relax_nml']['ens_size'] = ens_size
        nmld['inflate_ens_nml']['ens_size'] = ens_size
   
    # Set up time variables
    if date is not None:
        epoch = datetime(1601,1,1,0)
        darttime = date - epoch
        # Figure out the window
        start_window = date - timedelta(minutes=window_minutes)
        end_window = date + timedelta(minutes=window_minutes)
        first_ob = start_window - epoch
        last_ob = end_window - epoch

        # Date times
        nmld['schedule_nml']['first_bin_start'] = [start_window.year, start_window.month, start_window.day,\
                                                   start_window.hour, start_window.minute, start_window.second]
        nmld['schedule_nml']['first_bin_end'] = [end_window.year, end_window.month, end_window.day,\
                                                   end_window.hour, end_window.minute, end_window.second]
        nmld['schedule_nml']['last_bin_end'] = [end_window.year, end_window.month, end_window.day,\
                                                   end_window.hour, end_window.minute, end_window.second]
        nmld['schedule_nml']['bin_interval_days'] = 0
        nmld['schedule_nml']['bin_interval_seconds'] = (end_window-start_window).seconds
        nmld['obs_sequence_tool_nml']['filename_seq'] = '{:%Y%m%d%H%M%S}_obs_seq.prior'.format(date)
        nmld['obs_sequence_tool_nml']['first_obs_days'] = first_ob.days
        nmld['obs_sequence_tool_nml']['first_obs_seconds'] = first_ob.seconds
        nmld['obs_sequence_tool_nml']['last_obs_days'] = last_ob.days
        nmld['obs_sequence_tool_nml']['last_obs_seconds'] = last_ob.seconds
        nmld['obs_seq_coverage_nml']['first_analysis'] = [start_window.year, start_window.month, start_window.day,\
                                                   start_window.hour, start_window.minute, start_window.second]
        nmld['obs_seq_coverage_nml']['last_analysis'] = [end_window.year, end_window.month, end_window.day,\
                                                   end_window.hour, end_window.minute, end_window.second]




    write_namelist(nmld, outname)


if __name__ == '__main__':
    nmld = read_namelist('input.nml')
    write_dart_namelist(nmld,date=datetime(2014,8,12,12),outname='new_input.nml') 
    #print(nmld['model_nml'])
