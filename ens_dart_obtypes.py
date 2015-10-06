#!/usr/bin/python

# List of which ob types to use in DART assimilation
# NEW FORMAT -- Luke Madaus (7/5/2012)
#  -- Use Python's internal "True" and "False" values
#  -- After each sub-ob-type, there is a tuple of two
#     values--the first is the option to ASSIMILATE this
#     type, the second is the option to EVALUATE this type.
#     e.g.  To evaulate but not assimilate sfc_altimeter:
#     sfc_altimeter = (False, True)
#     The make_namelist_dart script is designed to override
#     the evaluate choice if the ob is being the assimilated.
#     Therefore, a value of (True, True) is the same as a value
#     of (True, False) as the ob HAS to be evaluated to be assimilated.

# FORMAT  --  sfc_wind = (ASSIMILATE?, EVALUATE?)

### SOUNDINGS ###
use_obs_soundings = False 

rawin_sfc_pres      = (False, False) 
rawin_sfc_altimeter = (False, False) 
rawin_wind          = (False, False) 
rawin_temp          = (False, False) 
rawin_humidity      = (False, False) 


### SURFACE OBS ###
use_obs_surface = True

sfc_altimeter    = (False, False) 
sfc_alt_tendency = (False, False)
sfc_verify_alt   = (False, False)
sfc_pressure     = (False, False) 
sfc_wind         = (False, False) 
sfc_temp         = (False, False)
sfc_humidity     = (False, False) 
ideal_sfc_pressure = (True, False) 
ideal_sfc_wind     = (True, False) 
ideal_sfc_temp     = (True, False) 
ideal_sfc_humidity = (True, False) 
wunderground_alt = (False, False)
wunderground_temp = (False, False)
wunderground_wind = (False, False)
wunderground_dewp = (False, False)
phone_altimeter  = (False, False)



### Odd PACNW option ###
use_obs_pacnw=False

pacnw_wind      = (False, False) 
pacnw_altimeter = (False, False)
pacnw_temp      = (False, False) 
pacnw_humidity  = (False, False) 

### AIRCRAFT ###
use_obs_aircraft = False 

acars_wind     = (False, False) 
acars_temp     = (False, False) 
acars_humidity = (False, False) 
tamdar_wind     = (False, False) 
tamdar_temp     = (False, False) 
tamdar_humidity = (False, False) 


### CLOUD-TOP WINDS ###
use_obs_ctw = False 

sat_ctw = (False, False)

### GPS REFRACTIVITY ###
use_obs_gpsro = False 

gpsro = (False, False)


### NEXRAD RADAR DATA ###
use_nexrad_data = False

radar_vad_winds          = (False, False)
radar_reflectivity       = (False, False)
radar_radial_velocity    = (False, False)
radar_clear_reflectivity = (False, False)
radar_drop_fall_speed    = (False, False)


### DROPSONDES ###
use_obs_dropsondes = False 

drop_pres     = (False, False)
drop_wind     = (False, False) 
drop_temp     = (False, False) 
drop_humidity = (False, False) 


### DRIFTSONDES ###
use_obs_driftsondes = False 

drift_wind     = (False, False) 
drift_temp     = (False, False) 
drift_humidity = (False, False) 


### TC INFO ###
use_obs_tcinfo = False

vortex_loc  = (False, False) 
vortex_pmin = (False, False) 
vortex_wmax = (False, False) 

