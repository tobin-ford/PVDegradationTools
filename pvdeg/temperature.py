"""
Collection of classes and functions to calculate different temperatures.
"""

import pvlib

def module(
    weather_df, 
    poa,
    temp_model='sapm', 
    conf='open_rack_glass_polymer',
    wind_speed_factor=1):

    """
    Calculate module temperature based on weather data from the National Solar Radiation 
    Database (NSRDB) for a given location (gid).
    
    Parameters
    ----------
    nsrdb_file : str
        The file path to the NSRDB file.
    gid : int
        The geographical location ID in the NSRDB file.
    poa : pandas.DataFrame
         Contains keys/columns 'poa_global', 'poa_direct', 'poa_diffuse', 
         'poa_sky_diffuse', 'poa_ground_diffuse'.
    temp_model : str, optional
        The temperature model to use, 'sapm' from pvlib by default.
    conf : str, optional
        The configuration of the PV module architecture and mounting configuration.
        Options: 'open_rack_glass_polymer' (default), 'open_rack_glass_glass',
        'close_mount_glass_glass', 'insulated_back_glass_polymer'
    
    Returns
    -------
    module_temperature : pandas.DataFrame
        The module temperature in degrees Celsius at each time step.
    """
    parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS[temp_model][conf]
    module_temperature = pvlib.temperature.sapm_module(
        poa_global=poa['poa_global'], 
        temp_air=weather_df['air_temperature'], 
        wind_speed=weather_df['wind_speed']*wind_speed_factor,
        a=parameters['a'],
        b=parameters['b'])

    return module_temperature