"""
Collection of classes and functions for standard development.
"""

import numpy as np
import pandas as pd
import pvlib
from rex import NSRDBX
from rex import Outputs
from pathlib import Path
from random import random
from concurrent.futures import ProcessPoolExecutor, as_completed
#from gaps import ProjectPoints

from . import temperature
from . import spectral
from . import utilities
from . import weather

def eff_gap(T_0, T_inf, level=1, x_0=6.1):
    '''
    Calculate a minimum installation distance for rooftop  mounded PV systems.
    
    Parameters
    ----------
    level : int, optional
        Options 1, or 2. Temperature level 1 or level 2 according to IEC TS 63216.
    x0 : float, optional
        Thermal decay constant (cm), [Kempe, PVSC Proceedings 2023]

    Returns
    -------
    x : float
        Recommended installation distance in centimeter per IEC TS 63126.
        Effective gap "x" for the lower limit for Level 1 or Level 0 modules (IEC TS 63216)

    References
    ----------
    M. Kempe, et al. Close Roof Mounted System Temperature Estimation for Compliance 
    to IEC TS 63126, PVSC Proceedings 2023
    '''

    if level == 1:
        T98 = 70
    if level == 2:
        T98 = 80

    T98_0 = T_0.quantile(q=0.98, interpolation='linear')
    T98_inf = T_inf.quantile(q=0.98, interpolation='linear')

    x = -x_0 * np.log(1-(T98_0-T98)/(T98_0-T98_inf))

    return x, T98_0, T98_inf


def calc_standoff(
    weather_df,
    meta,
    tilt=None,
    azimuth=180,
    sky_model='isotropic',
    temp_model='sapm',
    module_type='glass_polymer', # self.module
    level=1,
    x_0=6.1,
    wind_speed_factor=1):
    '''
    Calculate a minimum standoff distance for roof mounded PV systems.
    
    Parameters
    ----------
    weather_df : pd.DataFrame
        Weather data for a single location.
    meta : pd.DataFrame
        Meta data for a single location.
    tilt : float, optional
        Tilt angle of PV system relative to horizontal.
    azimuth : float, optional
        Azimuth angle of PV system relative to north.
    sky_model : str, optional
        Options: 'isotropic', 'klucher', 'haydavies', 'reindl', 'king', 'perez'.
    temp_model : str, optional
        Options: 'sapm', 'pvsyst', 'faiman', 'sandia'.
    module_type : str, optional
        Options: 'glass_polymer', 'glass_glass'.
    level : int, optional
        Options 1, or 2. Temperature level 1 or level 2 according to IEC TS 63216.
    x0 : float, optional
        Thermal decay constant (cm), [Kempe, PVSC Proceedings 2023]
    wind_speed_factor : float, optional
        Wind speed correction factor to account for different wind speed measurement heights 
        between weather database (e.g. NSRDB) and the tempeature model (e.g. SAPM)
    Returns
    -------
    x : float
        Minimum installation distance in centimeter per IEC TS 63126.
        Effective gap "x" for the lower limit for Level 1 or Level 0 modules (IEC TS 63216)

    References
    ----------
    M. Kempe, et al. Close Roof Mounted System Temperature Estimation for Compliance 
    to IEC TS 63126, PVSC Proceedings 2023
    '''
    if module_type == 'glass_polymer':
        conf_0 = 'insulated_back_glass_polymer'
        conf_inf = 'open_rack_glass_polymer'
    elif module_type == 'glass_glass':
        conf_0 = 'close_mount_glass_glass'
        conf_inf = 'open_rack_glass_glass'
  
    solar_position = spectral.solar_position(weather_df, meta)
    poa = spectral.poa_irradiance(weather_df, meta, solar_position, tilt, azimuth, sky_model)
    T_0 = temperature.module(weather_df, poa, temp_model, conf_0, wind_speed_factor)
    T_inf = temperature.module(weather_df, poa, temp_model, conf_inf, wind_speed_factor)
    x, T98_0, T98_inf = eff_gap(T_0, T_inf, level, x_0)

    return {'x':x, 'T98_0':T98_0, 'T98_inf':T98_inf}


def run_calc_standoff(
    project_points, 
    out_dir, 
    tag,
    weather_db,
    weather_satellite,
    weather_names,
    max_workers=None,
    tilt=None,
    azimuth=180,
    sky_model='isotropic',
    temp_model='sapm',
    module_type='glass_polymer',
    level=1,
    x_0=6.1,
    wind_speed_factor=1
):

    """
    parallelization utilizing gaps     #TODO: write docstring
    """

    #inputs
    weather_arg = {}
    weather_arg['satellite'] = weather_satellite
    weather_arg['names'] = weather_names
    weather_arg['NREL_HPC'] = True  #TODO: add argument or auto detect
    weather_arg['attributes'] = [
        'air_temperature', 
        'wind_speed', 
        'dhi', 'ghi', 
        'dni','relative_humidity'
        ]

    all_fields = ['x', 'T98_0', 'T98_inf']

    out_fp = Path(out_dir) / f"out_standoff{tag}.h5"
    shapes = {n : (len(project_points), ) for n in all_fields}
    attrs = {'x' : {'units': 'cm'},
             'T98_0' : {'units': 'Celsius'},
             'T98_inf' : {'units': 'Celsius'}}
    chunks = {n : None for n in all_fields}
    dtypes = {n : "float32" for n in all_fields}

    #TODO: is there a better way to add the meta data?
    nsrdb_fnames, hsds  = weather.get_NSRDB_fnames(
        weather_arg['satellite'], 
        weather_arg['names'], 
        weather_arg['NREL_HPC'])
    
    with NSRDBX(nsrdb_fnames[0], hsds=hsds) as f:
        meta = f.meta[f.meta.index.isin(project_points.gids)]
 
    Outputs.init_h5(
        out_fp,
        all_fields,
        shapes,
        attrs,
        chunks,
        dtypes,
        meta=meta.reset_index()
    )

    future_to_point = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for point in project_points:
            gid = int(point.gid)
            weather_df, meta = weather.load(
                database = weather_db, 
                id = gid, 
                **weather_arg)
            future = executor.submit(
                calc_standoff,
                weather_df, 
                meta,
                tilt, 
                azimuth, 
                sky_model,
                temp_model, 
                module_type, 
                level, 
                x_0,
                wind_speed_factor
            )
            future_to_point[future] = gid

        with Outputs(out_fp, mode="a") as out:
            for future in as_completed(future_to_point):
                result = future.result()
                gid = future_to_point.pop(future)

                ind = project_points.index(gid)
                for dset, data in result.items():
                    out[dset,  ind] = np.array([data])

    return out_fp.as_posix()
