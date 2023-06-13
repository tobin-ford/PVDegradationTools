"""Collection of functions for PV module design considertations.
"""

import numpy as np
import pandas as pd
from numba import jit
from rex import NSRDBX
from rex import Outputs
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from . import humidity


def edge_seal_ingress_rate(avg_psat):
        """
        This function generates a constant k, relating the average moisture ingress rate through a
        specific edge seal, Helioseal 101. Is an emperical estimation the rate of water ingress of
        water through edge seal material. This function was determined from numerical calculations
        from several locations and thus produces typical responses. This simplification works
        because the environmental temperature is not as important as local water vapor pressure.
        For the same environmental water concentration, a higher temperature results in lower
        absorption in the edge seal but lower diffusivity through the edge seal. In practice, these
        effects nearly cancel out makeing absolute humidity the primary parameter determining
        moisture ingress through edge seals.
        
        See: Kempe, Nobles, Postak Calderon,"Moisture ingress prediction in polyisobutylene‐based
        edge seal with molecular sieve desiccant", Progress in Photovoltaics, DOI: 10.1002/pip.2947

        Parameters
        -----------
        avg_psat : float
            Time averaged time averaged saturation point for an environment in kPa. 
            When looking at outdoor data, one should average over 1 year

        Returns
        -------
        k : float [cm/h^0.5]
            Ingress rate of water through edge seal.
            Specifically it is the ratio of the breakthrough distance X/t^0.5.
            With this constant, one can determine an approximate estimate of the ingress distance
            for a particular climate without more complicated numerical methods and detailed
            environmental analysis.

        """

        k = .0013 * (avg_psat)**.4933

        return k

def edge_seal_width(weather_df, meta,
                    k=None,
                    years=25):
    """
    Determine the width of edge seal required for given number of years water ingress.

    Parameters
    ----------
    weather_df : pd.DataFrame
        must be datetime indexed and contain at least temp_air
    meta : dict
        location meta-data (from weather file)
    k: float
        Ingress rate of water through edge seal. [cm/h^0.5]
        Specifically it is the ratio of the breakthrough distance X/t^0.5.
        See the function design.edge_seal_ingress_rate()
    years : integer, default = 25
        Integer number of years under water ingress
    Returns
    ----------
    width : float 
        Width of edge seal required for input number of years water ingress. [cm]
    """

    if k is None:
         avg_psat = humidity.psat(weather_df['temp_air'])[1]
         k = edge_seal_ingress_rate(avg_psat)
    
    width = k * (years * 365.25 * 24)**.5

    return width

#TODO: Where is dew_pt_temp coming from?
def edge_seal_from_dew_pt(weather_df, meta,
                          dew_pt_temp=None,
                          years=25,
                          full_results=False):
    """
    Compute the edge seal width required for 25 year water ingress directly from
    dew pt tempterature.

    Parameters
    ----------
    weather_df : pd.DataFrame
        must be datetime indexed and contain at least 'temp_air' and 'Dew Point'
    meta : dict
        location meta-data (from weather file)
    dew_pt_temp : float, or float series
        Dew Point Temperature [C]
    years : int, optional
        Number of years for water ingress. Default = 25
    full_results : boolean
        If true, returns all calculation steps: psat, avg_psat, k, edge seal width
        If false, returns only edge seal width

    Returns
    ----------
    edge_seal_width: float
        Width of edge seal [mm] required for 25 year water ingress

    Optional Returns
    ----------
    psat : series
        Hourly saturation point
    avg_psat : float
        Average saturation point over sample times
    k : float
        Ingress rate of water vapor
    """
    
    if dew_pt_temp is None:
         dew_pt_temp = weather_df['Dew Point']

    psat, avg_psat = humidity.psat(dew_pt_temp)

    k = .0013 * (avg_psat)**.4933

    width = edge_seal_width(weather_df, meta, k, years)

    res = {'psat':psat,
           'avg_psat':avg_psat,
           'k':k,
           'edge_seal_width':width}

    if full_results:
        return res
    else:
        return width


#TODO: Include gaps functionality