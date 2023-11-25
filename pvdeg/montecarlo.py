"""Colection of functions for monte carlo simulations.
"""

import numpy as np
import pandas as pd
from numba import njit
import pvlib # is this ok to include
from scipy.linalg import cholesky
from scipy import stats

from . import spectral
from . import temperature

# Quasi-Monte Carlo could be used to generate samples matrix
# best way of doing this might be to have a function for each of the equations in the degredation database, we dont need every solution

class ArugmentError(Exception):
    pass

def _symettric_correlation_matrix(*args):
    valid_corr = ['corr_Ea_X', 'corr_Ea_LnR0', 'corr_X_LnR0']
    
    for arg in args:
        if arg not in valid_corr:
            raise ArugmentError(f"invalid argument: {arg}. valid arguments are {valid_corr}") 
        if (arg >= -1 and arg <= 1):
            raise ValueError(f"invalid value: {arg}. valid arguments include [-1, 1]")

    identity_matrix = np.eye(len(args))
    
    # fill in lower triangular portion of matrix
    # following sequence of triangular numbers as follows
    # 0, 1, 3, 6, 10, 15
    
    # walks the matrix
    i = 1
    j = 0
    count = 0
    while i < len(args):
        if identity_matrix[i, j] != 1:
            identity_matrix[i, j] = args[count]
            count += 1
            j += 1
        else:
            i += 1
            j = 0

    # rename identity_matrix to accurate 
    symmetric_matrix = np.triu(identity_matrix) + np.triu(identity_matrix, k = 1).T
    return symmetric_matrix


def correlated(weather_df, meta,
               function,
               MC_params
               ):
    """
    Parameters
    ----------
    weather_df : pd.dataframe
        Dataframe containing at least dni, dhi, ghi, temperature, wind_speed
    meta : dict
        Location meta-data containing at least latitude, longitude, altitude
    MC_params : dict
        Dictionary with key-value pair of string to indicate modeling constant,
        value of list containing appropritate mean, standard deviation as floats
        could look like  {'ea' : [mean, stdev] }
    """