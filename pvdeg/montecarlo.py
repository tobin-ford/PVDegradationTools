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

class ArugmentError(Exception):
    pass


def _symettric_correlation_matrix(*args):
    """
    Helper function. creates a symettric correlation matrix from correlation coefficients

    Parameters
    ----------
    ### replace this but I didn't know how to write it better
    ### takes correlation coefficients as arguments = float
    ### ex : corr_X_Y = float
    """

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


# this probably isnt supposed to be *args
# it should probably just be the dictionary as shown in correlated()
def correlated_data_matrix(iterations, *args):
    """
    Helper function. Generates a matrix of samples for all modeling constants
    
    ----------
    iterations : integer
        number of samples to run for monte carlo simulation
    modeling_constants : dict, key : string, value : list[float]
        contains appropriate modeling constants and their respective mean and standard deviation
        ### I think *args should just become the dictionary

    Returns
    ----------
    correlated_df : pd.DataFrame
        returns a tall DataFrame of correlated data to be used in monte carlo trials with appropriate columns named as modeling constants 
    """
    valid_modeling_constants = ['R_0', 'lnR_0', 'R_D', 'E_a', 't_fail', 'A', 'FF', 'B', 'E'] # stopped at LeTID 'v_ab'

    for arg in args:
        if arg not in valid_modeling_constants:
            raise ArugmentError(f"invalid argument: {arg}. valid arguments are {valid_modeling_constants}") 
        if (arg >= -1 and arg <= 1):
            raise ValueError(f"invalid value: {arg}. valid arguments include [-1, 1]")

    uncorrelated_samples_matrix = np.random.normal(loc=0, scale=1, size=(len(args), iterations))
    coefficient_matrix = _symettric_correlation_matrix(args)

    # the decomposition is not in here
    # seems like it needs to be
    # look at jupyter notebook

    # i think if should go right here 

    correlated_samples = np.matmul(coefficient_matrix, uncorrelated_samples_matrix)
    correlated_samples = np.matrix(correlated_samples)

    keys_list = list(args.keys())
    mean_list = [value[0] for value in args.values()]
    stdev_list = [value[1] for value in args.values()]

    stdev_matrix = np.transpose(np.matrix(stdev_list))
    result = np.multiply(correlated_samples, stdev_matrix) + np.transpose(np.matrix([mean_list]))

    correlated_df = pd.DataFrame(np.transpose(result), columns=keys_list)
    return correlated_df



def correlated(weather_df, meta,
               function,
               MC_params, samples, corr
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

    # need to figure out the right way to call these
    A = _symettric_correlation_matrix(corr)
    decomposition = cholesky(A, lower = True)

    # the way I am passing MC_params to the function may be wrong
    correlated_df = correlated_data_matrix(20000, MC_params)