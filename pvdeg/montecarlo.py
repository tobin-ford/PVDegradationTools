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

'''
def cholesky_correlated(weather_df, meta,
                        
)'''