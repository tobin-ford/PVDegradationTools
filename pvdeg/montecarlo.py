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

"""corrlation class
created to be able to access the name attributes to name rows and columns in matrix
seems like the easiest way to access the string data for naming that I need at runtime
"""
class Corr:
    """modeling constants"""
    mc_1 = ''
    mc_2 = ''
    correlation = 0

    def __init__(self, mc_1_string, mc_2_string, corr): # -> None: // what does this -> None do
        self.mc_1 = mc_1_string
        self.mc_2 = mc_2_string
        self.correlation = corr
    
    def getModelingConstants(self)->list:
        return [self.mc_1, self.mc_2]

# size = # correlation coeff
def _symettric_correlation_matrix(corr: list[Corr])->pd.DataFrame:

    # unpack individual modeling constants
    modeling_constants = [mc for i in corr for mc in i.getModelingConstants()]

    uniques = np.unique(modeling_constants)

    # setting up identity matrix, correct (symmetric) labels for columns and rows
    identity_matrix = np.eye(len(corr))
    identity_df = pd.DataFrame(identity_matrix, columns = uniques, index=uniques)

    # still need to walk the matrix
    # -> fill in values
    # make this standalone function if bigger function cannot handle @njit
    for i in range(len(corr)): # because we want to start on the second row
       for j in range(len(corr)): # loops columns

            # diagonal entry case -> skips to start of next row
            if identity_df.iat[i, j] == 1:
                break

            # entries under diagonal {lower triangular - 1}
            else:
                # gets index and column name to check against coeff
                [x, y] = [identity_df.index[i], identity_df.index[j]] # SOMETHING WEIRD IS HAPPENING HERE

                # checks each correlation coefficients attributes to see if it matches the one we want to fill in at the given index
                for relation in corr:  
                    # skip to next correlation coefficient in list
                    if [x, y] != relation.getModelingConstants():
                        pass
                    # fills in appropriate value
                    else:
                        identity_df.iat[i, j] = relation.correlation
                
            
    # mirror the matrix
    # // skip this for now


    # identity_df should be renamed more appropriately 
    return identity_df


corr_Ea_X = Corr('Ea', 'X', 0.0269)
corr_Ea_LnR0 = Corr('Ea', 'LnR0', -0.9995)
corr_X_LnR0 = Corr('X', 'LnR0', -0.0400)

corr_coeff = [corr_Ea_X, corr_Ea_LnR0, corr_X_LnR0]

setup_matrix = _symettric_correlation_matrix(corr_coeff)
print(setup_matrix)