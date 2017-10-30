# Numerical schemes for simulating diffusion for outer code diffusion.py

from __future__ import absolute_import, division, print_function
import numpy as np

# The linear algebra package for BTCS (for solving the matrix equation)
import scipy.linalg as la

def FTCS(phiOld, d, nt):
    """
    Diffusion of profile in phiOld using FTCS using FTCS using non-dimensional
    diffusion coeficient, d
    :param: phiOld: the inital values of phi
    :param: d: given diffusion coefficeient
    :param: nt: number of time intervals
    """
    
    nx = len(phiOld)
    # new time-step array for phi
    phi = phiOld.copy()
    # for each time iteration
    for t in xrange(int(nt)):
 
        # endpoints
        phi[0] = phi[0] + d*(phi[1] - 2 * phi[0])
        phi[-1] = phi[-1] + d*(-2 * phi[-1] + phi[-2])
        for index in xrange(1, nx-1):
    
             # midpoints
            phi[index] = phi[index] + d*(phi[index + 1] - 2 * phi[index] + phi[index - 1])

    return phi


def BTCS(phi, d, nt):
    """
    Diffusion of profile in phi using BTCS using non-dimensional
    diffusion coefficent, d assuming fixed value boundary consditions
    """
    
    nx = len(phi)
    
    # array representing BTCS
    M = np.zeros([nx,nx])
    # Zero gradient boundary conditions
    M[0,0] = 1.    
    M[0,1] = -1.
    M[-1,-1] = 1.
    M[-1,-2] = -1.
    for i in xrange(1, nx-1):
        M[i,i-1] = -d
        M[i,i] = 1+2*d
        M[i,i+1] = -d
    
    # BTCS for all timesteps
    for it in xrange(int(nt)):
        # RHS for zero gradient boundary conditions
        phi[0] = 0
        phi[-1] = 0
        
        phi = la.solve(M,phi)
    
    return phi
        