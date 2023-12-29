import numpy as np
import scipy.sparse.linalg as spla

# import jax.scipy.sparse.linalg as jaxspla

def linsolve_direct(A, b):
    """Scipy sparse linear solver"""
    x = spla.spsolve(A, b)
    return x

def linsolve_gmres(A, b, x0, M):
    """Use Generalized Minimal RESidual iteration to solve Ax = b"""
    x, info =  spla.gmres(A=A, b=b, x0=x0, tol=1E-10, M=M, maxiter=1000)
    return x, info

# def linsolve_Jaxgmres(A, b, x0, tol, M):
#     """Use Generalized Minimal RESidual iteration to solve Ax = b"""
#     x, info =  jaxspla.gmres(A=A, b=b, x0=x0, tol=tol, M=M)
#     return x, info

def eigsolve_eigsh(A, M, k):
    """Find k eigenvalues and eigenvectors of the real symmetric square matrix or complex Hermitian matrix A"""
    w, v = spla.eigsh(A, k, M, sigma=1, which="LM",)
    return w, v
    
