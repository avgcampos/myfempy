# import cython
import numpy as np
from scipy import sparse

__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""


def getVectorization(ith, jth, val, loc, matrix, ee, elemdof):

    ith[(elemdof * elemdof) * ee : (elemdof * elemdof) * (ee + 1)] = np.tile(
        loc.reshape(1, elemdof).T, (1, elemdof)
    ).flatten("F")
    jth[(elemdof * elemdof) * ee : (elemdof * elemdof) * (ee + 1)] = np.transpose(
        np.tile(loc.reshape(1, elemdof).T, (1, elemdof))
    ).flatten("F")
    val[(elemdof * elemdof) * ee : (elemdof * elemdof) * (ee + 1)] = matrix.flatten("F")

    return ith, jth, val


def getLoadAssembler(loadaply, nodetot, nodedof):
    """
    getLoadAssembler Assembler module <ConcreteClassService>

    Returns:
        vector sparse
    """

    steps = len(np.unique(loadaply[:, 3]))
    if steps == 0:
        steps = 1
    else:
        pass

    forcevec = np.zeros((nodedof * nodetot, steps))

    for fstep in range(len(np.unique(loadaply[:, 3]))):
        forceaply = loadaply[np.where(loadaply[:, 3] == fstep + 1), :][0]
        nload = forceaply.shape[0]

        if nodedof == 1:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]
                else:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

        elif nodedof == 2:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 2:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

        elif nodedof == 3:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 2:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 3:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

        elif nodedof == 6:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 2:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 3:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 4:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 3))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 5:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 4))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 6:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 5))
                    forcevec[gdlload, fstep] += forceaply[ii, 2]
        else:
            pass
    # forcevec = sparse.csc_matrix(forcevec)
    return forcevec


# Dirichlet Homogeneous https://en.wikipedia.org/wiki/Dirichlet_boundary_condition
def getConstrains(constrains, nodetot, nodedof):
    """
    getConstrains Constrain module <ConcreteClassService>

    Returns:
        _description_
    """

    ntbc = len(constrains)
    fixedof = np.zeros((nodedof * nodetot, 1), dtype=int)
    constdof = np.zeros((nodedof * nodetot, 1), dtype=int)

    if nodedof == 1:
        for ii in range(ntbc):
            no = int(constrains[ii, 0])
            if int(constrains[ii, 1]) == 1:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 1, 0] = nodedof * no
            else:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 1, 0] = nodedof * no

    elif nodedof == 2:
        for ii in range(ntbc):
            no = int(constrains[ii, 0])
            if int(constrains[ii, 1]) == 1:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                else:
                    constdof[nodedof * no - 2, 0] = nodedof * no - 1
            elif int(constrains[ii, 1]) == 2:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 1, 0] = nodedof * no
            else:  # int(constrains[ii, 1]) == 0:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[0, nodedof * no - 2, 0] = nodedof * no - 1
                    constdof[0, nodedof * no - 1, 0] = nodedof * no

    elif nodedof == 3:
        for ii in range(ntbc):
            no = int(constrains[ii, 0])
            if int(constrains[ii, 1]) == 1:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 3, 0] = nodedof * no - 2
                else:
                    constdof[nodedof * no - 3, 0] = nodedof * no - 2
            elif int(constrains[ii, 1]) == 2:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                else:
                    constdof[nodedof * no - 2, 0] = nodedof * no - 1
            elif int(constrains[ii, 1]) == 3:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 1, 0] = nodedof * no
            else:  # int(constrains[ii, 1]) == 0:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 3, 0] = nodedof * no - 2
                    fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 3, 0] = nodedof * no - 2
                    constdof[nodedof * no - 2, 0] = nodedof * no - 1
                    constdof[nodedof * no - 1, 0] = nodedof * no

    elif nodedof == 6:
        for ii in range(ntbc):
            no = int(constrains[ii, 0])
            if int(constrains[ii, 1]) == 1:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 6, 0] = nodedof * no - 5
                else:
                    constdof[nodedof * no - 6, 0] = nodedof * no - 5

            elif int(constrains[ii, 1]) == 2:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 5, 0] = nodedof * no - 4
                else:
                    constdof[nodedof * no - 5, 0] = nodedof * no - 4

            elif int(constrains[ii, 1]) == 3:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 4, 0] = nodedof * no - 3
                else:
                    constdof[nodedof * no - 4, 0] = nodedof * no - 3

            elif int(constrains[ii, 1]) == 4:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 3, 0] = nodedof * no - 2
                else:
                    constdof[nodedof * no - 3, 0] = nodedof * no - 2

            elif int(constrains[ii, 1]) == 5:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                else:
                    constdof[nodedof * no - 2, 0] = nodedof * no - 1

            elif int(constrains[ii, 1]) == 6:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 1, 0] = nodedof * no

            else:  # int(constrains[ii, 1]) == 0:
                if constrains[ii, 2] == 0.0:
                    fixedof[nodedof * no - 6, 0] = nodedof * no - 5
                    fixedof[nodedof * no - 5, 0] = nodedof * no - 4
                    fixedof[nodedof * no - 4, 0] = nodedof * no - 3
                    fixedof[nodedof * no - 3, 0] = nodedof * no - 2
                    fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                    fixedof[nodedof * no - 1, 0] = nodedof * no
                else:
                    constdof[nodedof * no - 6, 0] = nodedof * no - 5
                    constdof[nodedof * no - 5, 0] = nodedof * no - 4
                    constdof[nodedof * no - 4, 0] = nodedof * no - 3
                    constdof[nodedof * no - 3, 0] = nodedof * no - 2
                    constdof[nodedof * no - 2, 0] = nodedof * no - 1
                    constdof[nodedof * no - 1, 0] = nodedof * no

    fixedof = fixedof[np.nonzero(fixedof)]
    fixedof = fixedof - np.ones_like(fixedof)
    constdof = constdof[np.nonzero(constdof)]
    constdof = constdof - np.ones_like(constdof)
    alldof = np.arange(0, nodedof * nodetot, 1, int)
    freedof = np.setdiff1d(alldof, fixedof)
    freedof = np.setdiff1d(freedof, constdof)
    return freedof, fixedof, constdof


# Dirichlet Non-Homogeneous
def getDirichletNH(constrains, nodetot, nodedof):
    """
    getLoadAssembler Assembler module <ConcreteClassService>

    Returns:
        vector sparse
    """

    steps = len(np.unique(constrains[:, 3][constrains[:, 3] != 0]))
    if steps == 0:
        steps = 1
    else:
        pass
    
    Uc = np.zeros((nodedof * nodetot, steps), dtype=np.float64)  # solution constrains

    for cstep in range(len(np.unique(constrains[:, 3][constrains[:, 3] != 0]))):
        forceaply = constrains[np.where(constrains[:, 3] == cstep + 1), :][0]
        nload = forceaply.shape[0]

        if nodedof == 1:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

        elif nodedof == 2:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 2:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

        elif nodedof == 3:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 2:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 3:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

        elif nodedof == 6:
            for ii in range(nload):
                if int(forceaply[ii, 1]) == 1:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 2:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 3:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 4:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 3))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 5:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 4))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

                elif int(forceaply[ii, 1]) == 6:
                    gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 5))
                    Uc[gdlload, cstep] = forceaply[ii, 2]

        else:
            pass
    return Uc
