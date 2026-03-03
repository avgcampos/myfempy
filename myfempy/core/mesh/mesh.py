import os
from abc import ABC, abstractmethod

import numpy as np


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


class Mesh(ABC):
    """
    Element API Class <ClassService>
    """
    @abstractmethod
    def getElementConection():
        pass

    @abstractmethod
    def getNodesCoord():
        pass

    @abstractmethod
    def getElementList():
        pass


class MeshMANUAL(Mesh):
    """Mesh Add Class <ConcreteClassService>"""

    def getNodesCoord(set_mesh):
        return np.asarray(set_mesh["COORD"])

    def getElementConection(set_mesh):
        return np.asarray(set_mesh["INCI"])

    def getElementList(conec, meshset, modeldata):
        elemlist = [[None] * 3]
        for ee in range(len(conec)):
            elemlist.append(
                [
                    int(conec[ee, 0]),
                    meshset,
                    modeldata["MATERIAL"]["PROPMAT"][conec[ee, 1] - 1]["NAME"],
                    modeldata["GEOMETRY"]["PROPGEO"][conec[ee, 2] - 1]["NAME"],
                    conec[ee, 3:].astype(int).tolist(),
                ]
            )
        elemlist = elemlist[1::][::]
        return elemlist
