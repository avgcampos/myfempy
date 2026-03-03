import numpy as np

from myfempy.core.mesh.mesh import Mesh


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


class LegacyLine2(Mesh):
    """Mesh Line Class <ConcreteClassService>"""

    def getElementConection(set_mesh):
        """get a linear 2 nodes mesh

        (i)------{1}------(j)

        """

        nel = set_mesh["NX"]

        conec = np.zeros((nel, 3), dtype=np.int64)
        conec[:, 0] = np.arange(1, nel + 1, 1).tolist()
        conec[:, 1] = np.arange(1, nel + 1, 1).tolist()
        conec[:, 2] = np.arange(2, nel + 2, 1).tolist()
        return conec

    def getNodesCoord(set_mesh):
        nel = set_mesh["NX"]
        lx = set_mesh["LX"]

        step = lx / nel

        coord = np.zeros((nel + 1, 4), dtype=np.float64)
        coord[:, 0] = np.arange(1, nel + 2, 1).tolist()
        coord[:, 1] = np.arange(0, lx + step, step).tolist()
        return coord

    def getElementList(conec, meshset, modeldata):
        elemlist = [[None] * 3]
        for ee in range(len(conec)):
            elemlist.append(
                [
                    int(conec[ee, 0]),
                    meshset,
                    modeldata["MATERIAL"]["PROPMAT"][0]["NAME"],
                    modeldata["GEOMETRY"]["PROPGEO"][0]["NAME"],
                    conec[ee, 1:].astype(int).tolist(),
                ]
            )
        elemlist = elemlist[1::][::]
        return elemlist
