from __future__ import annotations

import numpy as np

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import get_nodes_from_list


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


class BoundCondThermal(Structural):
    """Structural Load Class <ConcreteClassService>"""

    def getBCApply(Model, bclist):
        boncdnodeaply = np.zeros((1, 4))

        if bclist["TYPE"] == "insulated":
            bcapp = BoundCondThermal.__BCFixed(Model, bclist)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        elif bclist["TYPE"] == "temperature":
            bcapp = BoundCondThermal.__BCDispl(Model, bclist)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        # elif bclist['TYPE'] == "csymm":
        #     bcapp = BoundCondThermal.__BCCS(modelinfo, bclist)
        #     boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        else:
            pass
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    def __BCFixed(Model, bclist):
        boncdnodeaply = np.zeros((1, 4))

        nodelist = [
            bclist["DIR"],
            bclist["LOCX"],
            bclist["LOCY"],
            bclist["LOCZ"],
            bclist["TAG"],
            bclist["MESHNODE"],
        ]
        node_list_bc, dir_fc = get_nodes_from_list(
            nodelist, Model.coord, Model.regions
        )

        if bclist["DOF"] == "full":
            bcdof = 0
        else:
            bcdof = Model.modelinfo["dofs"]["d"][bclist["DOF"]]

        for j in range(len(node_list_bc)):
            bcapp = np.array([[int(node_list_bc[j]), bcdof, 0.0, int(bclist["STEP"])]])
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    def __BCDispl(Model, bclist):
        boncdnodeaply = np.zeros((1, 4))

        nodelist = [
            bclist["DIR"],
            bclist["LOCX"],
            bclist["LOCY"],
            bclist["LOCZ"],
            bclist["TAG"],
            bclist["MESHNODE"],
        ]
        node_list_bc, dir_fc = get_nodes_from_list(
            nodelist, Model.coord, Model.regions
        )

        bcdof = Model.modelinfo["dofs"]["d"][bclist["DOF"]]

        for j in range(len(node_list_bc)):
            bcapp = np.array(
                [
                    [
                        int(node_list_bc[j]),
                        bcdof,
                        float(bclist["VAL"]),
                        int(bclist["STEP"]),
                    ]
                ]
            )
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    # def __BCCS(modelinfo, bclist):
    #     boncdnodeaply = np.zeros((1, 4))

    #     nodelist = [bclist['DIR'], bclist['LOCX'], bclist['LOCY'], bclist['LOCZ'], bclist['TAG']]
    #     node_list_bc, dir_fc = get_nodes_from_list(
    #         nodelist, modelinfo["coord"], modelinfo["regions"]
    #     )

    #     if bclist['DOF'] == "left":
    #         bcdof = 11
    #     elif bclist['DOF'] == "right":
    #         bcdof = 12
    #     else:
    #         bcdof = 0

    #     for j in range(len(node_list_bc)):
    #         bcapp = np.array([[int(node_list_bc[j]), bcdof, 0.0, int(bclist['STEP'])]])
    #         boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

    #     boncdnodeaply = boncdnodeaply[1::][::]
    #     return boncdnodeaply
