import numpy as np

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.material.material import Material


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


class HeatSolid(Material):
    """Heat Solid Isotropic Material Class <ConcreteClassService>"""

    def getMaterialSet():
        matset = {
            "mat": "heatsolid",
            "type": "isotropic",
        }
        return matset

    def getElasticTensor(tabmat, inci, element_number, Model=None):
        Kxx = tabmat[int(inci[element_number, 2]) - 1]["KXX"]
        Kyy = tabmat[int(inci[element_number, 2]) - 1]["KYY"]
        Kzz = tabmat[int(inci[element_number, 2]) - 1]["KZZ"]
        D = np.zeros((3, 3), dtype=FLT64)
        D[0, 0] = Kxx
        D[1, 1] = Kyy
        D[2, 2] = Kzz
        return D

    def getElementGradTemp(Model, U, ptg, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])

        nodelist = Model.shape.getNodeList(Model.inci, element_number)

        loc = Model.shape.getLocKey(nodelist, nodedof)

        elementcoord = Model.shape.getNodeCoord(Model.coord, nodelist)
        
        diffN = Model.shape.getDiffShapeFuntion(np.array([ptg, ptg, ptg]), nodedof)
        
        invJ = Model.shape.getinvJacobi(np.array([ptg, ptg, ptg]), elementcoord, nodedof)
        
        B = Model.element.getB(diffN, invJ)

        N = Model.shape.getShapeFunctions(np.array([ptg, ptg, ptg]), nodedof)

        epsilon = np.dot(B, U[loc])  # B @ (U[loc])

        epsilon_T = np.dot(N, U[loc])

        strn_elm_xx = epsilon[0]
        
        strn_elm_yy = epsilon[1]
        
        strn_elm_zz = epsilon[2]

        # strn_elm_vm = np.sqrt(epsilon[0]**2 + epsilon[1]**2)

        strain = [epsilon_T[0], strn_elm_xx, strn_elm_yy, strn_elm_zz]

        return epsilon, strain

    def getTitleGradTemp():
        title = ["GRADTEMP", "GRADTEMP_XX", "GRADTEMP_YY", "GRADTEMP_ZZ"]
        return title

    def getElementHeatFlux(Model, epsilon, element_number):

        C = Model.material.getElasticTensor(Model.tabmat, Model.inci,  element_number)

        sigma = -1 * np.dot(C, epsilon)

        strs_elm_xx = sigma[0]
        
        strs_elm_yy = sigma[1]
        
        strs_elm_zz = sigma[2]

        strs_elm_vm = np.sqrt(sigma[0] ** 2 + sigma[1] ** 2 + sigma[2] ** 2)

        stress = [strs_elm_vm, strs_elm_xx, strs_elm_yy, strs_elm_zz]

        return sigma, stress

    def getTitleHeatFlux():
        title = ["HEATFLUX_MAG", "HEATFLUX_XX", "HEATFLUX_YY", "HEATFLUX_ZZ"]
        return title
    
