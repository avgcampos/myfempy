import numpy as np

from myfempy.core.geometry.geometry import Geometry


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


class CSection(Geometry):
    """C ("U/ Channel") Section Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "isection", "idgeo": 50}
        return geoset

    def getSectionProp(dim_sec):
        b = dim_sec["b"]
        h = dim_sec["h"]
        t = dim_sec["t"]
        d = dim_sec["d"]

        A = t*h + 2*d*(b-t)
        Izz = 1/12 * b*h**3 - 1/12 * (b-t)*(h-2*d)**3
        Iyy = 1/3 * (h-2*d)*t**3 + 2/3 * d*b**3 - t*h + 2*d*(b-t) * ((1/(t*h + 2*d*(b-t))) * (((h - 2*d)*t**3)/2 + d*b**2))**2
        Jxx = Iyy + Izz

        sect_prop = {
            "areacs": A,
            "inerzz": Izz,
            "ineryy": Iyy,
            "inerxx": Jxx,
            "thickn": 0.0,
        }
        return sect_prop

    def getCGCoord(tabgeo, inci, element_number):
        b = tabgeo[int(inci[element_number, 3] - 1)]["B"]
        h = tabgeo[int(inci[element_number, 3] - 1)]["H"]
        d = tabgeo[int(inci[element_number, 3] - 1)]["D"]
        t = tabgeo[int(inci[element_number, 3] - 1)]["T"]
        AREA = tabgeo[int(inci[element_number, 3] - 1)]["AREACS"]
        IZZ = tabgeo[int(inci[element_number, 3] - 1)]["INERZZ"]


        y_max = 0.5 * h
        y_min = -0.5 * h
        z_max = b - (1/AREA) * (((h - 2*d)*t**3)/2 + d*b**2)
        z_min = -(1/AREA) * (((h - 2*d)*t**3)/2 + d*b**2)
        r_max = np.sqrt(IZZ/AREA)

        cg = {
            "y_max": y_max,
            "y_min": y_min,
            "z_max": z_max,
            "z_min": z_min,
            "r_max": r_max,
        }

        return cg
