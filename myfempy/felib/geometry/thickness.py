import numpy as np

from myfempy.felib.geometry.geometry import Geometry


class Thickness(Geometry):
    '''Thickness Geoemtry Class <ConcreteClassService>'''

    def GeometrySet():
        geoset = {
            'geo': "thickness",
            'idgeo': 90
        }
        return geoset
    
    def getSectionProp(dim_sec):
        
        b = dim_sec["b"]
        h = dim_sec["h"]
        t = dim_sec["t"]
        d = dim_sec["d"]
        
        A = 0
        Izz = 0
        Iyy = 0
        Jxx = 0

        sect_prop = {
            "areacs": A,
            "inerzz": Izz,
            "ineryy": Iyy,
            "inerxx": Jxx,
            "thickn": t,
        }
        return sect_prop
    
    def getCGCoord(tabgeo, inci, element_number):
        b = tabgeo[int(inci[element_number, 3]) - 1, 5]
        h = tabgeo[int(inci[element_number, 3]) - 1, 6]
        t = tabgeo[int(inci[element_number, 3]) - 1, 7]
        d = tabgeo[int(inci[element_number, 3]) - 1, 8]
        
        y_max = 0
        y_min = 0
        z_max = t/2
        z_min = -t/2
        r_max = 0

        cg = {
            'y_max': y_max,
            'y_min': y_min,
            'z_max': z_max,
            'z_min': z_min,
            'r_max': r_max,
        }

        return cg