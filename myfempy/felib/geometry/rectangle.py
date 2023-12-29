import numpy as np
from geometry import Geometry


class Rectangle(Geometry):
    '''Rectangle Geoemtry Class <ConcreteClassService>'''

    def GeometrySet():
        geoset = {
            'geo': ["rectangle", 10],
        }
        return geoset
    
    def getSectionProp(dim_sec):
        
        b = dim_sec["b"]
        h = dim_sec["h"]
        t = dim_sec["t"]
        d = dim_sec["d"]
        
        A = b * h
        Izz = (1 / 12) * b * h**3
        Iyy = (1 / 12) * h * b**3
        Jxx = Iyy + Izz

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
        
        y_max = h * 0.5
        y_min = -h * 0.5
        z_max = b * 0.5
        z_min = -b * 0.5
        r_max = y_max

        cg = {
            'y_max': y_max,
            'y_min': y_min,
            'z_max': z_max,
            'z_min': z_min,
            'r_max': r_max,
        }

        return cg