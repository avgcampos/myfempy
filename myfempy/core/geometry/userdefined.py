import numpy as np

from myfempy.core.geometry.geometry import Geometry


class UserDefined(Geometry):
    """User Defined Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset =  {"geo": "userdefined", "idgeo": 99}
        return geoset
    
    def getCGCoord(tabgeo, inci, element_number):

        cg = {
            "y_max": 1.0,
            "y_min": -1.0,
            "z_max": 1.0,
            "z_min": -1.0,
            "r_max": 1.0,
        }
        return cg
