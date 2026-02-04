import numpy as np

from myfempy.core.geometry.geometry import Geometry


class UserDefined(Geometry):
    """User Defined Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "userdefined", "idgeo": 99}
        return geoset

    def getCGCoord(tabgeo, inci, element_number):

        cg = {
            "y_max": tabgeo[int(inci[element_number, 3] - 1)]["YMAX"],
            "y_min": tabgeo[int(inci[element_number, 3] - 1)]["YMIN"],
            "z_max": tabgeo[int(inci[element_number, 3] - 1)]["ZMAX"],
            "z_min": tabgeo[int(inci[element_number, 3] - 1)]["ZMIN"],
            "r_max": tabgeo[int(inci[element_number, 3] - 1)]["RMAX"],
        }
        return cg
