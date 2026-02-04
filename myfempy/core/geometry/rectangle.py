import numpy as np

from myfempy.core.geometry.geometry import Geometry


class Rectangle(Geometry):
    """Rectangle ("Solid") Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "rectangle", "idgeo": 10}
        return geoset

    def getSectionProp(dim_sec):
        b = dim_sec["b"]
        h = dim_sec["h"]

        A = b * h
        Izz = 0.08333333333333 * b * h**3
        Iyy = 0.08333333333333 * h * b**3
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

        y_max = h * 0.5
        y_min = -h * 0.5
        z_max = b * 0.5
        z_min = -b * 0.5
        r_max = 0.2887 * b
                                                                          
        cg = {
            "y_max": y_max,
            "y_min": y_min,
            "z_max": z_max,
            "z_min": z_min,
            "r_max": r_max,
        }

        return cg
