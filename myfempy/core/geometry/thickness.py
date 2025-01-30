import numpy as np

from myfempy.core.geometry.geometry import Geometry


class Thickness(Geometry):
    """Thickness Geoemtry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "thickness", "idgeo": 1}
        return geoset

    def getSectionProp(dim_sec):
        t = dim_sec["t"]

        sect_prop = {
            "areacs": 0.0,
            "inerzz": 0.0,
            "ineryy": 0.0,
            "inerxx": 0.0,
            "thickn": t,
        }
        return sect_prop

    def getCGCoord(tabgeo, inci, element_number):
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]

        y_max = 0
        y_min = 0
        z_max = 0.5 * t
        z_min = -0.5 * t
        r_max = 0

        cg = {
            "y_max": y_max,
            "y_min": y_min,
            "z_max": z_max,
            "z_min": z_min,
            "r_max": r_max,
        }

        return cg
