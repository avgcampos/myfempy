import numpy as np

from myfempy.core.geometry.geometry import Geometry


class RectangleTube(Geometry):
    """Rectangle Tube Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "rectangle_tube", "idgeo": 11}
        return geoset

    def getSectionProp(dim_sec):
        b = dim_sec["b"]
        h = dim_sec["h"]
        t = dim_sec["t"]

        A = b * h - ((b - 2 * t) * (h - 2 * t))
        Izz = (1 / 12) * (b * h**3) - (1 / 12) * ((b - 2 * t) * (h - 2 * t) ** 3)
        Iyy = (1 / 12) * (h * b**3) - (1 / 12) * ((h - 2 * t) * (b - 2 * t) ** 3)
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
        r_max = y_max

        cg = {
            "y_max": y_max,
            "y_min": y_min,
            "z_max": z_max,
            "z_min": z_min,
            "r_max": r_max,
        }

        return cg
