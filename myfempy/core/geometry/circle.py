import numpy as np

from myfempy.core.geometry.geometry import Geometry


class Circle(Geometry):
    """Circle ("Solid") Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "circle", "idgeo": 20}
        return geoset

    def getSectionProp(dim_sec):
        d = dim_sec["d"]

        A = 0.25 * np.pi * d**2
        Izz = 0.015625 * np.pi * d**4
        Iyy = Izz
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
        d = tabgeo[int(inci[element_number, 3] - 1)]["D"]

        y_max = d * 0.5
        y_min = -d * 0.5
        z_max = d * 0.5
        z_min = -d * 0.5
        r_max = d * 0.5

        cg = {
            "y_max": y_max,
            "y_min": y_min,
            "z_max": z_max,
            "z_min": z_min,
            "r_max": r_max,
        }

        return cg
