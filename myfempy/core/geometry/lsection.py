import numpy as np

from myfempy.core.geometry.geometry import Geometry


class LSection(Geometry):
    """L ("Unequal-Legged Angle") Section Geometry Class <ConcreteClassService>"""

    def GeometrySet():
        geoset = {"geo": "isection", "idgeo": 60}
        return geoset

    def getSectionProp(dim_sec):
        b = dim_sec["b"]
        h = dim_sec["h"]
        t = dim_sec["t"]
        # d = dim_sec["d"]

        A = (h + b - t)*t
        Izz = t/3  *(b*t**2 + h**3 - t**3) - ((h + b - t)*t) * (t/(2*((h + b - t)*t)) * (h**2 + b*t -t**2))**2
        Iyy = t/3 * (h*t**2 + b**3 - t**3) - ((h + b - t)*t) * (t/(2*((h + b - t)*t)) * (b**2 + h*t -t**2))**2
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
        # d = tabgeo[int(inci[element_number, 3] - 1)]["D"]
        t = tabgeo[int(inci[element_number, 3] - 1)]["T"]
        AREA = tabgeo[int(inci[element_number, 3] - 1)]["AREACS"]
        IZZ = tabgeo[int(inci[element_number, 3] - 1)]["INERZZ"]


        y_max = h - t/(2*AREA) * (h**2 + b*t -t**2)
        y_min = -t/(2*AREA) * (h**2 + b*t -t**2)
        z_max = b - t/(2*AREA) * (b**2 + h*t -t**2)
        z_min = -t/(2*AREA) * (b**2 + h*t -t**2)
        r_max = np.sqrt(IZZ/AREA)

        cg = {
            "y_max": y_max,
            "y_min": y_min,
            "z_max": z_max,
            "z_min": z_min,
            "r_max": r_max,
        }

        return cg
