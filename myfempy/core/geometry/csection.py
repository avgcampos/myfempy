import numpy as np

from myfempy.core.geometry.geometry import Geometry


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
