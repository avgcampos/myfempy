#!/usr/bin/env python
__doc__ = """
cross section
"""
import numpy as np


def sec_def(keysecdef: str):
    """cross section def

    Arguments:
        keysecdef:str -- key section def

    Returns:
        idsecdef:int  -- id number of cross section 
    """
    matdef = {
        "rectangle": 10,
        "R": 10,
        "rectangletube": 11,
        "Rt": 11,
        "circle": 20,
        "Ci": 20,
        "circletube": 21,
        "Ct": 21,
        "isection": 30,
        "I": 30,
        "spring": 40,
        "Sp": 40,
    }
    idsecdef = matdef[keysecdef]
    return idsecdef


def sect_prop(sec_set: str, dim_sec: dict):
    """cross section property

    Arguments:
        sec_set:str -- section setting
        dim_sec:dict{} -- section's dimensions

    Returns:
        A:float   -- area
        Izz:float -- inercia zz
        Iyy:float -- inercia yy
        Jxx:float -- inercia xx
    """
    b = dim_sec["b"]
    h = dim_sec["h"]
    t = dim_sec["t"]
    d = dim_sec["d"]
    if (sec_set == "rectangle") or (sec_set == "R"):
        A = b * h
        Izz = (1 / 12) * b * h**3
        Iyy = (1 / 12) * h * b**3
        Jxx = Iyy + Izz
    elif (sec_set == "rectangletube") or (sec_set == "Rt"):
        A = b * h - ((b - 2 * t) * (h - 2 * t))
        Izz = (1 / 12) * (b * h**3) - (1 / 12) * ((b - 2 * t) * (h - 2 * t) ** 3)
        Iyy = (1 / 12) * (h * b**3) - (1 / 12) * ((h - 2 * t) * (b - 2 * t) ** 3)
        Jxx = Iyy + Izz
    elif (sec_set == "circle") or (sec_set == "Ci"):
        A = (1 / 4) * np.pi * d**2
        Izz = (1 / 64) * np.pi * d**4
        Iyy = Izz
        Jxx = Iyy + Izz
    elif (sec_set == "circletube") or (sec_set == "Ct"):
        A = (1 / 4) * np.pi * (d**2 - (d - 2 * t) ** 2)
        Izz = (1 / 64) * np.pi * (d**4 - (d - 2 * t) ** 4)
        Iyy = Izz
        Jxx = Iyy + Izz
    elif (sec_set == "isection") or (sec_set == "I"):
        A = 2 * b * d + t * (h - 2 * d)
        Izz = (b * h**3) / 12 - ((b - t) * (h - 2 * d) ** 3) / 12
        Iyy = ((h - 2 * d) * t**3) / 12 + 2 * (d * b**3) / 12
        Jxx = Iyy + Izz
    elif (sec_set == "spring") or (sec_set == "S"):
        A = 1
        Izz = 1
        Iyy = 1
        Jxx = 1
    return A, Izz, Iyy, Jxx


def cg_coord(tabgeo: np.ndarray, inci: np.ndarray, num_elm: int):
    """coord cg compute

    Arguments:
        tabgeo:list[]  -- table of geometry prop.
        inci:list[]    -- elements conection and prop. list
        num_elm:int    -- element(in mesh) number

    Returns:
        CG:np.array  -- coord of CG
    """
    b = tabgeo[int(inci[num_elm, 3]) - 1, 5]
    h = tabgeo[int(inci[num_elm, 3]) - 1, 6]
    t = tabgeo[int(inci[num_elm, 3]) - 1, 7]
    d = tabgeo[int(inci[num_elm, 3]) - 1, 8]
    y_max = 1
    y_min = 1
    z_max = 1
    z_min = 1
    r_max = 1
    if tabgeo[int(inci[num_elm, 3]) - 1, -1] == 10:
        y_max = h * 0.5
        y_min = -h * 0.5
        z_max = b * 0.5
        z_min = -b * 0.5
        r_max = y_max
    elif tabgeo[int(inci[num_elm, 3]) - 1, -1] == 11:
        y_max = h * 0.5
        y_min = -h * 0.5
        z_max = b * 0.5
        z_min = -b * 0.5
        r_max = y_max
    elif tabgeo[int(inci[num_elm, 3]) - 1, -1] == 20:
        y_max = d * 0.5
        y_min = -d * 0.5
        z_max = d * 0.5
        z_min = -d * 0.5
        r_max = y_max
    elif tabgeo[int(inci[num_elm, 3]) - 1, -1] == 21:
        y_max = d * 0.5
        y_min = -d * 0.5
        z_max = d * 0.5
        z_min = -d * 0.5
        r_max = y_max
    elif tabgeo[int(inci[num_elm, 3]) - 1, -1] == 30:
        y_max = h * 0.5
        y_min = -h * 0.5
        z_max = b * 0.5
        z_min = -b * 0.5
        r_max = y_max
    elif tabgeo[int(inci[num_elm, 3]) - 1, -1] == 40:
        y_max = d * 0.5
        y_min = -d * 0.5
        z_max = d * 0.5
        z_min = -d * 0.5
        r_max = y_max
    CG = np.array([y_max, y_min, z_max, z_min, r_max])
    return CG
