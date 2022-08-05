#!/usr/bin/env python
__doc__ = """
Finite Elements Setting
"""


def get_elemset(keyelem):
    if (keyelem == "spring21") or (keyelem == 110):
        from myfempy.felib.struct.spring21 import Spring21

        return Spring21
    elif (keyelem == "truss21") or (keyelem == 120):
        from myfempy.felib.struct.truss21 import Truss21

        return Truss21
    elif (keyelem == "beam21") or (keyelem == 130):
        from myfempy.felib.struct.beam21 import Beam21

        return Beam21
    elif (keyelem == "frame21") or (keyelem == 140):
        from myfempy.felib.struct.frame21 import Frame21

        return Frame21
    elif (keyelem == "frame22") or (keyelem == 141):
        from myfempy.felib.struct.frame22 import Frame22

        return Frame22
    elif (keyelem == "plane31") or (keyelem == 210):
        from myfempy.felib.struct.plane31 import Plane31

        return Plane31
    elif (keyelem == "plane41") or (keyelem == 220):
        from myfempy.felib.struct.plane41 import Plane41

        return Plane41
    elif (keyelem == "solid41") or (keyelem == 310):
        from myfempy.felib.struct.solid41 import Solid41

        return Solid41
    elif (keyelem == "solid81") or (keyelem == 320):
        from myfempy.felib.struct.solid81 import Solid81

        return Solid81
    else:
        print("erro idelem not defined")
