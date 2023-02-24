#!/usr/bin/env python
__doc__ = """
Material Setting
"""
import numpy as np

def mat_def(keymatdef: str):
    """material def

    Arguments:
       keymatdef:str -- key material def

    Returns:
       idmatdef:int  -- id number of cross section 
    """
    matdef = {
        "springlinear": 10,
        "springnonlin": 11,
        "isotropic": 20,
        "orthotropic": 30,
    }
    idmatdef = matdef[keymatdef]
    return idmatdef


def mat_beh(keymatbeh: str):
    """material behavior

    Arguments:
        keymatbeh:str -- key material behavior

    Returns:
        idmatbeh:int  -- id mat. beh.
    """
    matbeh = {
        "lumped": 1,
        "axial": 2,
        "planestress": 3,
        "planestrain": 4,
        "solid": 5,
    }
    idmatbeh = matbeh[keymatbeh]
    return idmatbeh


def get_elasticity(tabmat: np.ndarray, inci: np.ndarray, num_elm: int):
    """get elasticity matrix D

    Arguments:
        tabmat:list[] -- table of material prop.
        inci:list[]   -- elements conection and prop. list
        num_elm:int   -- element(in mesh) number

    Returns:
       elasticity class
    """
    if tabmat[int(inci[num_elm, 2]) - 1, -1] == 1:
        if tabmat[int(inci[num_elm, 2]) - 1, -2] == 10:
            from myfempy.felib.materials.lumped import Elasticity

            mat = Elasticity(tabmat, inci, num_elm)
            return mat.springlinear()
        elif tabmat[int(inci[num_elm, 2]) - 1, -2] == 11:
            print("Not implemented yet")
    elif tabmat[int(inci[num_elm, 2]) - 1, -1] == 2:
        if tabmat[int(inci[num_elm, 2]) - 1, -2] == 20:
            from myfempy.felib.materials.axial import Elasticity

            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2]) - 1, -2] == 30:
            print("Not implemented yet")
    elif tabmat[int(inci[num_elm, 2]) - 1, -1] == 3:
        if tabmat[int(inci[num_elm, 2]) - 1, -2] == 20:
            from myfempy.felib.materials.planestress import Elasticity

            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2]) - 1, -2] == 30:
            print("Not implemented yet")
    elif tabmat[int(inci[num_elm, 2]) - 1, -1] == 4:
        if tabmat[int(inci[num_elm, 2]) - 1, -2] == 20:
            from myfempy.felib.materials.planestrain import Elasticity

            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2]) - 1, -2] == 30:
            print("Not implemented yet")
    elif tabmat[int(inci[num_elm, 2]) - 1, -1] == 5:
        if tabmat[int(inci[num_elm, 2]) - 1, -2] == 20:
            from myfempy.felib.materials.solid import Elasticity

            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2]) - 1, -2] == 30:
            print("Not implemented yet")
    else:
        print('erro')