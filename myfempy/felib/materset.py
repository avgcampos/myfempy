#!/usr/bin/env python
"""
Material Setting
"""

def mat_def(keymatdef):
    matdef = {
        'springlinear': 10,
        'springnonlin': 11,
        'isotropic': 20,
        'orthotropic': 30,
    }
    idmatdef = matdef[keymatdef]
    return idmatdef


def mat_beh(keymatbeh):
    matbeh = {
        'lumped': 1,
        'axial': 2,
        'planestress': 3,
        'planestrain': 4,
        'solid': 5,
    }
    idmatbeh = matbeh[keymatbeh]
    return idmatbeh


def get_elasticity(tabmat, inci, num_elm):
    if tabmat[int(inci[num_elm, 2])-1, -1] == 1:
        if tabmat[int(inci[num_elm, 2])-1, -2] == 10:
            from myfempy.felib.materials.lumped import Elasticity
            mat = Elasticity(tabmat, inci, num_elm)
            return mat.springlinear()
        elif tabmat[int(inci[num_elm, 2])-1, -2] == 11:
            print('Not implemented yet')
    elif tabmat[int(inci[num_elm, 2])-1, -1] == 2:
        if tabmat[int(inci[num_elm, 2])-1, -2] == 20:
            from myfempy.felib.materials.axial import Elasticity
            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2])-1, -2] == 30:
            print('Not implemented yet')
    elif tabmat[int(inci[num_elm, 2])-1, -1] == 3:
        if tabmat[int(inci[num_elm, 2])-1, -2] == 20:
            from myfempy.felib.materials.planestress import Elasticity
            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2])-1, -2] == 30:
            print('Not implemented yet')
    elif tabmat[int(inci[num_elm, 2])-1, -1] == 4:
        if tabmat[int(inci[num_elm, 2])-1, -2] == 20:
            from myfempy.felib.materials.planestrain import Elasticity
            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2])-1, -2] == 30:
            print('Not implemented yet')
    elif tabmat[int(inci[num_elm, 2])-1, -1] == 5:
        if tabmat[int(inci[num_elm, 2])-1, -2] == 20:
            from myfempy.felib.materials.solid import Elasticity
            mat = Elasticity(tabmat, inci, num_elm)
            return mat.isotropic()
        elif tabmat[int(inci[num_elm, 2])-1, -2] == 30:
            print('Not implemented yet')
