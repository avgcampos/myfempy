#!/usr/bin/env python
"""
lumped.py: SpringLinear material
"""
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"


class Elasticity:
    def __init__(self, tabmat, inci, num_elm):
        self.S = tabmat[int(inci[num_elm, 2])-1, 7]  # spring stiffness
        self.D = tabmat[int(inci[num_elm, 2])-1, 8]  # spring dampe

    def springlinear(self):
        S = self.S
        D = self.D
        E = [S, D]
        return E
