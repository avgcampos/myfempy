#!/usr/bin/env python
"""
Quadrature
"""
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"
import numpy as np
sqrt3 = np.sqrt(3.0)


class Quadrature:
    def no_interpol(npp):
        if npp == 1:
            xp = np.array([[1.0], [1.0], [1.0]])
            wp = np.array([1])
            return xp, wp

    def gaussian(npp):
        if npp == 1:
            xp = np.array([[1/3], [1/3]])
            wp = np.array([1])
            return xp, wp
        elif npp == 2:
            xp = np.array([[-0.577350269], [0.577350269]])
            wp = np.array([1])
            return xp, wp
        elif npp == 3:
            xp = np.array([[2/3, 1/6, 1/6],
                           [1/6, 1/6, 2/3]])
            wp = np.array([1/3, 1/3, 1/3])
            return xp, wp
        elif npp == 4:
            xp = np.array([[-1/sqrt3, 1/sqrt3, 1/sqrt3, -1/sqrt3],
                           [-1/sqrt3, -1/sqrt3, 1/sqrt3, 1/sqrt3]])
            wp = np.array([1, 1, 1, 1])
            return xp, wp
        elif npp == 8:
            xp = np.array([[-1/sqrt3, 1/sqrt3, 1/sqrt3, -1/sqrt3, -1/sqrt3, 1/sqrt3, 1/sqrt3, -1/sqrt3],
                           [-1/sqrt3, -1/sqrt3, 1/sqrt3, 1/sqrt3, -
                               1/sqrt3, -1/sqrt3, 1/sqrt3, 1/sqrt3],
                           [-1/sqrt3, -1/sqrt3, -1/sqrt3, -1/sqrt3, 1/sqrt3, 1/sqrt3, 1/sqrt3, 1/sqrt3]])
            wp = np.array([1, 1, 1, 1, 1, 1, 1, 1])
            return xp, wp
        else:
            print('erro pontos de integração Gauss')
