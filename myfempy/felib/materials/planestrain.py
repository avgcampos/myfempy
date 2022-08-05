#!/usr/bin/env python
import numpy as np

__doc__ = """
planestrain.py: Plane Strain Isotropic material
"""


class Elasticity:
    def __init__(self, tabmat, inci, num_elm):
        self.E = tabmat[int(inci[num_elm, 2]) - 1, 0]  # material elasticity
        self.v = tabmat[int(inci[num_elm, 2]) - 1, 1]  # material poisson ratio

    def isotropic(self):
        D = np.zeros((3, 3))
        D[0, 0] = self.E * (1.0 - self.v) / ((1 + self.v) * (1.0 - 2.0 * self.v))
        D[0, 1] = D[0, 0] * self.v / (1.0 - self.v)
        D[1, 0] = D[0, 1]
        D[1, 1] = D[0, 0]
        D[2, 2] = D[0, 0] * 0.5 * (1.0 - 2.0 * self.v) / (1.0 - self.v)
        return D


if __name__ == "__main__":
    import doctest

    doctest.testmod()
