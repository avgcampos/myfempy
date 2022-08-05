#!/usr/bin/env python
__doc__ = """
lumped.py: SpringLinear material
"""


class Elasticity:
    def __init__(self, tabmat, inci, num_elm):
        self.S = tabmat[int(inci[num_elm, 2]) - 1, 7]  # spring stiffness
        self.D = tabmat[int(inci[num_elm, 2]) - 1, 8]  # spring dampe

    def springlinear(self):
        S = self.S
        D = self.D
        E = [S, D]
        return E


if __name__ == "__main__":
    import doctest

    doctest.testmod()
