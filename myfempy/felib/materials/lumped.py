#!/usr/bin/env python
__doc__ = """
lumped.py: SpringLinear material
"""


class Elasticity:
    """elasticity set class""" 

    def __init__(self, tabmat: np.ndarray, inci: np.ndarray, num_elm: int):
        self.K = tabmat[int(inci[num_elm, 2]) - 1, 7]  # spring stiffness
        self.C = tabmat[int(inci[num_elm, 2]) - 1, 8]  # spring dampe

    def springlinear(self):
        """spring material lumped

        Returns:
             D:list[]  -- elasticity matrix
        """
        K = self.K
        C = self.C
        D = [K, C]
        return D


if __name__ == "__main__":
    import doctest

    doctest.testmod()
