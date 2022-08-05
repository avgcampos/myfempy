#!/usr/bin/env python
__doc__ = """
Displacement Calculator
"""
import numpy as np
from myfempy.tools.tools import loading_bar_v1


class Deformation:
    def __init__(self, modelinfo):
        self.nodedof = modelinfo["nodedof"][0]
        self.nnode = len(modelinfo["coord"])
        self.coord = modelinfo["coord"]

    def ux(self, U):
        Udef = np.zeros((self.nnode, 3), dtype=float)
        Umag = np.zeros((self.nnode, 1), dtype=float)
        loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, self.nnode + 1):
            loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 0] = U[self.nodedof * nn - 2]
            Udef[nn - 1, 1] = U[self.nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(
                U[self.nodedof * nn - 2] ** 2 + U[self.nodedof * nn - 1] ** 2
            )
        result = np.concatenate((Umag, Udef), axis=1)
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return result

    def uy_rz(self, U):
        Udef = np.zeros((self.nnode, 3), dtype=float)
        Umag = np.zeros((self.nnode, 1), dtype=float)
        loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, self.nnode + 1):
            loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 1] = U[self.nodedof * nn - 2]
            # Rdef[nn-1,2] = U[self.nodedof*nn-1,0]
            Umag[nn - 1, 0] = Udef[nn - 1, 1]
        result = np.concatenate((Umag, Udef), axis=1)
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]

        return result

    def ux_uy_rz(self, U):
        Udef = np.zeros((self.nnode, 3), dtype=float)
        Umag = np.zeros((self.nnode, 1), dtype=float)
        loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, self.nnode + 1):
            loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 0] = U[self.nodedof * nn - 3]
            Udef[nn - 1, 1] = U[self.nodedof * nn - 2]
            # Rdef[nn-1,2] = U[self.nodedof*nn-1,0]
            Umag[nn - 1, 0] = np.sqrt(
                U[self.nodedof * nn - 3] ** 2 + U[self.nodedof * nn - 2] ** 2
            )
        result = np.concatenate((Umag, Udef), axis=1)
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return result

    def ux_uy_uz_rx_ry_rz(self, U):
        Udef = np.zeros((self.nnode, 3), dtype=float)
        Umag = np.zeros((self.nnode, 1), dtype=float)
        for nn in range(1, self.nnode + 1):
            Udef[nn - 1, 0] = U[self.nodedof * nn - 6]
            Udef[nn - 1, 1] = U[self.nodedof * nn - 5]
            Udef[nn - 1, 2] = U[self.nodedof * nn - 4]
            Umag[nn - 1, 0] = np.sqrt(
                U[self.nodedof * nn - 6] ** 2
                + U[self.nodedof * nn - 5] ** 2
                + U[self.nodedof * nn - 4] ** 2
            )
        result = np.concatenate((Umag, Udef), axis=1)
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return result

    def ux_uy(self, U):
        Udef = np.zeros((self.nnode, 3), dtype=float)
        Umag = np.zeros((self.nnode, 1), dtype=float)
        loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, self.nnode + 1):
            loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 0] = U[self.nodedof * nn - 2]
            Udef[nn - 1, 1] = U[self.nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(
                U[self.nodedof * nn - 2] ** 2 + U[self.nodedof * nn - 1] ** 2
            )
        result = np.concatenate((Umag, Udef), axis=1)
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return result

    def ux_uy_uz(self, U):
        Udef = np.zeros((self.nnode, 3), dtype=float)
        Umag = np.zeros((self.nnode, 1), dtype=float)
        loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, self.nnode + 1):
            loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 0] = U[self.nodedof * nn - 3]
            Udef[nn - 1, 1] = U[self.nodedof * nn - 2]
            Udef[nn - 1, 2] = U[self.nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(
                U[self.nodedof * nn - 3] ** 2
                + U[self.nodedof * nn - 2] ** 2
                + U[self.nodedof * nn - 1] ** 2
            )
        result = np.concatenate((Umag, Udef), axis=1)
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return result
