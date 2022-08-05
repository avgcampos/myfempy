#!/usr/bin/env python
from myfempy.felib.materset import get_elasticity
import numpy as np
__doc__ = """
frame21.py: Frame 2D 2-node linear Finite Element
"""


class Frame21:
    def __init__(self, modelinfo):

        self.dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        self.fulldof = modelinfo["nodedof"][0]*len(modelinfo["coord"])
        self.nodedof = modelinfo['nodedof'][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode = len(modelinfo["coord"])
        self.inci = modelinfo['inci']
        self.coord = modelinfo['coord']
        self.tabmat = modelinfo['tabmat']
        self.tabgeo = modelinfo['tabgeo']

    @staticmethod
    def elemset():
        dofelem = {'key': 'frame21',
                   'id': 140,
                   'def': 'struct 1D',
                   'dofs': ['ux', 'uy', 'rz'],
                   'nnodes': ['i', 'j'],
                   'tensor': ['sxx']}

        return dofelem

    def lockey(self, list_node):
        noi = list_node[0]
        noj = list_node[1]
        loc = np.array([self.nodedof*noi-3, self.nodedof*noi-2, self.nodedof*noi-1,
                        self.nodedof*noj-3, self.nodedof*noj-2, self.nodedof*noj-1])
        return loc

    def stiff_linear(self, ee):
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi-1, 1]
        noiy = self.coord[noi-1, 2]
        nojx = self.coord[noj-1, 1]
        nojy = self.coord[noj-1, 2]
        D = get_elasticity(self.tabmat, self.inci, ee)
        E = D[0]
        A = self.tabgeo[int(self.inci[ee, 3]-1), 0]
        I33 = self.tabgeo[int(self.inci[ee, 3]-1), 1]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        T = np.zeros((self.dofe, self.dofe))
        T[0, 0] = c
        T[0, 1] = s
        T[1, 0] = -s
        T[1, 1] = c
        T[2, 2] = 1.0
        T[3, 3] = c
        T[3, 4] = s
        T[4, 3] = -s
        T[4, 4] = c
        T[5, 5] = 1.0
        kef2 = np.zeros((self.dofe, self.dofe))
        kef2[0, 0] = (A*E)/L
        kef2[0, 3] = -(A*E)/L
        kef2[1, 1] = 12*(E*I33)/L**3
        kef2[1, 2] = 6*(E*I33)/L**2
        kef2[1, 4] = -12*(E*I33)/L**3
        kef2[1, 5] = 6*(E*I33)/L**2
        kef2[2, 1] = 6*(E*I33)/L**2
        kef2[2, 2] = 4*(E*I33)/L
        kef2[2, 4] = -6*(E*I33)/L**2
        kef2[2, 5] = 2*(E*I33)/L
        kef2[3, 0] = -(A*E)/L
        kef2[3, 3] = (A*E)/L
        kef2[4, 1] = -12*(E*I33)/L**3
        kef2[4, 2] = -6*(E*I33)/L**2
        kef2[4, 4] = 12*(E*I33)/L**3
        kef2[4, 5] = -6*(E*I33)/L**2
        kef2[5, 1] = 6*(E*I33)/L**2
        kef2[5, 2] = 2*(E*I33)/L
        kef2[5, 4] = -6*(E*I33)/L**2
        kef2[5, 5] = 4*(E*I33)/L
        kef2 = (E*I33/L**3)*kef2
        kef2t = np.dot(np.dot(np.transpose(T), kef2), T)
        list_node = [noi, noj]
        loc = Frame21.lockey(self, list_node)
        return kef2t, loc

    def mass(self, ee):
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi-1, 1]
        noiy = self.coord[noi-1, 2]
        nojx = self.coord[noj-1, 1]
        nojy = self.coord[noj-1, 2]
        R = self.tabmat[int(self.inci[ee, 2]-1), 6]
        A = self.tabgeo[int(self.inci[ee, 3]-1), 0]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        mef2d2 = np.zeros((self.dofe, self.dofe))
        mef2d2[0, 0] = 140
        mef2d2[0, 3] = 70
        mef2d2[1, 1] = 156
        mef2d2[1, 2] = 22*L
        mef2d2[1, 4] = 54
        mef2d2[1, 5] = -13*L
        mef2d2[2, 1] = 22*L
        mef2d2[2, 2] = 4*L**2
        mef2d2[2, 4] = 13*L
        mef2d2[2, 5] = -3*L**2
        mef2d2[3, 0] = 70
        mef2d2[3, 3] = 140
        mef2d2[4, 1] = 54
        mef2d2[4, 2] = 13*L
        mef2d2[4, 4] = 156
        mef2d2[4, 5] = -22*L
        mef2d2[5, 1] = -13*L
        mef2d2[5, 2] = -3*L**2
        mef2d2[5, 4] = -22*L
        mef2d2[5, 5] = 4*L**2
        mef2d2 = (R*A*L/420)*mef2d2
        list_node = [noi, noj]
        loc = Frame21.lockey(self, list_node)
        return mef2d2, loc

    def intforces(self, U, lines):
        Fint = np.zeros((self.fulldof), dtype=float)
        Nx = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        Vy = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        Mz = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        domL = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        for ee in range(self.nelem):
            noi = int(self.inci[ee, 4])
            noj = int(self.inci[ee, 5])
            noix = self.coord[noi-1, 1]
            noiy = self.coord[noi-1, 2]
            nojx = self.coord[noj-1, 1]
            nojy = self.coord[noj-1, 2]
            L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
            s = (nojy-noiy)/L
            c = (nojx-noix)/L
            T = np.zeros((self.dofe, self.dofe))
            T[0, 0] = c
            T[0, 1] = s
            T[1, 0] = -s
            T[1, 1] = c
            T[2, 2] = 1.0
            T[3, 3] = c
            T[3, 4] = s
            T[4, 3] = -s
            T[4, 4] = c
            T[5, 5] = 1.0
            kef2T, loc = Frame21.stiff_linear(self, ee)
            kef2 = np.dot(np.dot(np.transpose(T), kef2T), T)
            F = np.dot(kef2, T@U[loc])
            Fint[loc] = [-F[0], F[1], -F[2], F[3], -F[4], F[5]]
        for ed in range(len(lines)):
            nnodes = lines[ed][1]
            for nn in range(0, len(nnodes)):
                node = int(nnodes[nn])
                Nx[nn, ed] = Fint[self.nodedof*node-3]
                Vy[nn, ed] = Fint[self.nodedof*node-2]
                Mz[nn, ed] = Fint[self.nodedof*node-1]
                domL[nn, ed] = self.coord[node-1, 1]
            domLIdc = np.argsort(domL[:, ed], axis=0)
            domL[:, ed] = domL[domLIdc, ed]
            Mz[:, ed] = Mz[domLIdc, ed]
            Vy[:, ed] = Vy[domLIdc, ed]
            Nx[:, ed] = Nx[domLIdc, ed]
        ifb = {'le': domL, 'val': [Nx, Vy, Mz]}
        title = ["NX", "VY", "MZ"]
        return ifb, title

    def matrix_b(self, ee, csc):
        y = csc[0]
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi-1, 1]
        noiy = self.coord[noi-1, 2]
        nojx = self.coord[noj-1, 1]
        nojy = self.coord[noj-1, 2]
        D = get_elasticity(self.tabmat, self.inci, ee)
        E = D[0]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        T = np.array([[c, s, 0, 0, 0, 0],
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])
        coord_local = np.dot(T, np.array(
            [[noix], [noiy], [1], [nojx], [nojy], [1]]))
        x_mid = (coord_local[3][0] - coord_local[0][0])/2
        B = np.array([-E/L, -(y*E)*((12*x_mid)/(L**3) - 6/(L**2)),
                      -(y*E)*((6*x_mid)/(L**2) - 4/L),
                      E/L, -(y*E)*(-(12*x_mid)/(L**3) + 6/(L**2)),
                      -(y*E)*((6*x_mid)/(L**2) - 2/L)])
        return B, T


if __name__ == "__main__":
    import doctest
    doctest.testmod()
