# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""

import numpy as np
from myfempy.felib.materset import get_elasticity
from myfempy.felib.quadrature import Quadrature

# %%------------------------------------------------------------------------------


class Solid41:
    '''Tetrahedron Isoparametric Solid 8-node linear Finite Element'''

    def __init__(self, modelinfo):

        self.dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        self.nodecon = modelinfo['nodecon'][0]
        self.fulldof = modelinfo["nodedof"][0]*len(modelinfo["coord"])
        self.nodedof = modelinfo['nodedof'][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode = len(modelinfo["coord"])
        self.inci = modelinfo['inci']
        self.coord = modelinfo['coord']
        self.tabmat = modelinfo['tabmat']
        self.tabgeo = modelinfo['tabgeo']
        self.ntensor = modelinfo['ntensor'][0]

        if modelinfo['quadra'][0] == 1:
            self.npp = modelinfo['quadra'][1]
            self.quadra = Quadrature.gaussian(self.npp)

    def elemset(self):

        dofelem = {'key': 'solid41',
                   'id': 310,
                   'def': 'struct 3D',
                   'dofs': ['ux', 'uy', 'uz'],
                   'nnodes': ['i', 'j', 'k', 'l'],
                   'tensor': ['sxx', 'syy', 'szz', 'sxy', 'syz', 'szx']}

        return dofelem

    def lockey(self, nodelist):

        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]

        loc = np.array([self.nodedof*noi-3, self.nodedof*noi-2, self.nodedof*noi-1,
                        self.nodedof*noj-3, self.nodedof*noj-2, self.nodedof*noj-1,
                        self.nodedof*nok-3, self.nodedof*nok-2, self.nodedof*nok-1,
                        self.nodedof*nol-3, self.nodedof*nol-2, self.nodedof*nol-1])
        return loc

    def matriz_B(self, nodelist, intpl):

        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]

        xi = self.coord[noi-1, 1]
        yi = self.coord[noi-1, 2]
        zi = self.coord[noi-1, 3]
        xj = self.coord[noj-1, 1]
        yj = self.coord[noj-1, 2]
        zj = self.coord[noj-1, 3]
        xk = self.coord[nok-1, 1]
        yk = self.coord[nok-1, 2]
        zk = self.coord[nok-1, 3]
        xl = self.coord[nol-1, 1]
        yl = self.coord[nol-1, 2]
        zl = self.coord[nol-1, 3]

        C = np.array([[1, xi, yi, zi], [1, xj, yj, zj],
                     [1, xk, yk, zk], [1, xl, yl, zl]])

        V = (1/6)*abs(np.linalg.det(C))

        # old LOGAN
        alpi = np.linalg.det(
            np.array([[xj, yj, zj], [xk, yk, zk], [xl, yl, zl]]))
        beti = - \
            np.linalg.det(np.array([[1, yj, zj], [1, yk, zk], [1, yl, zl]]))
        gami = np.linalg.det(np.array([[1, xj, zj], [1, xk, zk], [1, xl, zl]]))
        deli = - \
            np.linalg.det(np.array([[1, xj, yj], [1, xk, yk], [1, xl, yl]]))

        alpj = - \
            np.linalg.det(np.array([[xi, yi, zi], [xk, yk, zk], [xl, yl, zl]]))
        betj = np.linalg.det(np.array([[1, yi, zi], [1, yk, zk], [1, yl, zl]]))
        gamj = - \
            np.linalg.det(np.array([[1, xi, zi], [1, xk, zk], [1, xl, zl]]))
        delj = np.linalg.det(np.array([[1, xi, yi], [1, xk, yk], [1, xl, yl]]))

        alpk = np.linalg.det(
            np.array([[xi, yi, zi], [xj, yj, zj], [xl, yl, zl]]))
        betk = - \
            np.linalg.det(np.array([[1, yi, zi], [1, yj, zj], [1, yl, zl]]))
        gamk = np.linalg.det(np.array([[1, xi, zi], [1, xj, zj], [1, xl, zl]]))
        delk = - \
            np.linalg.det(np.array([[1, xi, yi], [1, xj, yj], [1, xl, yl]]))

        alpl = - \
            np.linalg.det(np.array([[xi, yi, zi], [xj, yj, zj], [xk, yk, zk]]))
        betl = np.linalg.det(np.array([[1, yi, zi], [1, yj, zj], [1, yk, zk]]))
        gaml = - \
            np.linalg.det(np.array([[1, xi, zi], [1, xj, zj], [1, xk, zk]]))
        dell = np.linalg.det(np.array([[1, xi, yi], [1, xj, yj], [1, xk, yk]]))

        # new LIU
        # alpi = np.linalg.det(
        #     np.array([[xj, yj, zj], [xk, yk, zk], [xl, yl, zl]]))
        # beti = - \
        #     np.linalg.det(np.array([[1, yj, zj], [1, yk, zk], [1, yl, zl]]))
        # gami = - \
        #     np.linalg.det(np.array([[yj, 1, zj], [yk, 1, zk], [yl, 1, zl]]))
        # deli = - \
        #     np.linalg.det(np.array([[yj, zj, 1], [yk, zk, 1], [yl, zl, 1]]))

        # alpj = np.linalg.det(
        #     np.array([[xk, yk, zk], [xl, yl, zl], [xi, yi, zi]]))
        # betj = - \
        #     np.linalg.det(np.array([[1, yk, zk], [1, yl, zl], [1, yi, zi]]))
        # gamj = - \
        #     np.linalg.det(np.array([[yk, 1, zk], [yl, 1, zl], [yi, 1, zi]]))
        # delj = - \
        #     np.linalg.det(np.array([[yk, zk, 1], [yl, zl, 1], [yi, zi, 1]]))

        # alpk = np.linalg.det(
        #     np.array([[xl, yl, zl], [xi, yi, zi], [xj, yj, zj]]))
        # betk = - \
        #     np.linalg.det(np.array([[1, yl, zl], [1, yi, zi], [1, yj, zj]]))
        # gamk = - \
        #     np.linalg.det(np.array([[yl, 1, zl], [yi, 1, zi], [yj, 1, zj]]))
        # delk = - \
        #     np.linalg.det(np.array([[yl, zl, 1], [yi, zi, 1], [yj, zj, 1]]))

        # alpl = np.linalg.det(
        #     np.array([[xi, yi, zi], [xj, yj, zj], [xk, yk, zk]]))
        # betl = - \
        #     np.linalg.det(np.array([[1, yi, zi], [1, yj, zj], [1, yk, zk]]))
        # gaml = - \
        #     np.linalg.det(np.array([[yi, 1, zi], [yj, 1, zj], [yk, 1, zk]]))
        # dell = - \
        #     np.linalg.det(np.array([[yi, zi, 1], [yj, zj, 1], [yk, zk, 1]]))

        B = np.zeros((self.ntensor, self.dofe))

        B[0, 0] = beti
        B[0, 3] = betj
        B[0, 6] = betk
        B[0, 9] = betl

        B[1, 1] = gami
        B[1, 4] = gamj
        B[1, 7] = gamk
        B[1, 10] = gaml

        B[2, 2] = deli
        B[2, 5] = delj
        B[2, 8] = delk
        B[2, 11] = dell

        B[3, 0] = gami
        B[3, 1] = beti
        B[3, 3] = gamj
        B[3, 4] = betj
        B[3, 6] = gamk
        B[3, 7] = betk
        B[3, 9] = gaml
        B[3, 10] = betl

        B[4, 1] = deli
        B[4, 2] = gami
        B[4, 4] = delj
        B[4, 5] = gamj
        B[4, 7] = delk
        B[4, 8] = gamk
        B[4, 10] = dell
        B[4, 11] = gaml

        B[5, 0] = deli
        B[5, 2] = beti
        B[5, 3] = delj
        B[5, 5] = betj
        B[5, 6] = delk
        B[5, 8] = betk
        B[5, 9] = dell
        B[5, 11] = betl

        B = (1/(6*V))*B

        return B, V

    # montagem matriz triangular CST 2D Sparse
    def stiff_linear(self, ee):

        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nol = int(self.inci[ee, 7])

        nodelist = [noi, noj, nok, nol]

        D = get_elasticity(self.tabmat, self.inci, ee)

        intpl = 0.0

        B, V = Solid41.matriz_B(self, nodelist, intpl)
        ket4 = V*np.dot(np.dot(np.transpose(B), D), B)

        loc = Solid41.lockey(self, nodelist)

        return ket4, loc

    def mass(self, ee):

        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nol = int(self.inci[ee, 7])

        xi = self.coord[noi-1, 1]
        yi = self.coord[noi-1, 2]
        zi = self.coord[noi-1, 3]
        xj = self.coord[noj-1, 1]
        yj = self.coord[noj-1, 2]
        zj = self.coord[noj-1, 3]
        xk = self.coord[nok-1, 1]
        yk = self.coord[nok-1, 2]
        zk = self.coord[nok-1, 3]
        xl = self.coord[nol-1, 1]
        yl = self.coord[nol-1, 2]
        zl = self.coord[nol-1, 3]

        C = np.array([[1, xi, yi, zi],
                      [1, xj, yj, zj],
                      [1, xk, yk, zk],
                      [1, xl, yl, zl]])

        V = (1/6)*abs(np.linalg.det(C))

        # alpi = np.linalg.det(np.array([[xj,yj,zj],[xk,yk,zk],[xl,yl,zl]]))
        # beti = -np.linalg.det(np.array([[1,yj,zj],[1,yk,zk],[1,yl,zl]]))
        # gami = -np.linalg.det(np.array([[yj,1,zj],[yk,1,zk],[yl,1,zl]]))
        # deli = -np.linalg.det(np.array([[yj,zj,1],[yk,zk,1],[yl,zl,1]]))

        # alpj = np.linalg.det(np.array([[xk,yk,zk],[xl,yl,zl],[xi,yi,zi]]))
        # betj = -np.linalg.det(np.array([[1,yk,zk],[1,yl,zl],[1,yi,zi]]))
        # gamj = -np.linalg.det(np.array([[yk,1,zk],[yl,1,zl],[yi,1,zi]]))
        # delj = -np.linalg.det(np.array([[yk,zk,1],[yl,zl,1],[yi,zi,1]]))

        # alpk = np.linalg.det(np.array([[xl,yl,zl],[xi,yi,zi],[xj,yj,zj]]))
        # betk = -np.linalg.det(np.array([[1,yl,zl],[1,yi,zi],[1,yj,zj]]))
        # gamk = -np.linalg.det(np.array([[yl,1,zl],[yi,1,zi],[yj,1,zj]]))
        # delk = -np.linalg.det(np.array([[yl,zl,1],[yi,zi,1],[yj,zj,1]]))

        # alpl = np.linalg.det(np.array([[xi,yi,zi],[xj,yj,zj],[xk,yk,zk]]))
        # betl = -np.linalg.det(np.array([[1,yi,zi],[1,yj,zj],[1,yk,zk]]))
        # gaml = -np.linalg.det(np.array([[yi,1,zi],[yj,1,zj],[yk,1,zk]]))
        # dell = -np.linalg.det(np.array([[yi,zi,1],[yj,zj,1],[yk,zk,1]]))

        #

        # N1 = (alpi+beti*xi+gami*yi+deli*zi)
        # N2 = (alpj+betj*xj+gamj*yj+delj*zj)
        # N3 = (alpk+betk*xk+gamk*yk+delk*zk)
        # N4 = (alpl+betl*xl+gaml*yl+dell*zl)

        # N = (1/(6*V))*np.array([[N1,0,0,N2,0,0,N3,0,0,N4,0,0],\
        #                         [0,N1,0,0,N2,0,0,N3,0,0,N4,0],\
        #                         [0,0,N1,0,0,N2,0,0,N3,0,0,N4]])

        # met8 = np.dot(np.transpose(N),N)*R*V

        R = self.tabmat[int(self.inci[ee, 2]-1), 6]

        met4 = np.zeros((self.dofe, self.dofe))
        met4[0, 0] = 2.0
        met4[0, 3] = 1.0
        met4[0, 6] = 1.0
        met4[0, 9] = 1.0

        met4[1, 1] = 2.0
        met4[1, 4] = 1.0
        met4[1, 7] = 1.0
        met4[1, 10] = 1.0

        met4[2, 2] = 2.0
        met4[2, 5] = 1.0
        met4[2, 8] = 1.0
        met4[2, 11] = 1.0

        met4[3, 0] = 1.0
        met4[3, 3] = 2.0
        met4[3, 6] = 1.0
        met4[3, 9] = 1.0

        met4[4, 1] = 1.0
        met4[4, 4] = 2.0
        met4[4, 7] = 1.0
        met4[4, 10] = 1.0

        met4[5, 2] = 1.0
        met4[5, 5] = 2.0
        met4[5, 8] = 1.0
        met4[5, 11] = 1.0

        met4[6, 0] = 1.0
        met4[6, 3] = 1.0
        met4[6, 6] = 2.0
        met4[6, 9] = 1.0

        met4[7, 1] = 1.0
        met4[7, 4] = 1.0
        met4[7, 7] = 2.0
        met4[7, 10] = 1.0

        met4[8, 2] = 1.0
        met4[8, 5] = 1.0
        met4[8, 8] = 2.0
        met4[8, 11] = 1.0

        met4[9, 0] = 1.0
        met4[9, 3] = 1.0
        met4[9, 6] = 1.0
        met4[9, 9] = 2.0

        met4[10, 1] = 1.0
        met4[10, 4] = 1.0
        met4[10, 7] = 1.0
        met4[10, 10] = 2.0

        met4[11, 2] = 1.0
        met4[11, 5] = 1.0
        met4[11, 8] = 1.0
        met4[11, 11] = 2.0

        met4 = ((R*V)/20)*met4

        list_node = [noi, noj, nok, nol]
        loc = Solid41.lockey(self, list_node)

        return met4, loc
