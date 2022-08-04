#!/usr/bin/env python
__doc__ ="""
planestress.py: Plane Stress Isotropic and Elasticity Material
"""
import numpy as np
from myfempy.felib.felemset import get_elemset
from myfempy.felib.quadrature import Quadrature


class Elasticity:
    def __init__(self, tabmat, inci, ee):
        self.E = tabmat[int(inci[ee, 2])-1, 0]  # material elasticity
        self.v = tabmat[int(inci[ee, 2])-1, 1]  # material poisson ratio

    def isotropic(self):
        D = np.zeros((3, 3))
        D[0, 0] = self.E/(1.0-self.v*self.v)
        D[0, 1] = D[0, 0]*self.v
        D[1, 0] = D[0, 1]
        D[1, 1] = D[0, 0]
        D[2, 2] = self.E/(2.0*(1.0+self.v))
        return D


class Tensor:

    def __init__(self, modelinfo, U, ee):
        self.ee = ee
        self.U = U
        self.modelinfo = modelinfo
        self.dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        self.ntensor = modelinfo['ntensor'][0]
        self.inci = modelinfo['inci']
        self.tabmat = modelinfo['tabmat']
        element = get_elemset(int(self.inci[ee][1]))
        self.setelement = element(self.modelinfo)
        listnodes = self.inci[ee, 4:]
        self.nodelist = listnodes[listnodes.nonzero()].astype(int).tolist()
        self.npp = 1
        self.xpp = [0.0]  # xp[0,:]
        self.ypp = [0.0]  # xp[1,:]
        self.zpp = [0.0]  # xp[2,:]
        self.wp = [2.0]

    def strain(self):
        B = np.zeros((self.ntensor, self.dofe))
        for pp in range(0, self.npp):
            intpl = [self.xpp[pp], self.ypp[pp]]
            Bpp, detJ = self.setelement.matriz_b(self.nodelist, intpl)
            B += Bpp
        loc = self.setelement.lockey(self.nodelist)
        epsilon = B@(self.U[loc])
        strn_elm_xx = epsilon[0]
        strn_elm_yy = epsilon[1]
        strn_elm_xy = epsilon[2]
        strn_elm_vm = np.sqrt(
            epsilon[0]**2 - epsilon[0]*epsilon[1] + epsilon[1]**2 + 3*epsilon[2]**2)
        strain = [strn_elm_vm, strn_elm_xx, strn_elm_yy, strn_elm_xy]
        title = ["STRAIN_VM", "STRAIN_XX", "STRAIN_YY", "STRAIN_XY"]
        return epsilon, strain, title

    def stress(self, epsilon):
        M = Elasticity(self.tabmat, self.inci, self.ee)
        D = M.isotropic()
        sigma = np.dot(D, epsilon)
        strs_elm_xx = sigma[0]
        strs_elm_yy = sigma[1]
        strs_elm_xy = sigma[2]
        strs_elm_vm = np.sqrt(sigma[0]**2 - sigma[0]
                              * sigma[1] + sigma[1]**2 + 3*sigma[2]**2)
        stress = [strs_elm_vm, strs_elm_xx, strs_elm_yy, strs_elm_xy]
        title = ["STRESS_VM", "STRESS_XX", "STRESS_YY", "STRESS_XY"]
        return stress, title
