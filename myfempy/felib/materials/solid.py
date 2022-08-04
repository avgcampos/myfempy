#!/usr/bin/env python
__doc__ ="""
solid.py: Solid Isotropic and Elasticity Material
"""
import numpy as np
from myfempy.felib.felemset import get_elemset


class Elasticity:
    def __init__(self, tabmat, inci, num_elm):
        self.E = tabmat[int(inci[num_elm, 2])-1, 0]  # material elasticity
        self.v = tabmat[int(inci[num_elm, 2])-1, 1]  # material poisson ratio

    def isotropic(self):
        D = np.zeros((6, 6))
        fac = 1.0 / (2.0 * self.v * self.v + self.v - 1.0)
        D[0, 0] = fac * self.E * (self.v - 1.0)
        D[0, 1] = -1.0 * fac * self.E * self.v
        D[0, 2] = D[0, 1]
        D[1, 0] = D[0, 1]
        D[1, 1] = D[0, 0]
        D[1, 2] = D[0, 1]
        D[2, 0] = D[0, 1]
        D[2, 1] = D[0, 1]
        D[2, 2] = D[0, 0]
        D[3, 3] = self.E / (2.0 + 2.0 * self.v)
        D[4, 4] = D[3, 3]
        D[5, 5] = D[3, 3]
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
            intpl = [self.xpp[pp], self.ypp[pp], self.zpp[pp]]
            Bpp, detJ = self.setelement.matriz_b(self.nodelist, intpl)
            B += Bpp
        loc = self.setelement.lockey(self.nodelist)
        epsilon = B@(self.U[loc])
        strn_elm_xx = epsilon[0]
        strn_elm_yy = epsilon[1]
        strn_elm_zz = epsilon[2]
        strn_elm_xy = epsilon[3]
        strn_elm_yz = epsilon[4]
        strn_elm_zx = epsilon[5]
        strn_elm_eqv = np.sqrt(0.5*((epsilon[0]-epsilon[1])**2 + (epsilon[1] - epsilon[2])**2 + (
            epsilon[2]-epsilon[0])**2 + 6*(epsilon[3]**2 + epsilon[4]**2 + epsilon[5]**2)))
        strain = [strn_elm_eqv, strn_elm_xx, strn_elm_yy,
                  strn_elm_zz, strn_elm_xy, strn_elm_yz, strn_elm_zx]
        title = ["STRAIN_VM", "STRAIN_XX", "STRAIN_YY",
                 "STRAIN_ZZ", "STRAIN_XY", "STRAIN_YZ", "STRAIN_ZX"]
        return epsilon, strain, title

    def stress(self, epsilon):
        M = Elasticity(self.tabmat, self.inci, self.ee)
        D = M.isotropic()
        sigma = np.dot(D, epsilon)
        strs_elm_xx = sigma[0]
        strs_elm_yy = sigma[1]
        strs_elm_zz = sigma[2]
        strs_elm_xy = sigma[3]
        strs_elm_yz = sigma[4]
        strs_elm_zx = sigma[5]
        strs_elm_eqv = np.sqrt(0.5*((sigma[0]-sigma[1])**2 + (sigma[1]-sigma[2])**2 + (
            sigma[2]-sigma[0])**2 + 6*(sigma[3]**2 + sigma[4]**2 + sigma[5]**2)))
        stress = [strs_elm_eqv, strs_elm_xx, strs_elm_yy,
                  strs_elm_zz, strs_elm_xy, strs_elm_yz, strs_elm_zx]
        title = ["STRESS_VM", "STRESS_XX", "STRESS_YY",
                 "STRESS_ZZ", "STRESS_XY", "STRESS_YZ", "STRESS_ZX"]
        return stress, title
