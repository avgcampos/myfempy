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
from myfempy.felib.felemset import get_elemset
from myfempy.felib.crossec import cg_coord

# Material
#-----------------------------------------------------------------------------#


class Elasticity:
    '''Isotropic material'''

    def __init__(self, tabmat, inci, num_elm):

        self.E = tabmat[int(inci[num_elm, 2])-1, 0]  # material elasticity
        self.G = tabmat[int(inci[num_elm, 2])-1, 2]  # material shear

    def isotropic(self):

        E = self.E
        G = self.G
        D = [E, G]

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
        self.tabgeo = modelinfo['tabgeo']

        element = get_elemset(int(self.inci[ee][1]))
        self.setelement = element(self.modelinfo)

        # noi=int(self.inci[ee,4])
        # noj=int(self.inci[ee,5])
        # self.nodelist = [noi, noj]

        listnodes = self.inci[ee, 4:]
        self.nodelist = listnodes[listnodes.nonzero()].astype(int).tolist()

        cg = cg_coord(self.tabgeo, self.inci, ee)
        self.y_max = cg[0]
        self.y_min = cg[1]
        self.z_max = cg[2]
        self.z_min = cg[3]
        self.r_max = cg[4]

        # xp, wp = self.quadra
        self.xpp = [0.0]  # xp[0,:]
        self.ypp = [0.0]  # xp[1,:]
        self.zpp = [0.0]  # xp[2,:]
        self.wp = [2.0]

    def strain(self):

        csc_max = [self.y_max, self.z_max, self.r_max]
        Bmax, T = self.setelement.matrix_b(self.ee, csc_max)

        csc_min = [self.y_min, self.z_min, self.r_max]
        Bmin, T = self.setelement.matrix_b(self.ee, csc_min)

        loc = self.setelement.lockey(self.nodelist)

        epsilon_max = np.dot(Bmax, T@(self.U[loc]))
        epsilon_min = np.dot(Bmin, T@(self.U[loc]))
        strain = 0.0

        epsilon = [epsilon_max, epsilon_min]

        title = ["STRAIN_XX_MAX", "STRAIN_XX_MIN"]

        return epsilon, strain, title

    def stress(self, epsilon):

        stress_max = epsilon[0]
        stress_min = epsilon[1]

        stress = [stress_max, stress_min]

        title = ["STRESS_XX_MAX", "STRESS_XX_MIN"]

        return stress, title

 # class Ortotropic:
