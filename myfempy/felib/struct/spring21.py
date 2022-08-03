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

# %%------------------------------------------------------------------------------


class Spring21:
    '''Spring 2D 2-node linear Finite Element'''

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

        dofelem = {'key': 'spring21',
                   'id': 110,
                   'def': 'struct 1D',
                   'dofs': ['ux', 'uy'],
                   'nnodes': ['i', 'j'],
                   'tensor': ['None']}

        return dofelem

    def lockey(self, list_node):

        noi = list_node[0]
        noj = list_node[1]

        loc = np.array([self.nodedof*noi-2, self.nodedof*noi-1,
                        self.nodedof*noj-2, self.nodedof*noj-1])

        return loc

    # montagem matriz de rigidez mola 0D
    def stiff_linear(self, ee):

        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])

        noix = self.coord[noi-1, 1]
        noiy = self.coord[noi-1, 2]
        nojx = self.coord[noj-1, 1]
        nojy = self.coord[noj-1, 2]

        D = get_elasticity(self.tabmat, self.inci, ee)
        S = D[0]

        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L

        T = np.zeros((self.dofe, self.dofe))
        T[0, 0] = c
        T[0, 1] = s
        T[1, 0] = -s
        T[1, 1] = c
        T[2, 2] = c
        T[2, 3] = s
        T[3, 2] = -s
        T[3, 3] = c

        kes0 = np.zeros((self.dofe, self.dofe))
        kes0[0, 0] = 1.0
        kes0[0, 2] = -1.0
        kes0[2, 0] = -1.0
        kes0[2, 2] = 1.0

        kes2 = S*kes0
        kes2T = np.dot(np.dot(np.transpose(T), kes2), T)

        list_node = [noi, noj]
        loc = Spring21.lockey(self, list_node)

        return kes2T, loc

    # montagem matriz de massa LUMBED mola 0D
    def mass(self, ee):

        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])

        R = self.tabmat[int(self.inci[ee, 2]-1), 6]

        mes2 = R*np.eye(self.dofe)

        list_node = [noi, noj]
        loc = Spring21.lockey(self, list_node)

        return mes2, loc

# %% POST-PROCESS
