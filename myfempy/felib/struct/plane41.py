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

#%%------------------------------------------------------------------------------

class Plane41:
    '''Quatrangular Isoparametric Plane 4-node linear Finite Element'''

    def __init__(self, modelinfo):
                              
        self.dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        self.nodecon = modelinfo['nodecon'][0]
        self.fulldof = modelinfo["nodedof"][0]*len(modelinfo["coord"])
        self.nodedof = modelinfo['nodedof'][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode =len(modelinfo["coord"])
        self.inci = modelinfo['inci']
        self.coord = modelinfo['coord']
        self.tabmat = modelinfo['tabmat']
        self.tabgeo = modelinfo['tabgeo']
        self.ntensor = modelinfo['ntensor'][0]
        
        if modelinfo['quadra'][0] == 1:
            self.npp = modelinfo['quadra'][1]
            self.quadra = Quadrature.gaussian(self.npp)

        elif modelinfo['quadra'][0] == 0:
            self.npp = modelinfo['quadra'][1]
            self.quadra = Quadrature.no_interpol(self.npp)
    
    
    def elemset():
    
        dofelem = {'key':'plane41',
                   'id': 220,
                   'def':'struct 2D',
                   'dofs':['ux', 'uy'],
                   'nnodes': ['i', 'j', 'k', 'l'],
                   'tensor': ['sxx', 'syy', 'sxy']}
    
        return dofelem
    
    
    
    def lockey(self, nodelist):
        
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
        
        loc = np.array([self.nodedof*noi-2,self.nodedof*noi-1,\
                        self.nodedof*noj-2,self.nodedof*noj-1,\
                        self.nodedof*nok-2,self.nodedof*nok-1,\
                        self.nodedof*nol-2,self.nodedof*nol-1])
        
        return loc
    
    def matriz_B(self, nodelist, intpl):
        
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
                        
        xi=self.coord[noi-1,1]
        yi=self.coord[noi-1,2]
        xj=self.coord[noj-1,1]
        yj=self.coord[noj-1,2]
        xk=self.coord[nok-1,1]
        yk=self.coord[nok-1,2]
        xl=self.coord[nol-1,1]
        yl=self.coord[nol-1,2]
        
        matXY = np.array([[xi,yi],[xj,yj],[xk,yk],[xl,yl]])
                    
        x = intpl[0]
        y = intpl[1]
        
        N1x = -(1-y)
        N2x =  (1-y)
        N3x =  (1+y)
        N4x = -(1+y)
        N1y = -(1-x)
        N2y = -(1+x)
        N3y =  (1+x)
        N4y =  (1-x)
        
        # dN = (1/4)*np.array([[-(1-y[pp]),(1-y[pp]),(1+y[pp]),-(1+y[pp])],\
        #                      [-(1-x[pp]),-(1+x[pp]),(1+x[pp]),(1-x[pp])]])
        
        dN = (1/4)*np.array([[N1x,N2x,N3x,N4x],
                             [N1y,N2y,N3y,N4y]])
        
        J = np.dot(dN,matXY)
        # detJ = np.linalg.det(J)
        # invJ = np.linalg.inv(J)
        
        detJ = J[0,0]*J[1,1] - J[1,0]*J[0,1]
                    
        invJ = (1/detJ)*np.array([[J[1,1],-J[0,1],0.0,0.0],\
                                  [0.0,0.0,-J[1,0],J[0,0]],\
                                  [-J[1,0],J[0,0],J[1,1],-J[0,1]]])
        # invJ = np.linalg.inv(J)
                        
                    
        # pN = (1/4)*np.array([[-(1-y[pp]),0,(1-y[pp]),0,(1+y[pp]),0,-(1+y[pp]),0],\
        #                       [-(1-x[pp]),0,-(1+x[pp]),0,(1+x[pp]),0,(1-x[pp]),0],\
        #                       [0,-(1-y[pp]),0,(1-y[pp]),0,(1+y[pp]),0,-(1+y[pp])],\
        #                       [0,-(1-x[pp]),0,-(1+x[pp]),0,(1+x[pp]),0,(1-x[pp])]])
                                        
        
        
        pN = np.zeros((self.nodecon,self.dofe))
        pN[0,0] = N1x
        pN[0,2] = N2x
        pN[0,4] = N3x
        pN[0,6] = N4x
    
        pN[1,0] = N1y
        pN[1,2] = N2y
        pN[1,4] = N3y
        pN[1,6] = N4y
        
        pN[2,1] = N1x
        pN[2,3] = N2x
        pN[2,5] = N3x
        pN[2,7] = N4x
        
        pN[3,1] = N1y
        pN[3,3] = N2y
        pN[3,5] = N3y
        pN[3,7] = N4y
        
        pN = (1/4)*pN
            
        B = np.dot(invJ, pN)
            
    # B += Bpp*wp[pp]*wp[pp]
    # A += detJ
        
        return B, detJ
        
    # montagem matriz triangular CST 2D Sparse
    def stiff_linear(self, ee):
        
        noi=int(self.inci[ee,4])
        noj=int(self.inci[ee,5])
        nok=int(self.inci[ee,6])
        nol=int(self.inci[ee,7])
        
        nodelist = [noi, noj, nok, nol]
        # self.nodelist = nodelist

        D = get_elasticity(self.tabmat,self.inci,ee)
        
        L = self.tabgeo[int(self.inci[ee,3]-1),4]
        
        xp, wp = self.quadra
        xpp = xp[0,:]
        ypp = xp[1,:]
                
        A = 0.0
        # B = np.zeros((3,8))
        keq4 = np.zeros((self.dofe, self.dofe))
        for pp in range(0, self.npp):
            
            # self.x = 
            # self.y = ypp[pp]
            intpl = [xpp[pp], ypp[pp]]
            
            Bpp, detJ = Plane41.matriz_B(self, nodelist, intpl)
            A += detJ
            keq4 += np.dot(np.dot(np.transpose(Bpp),D),Bpp)*L*detJ*wp[pp]*wp[pp]
            
        loc = Plane41.lockey(self, nodelist)
    
        return keq4, loc
    
    # montagem matriz triangular CST 2D Sparse
    def mass(self, ee):

        noi=int(self.inci[ee,4])
        noj=int(self.inci[ee,5])
        nok=int(self.inci[ee,6])
        nol=int(self.inci[ee,7])
        
        xi=self.coord[noi-1,1]
        yi=self.coord[noi-1,2]
        xj=self.coord[noj-1,1]
        yj=self.coord[noj-1,2]
        xk=self.coord[nok-1,1]
        yk=self.coord[nok-1,2]
        xl=self.coord[nol-1,1]
        yl=self.coord[nol-1,2]
        
        matXY = np.array([[xi,yi],[xj,yj],[xk,yk],[xl,yl]])
                
        R = self.tabmat[int(self.inci[ee,2]-1),6]
        
        L = self.tabgeo[int(self.inci[ee,3]-1),4]
        
        xp, wp = self.quadra
        x = xp[0,:]
        y = xp[1,:]
        
        meq4 = np.zeros((self.dofe, self.dofe))
        for pp in range(0,self.npp):
            
            N1 = (1-x[pp])*(1-y[pp])
            N2 = (1+x[pp])*(1-y[pp])
            N3 = (1+x[pp])*(1+y[pp])
            N4 = (1-x[pp])*(1+y[pp])
            
            N = (1/4)*np.array([[N1,0,N2,0,N3,0,N4,0],[0,N1,0,N2,0,N3,0,N4]])
                        
            N1x = -(1-y[pp])
            N2x =  (1-y[pp])
            N3x =  (1+y[pp])
            N4x = -(1+y[pp])
            N1y = -(1-x[pp])
            N2y = -(1+x[pp])
            N3y =  (1+x[pp])
            N4y =  (1-x[pp])           
            
            dN = (1/4)*np.array([[N1x,N2x,N3x,N4x],[N1y,N2y,N3y,N4y]])
            
            J = np.dot(dN,matXY)
            
            # detJ = J[0,0]*J[1,1] - J[1,0]*J[0,1]
            detJ = np.linalg.det(J)
            
            meq4 += np.dot(np.transpose(N),N)*R*L*detJ*wp[pp]*wp[pp]              
        
        list_node = [noi, noj, nok, nol]
        loc = Plane41.lockey(self, list_node)
        
        return meq4, loc
