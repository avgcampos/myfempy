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

import sys
import numpy as np
import scipy.sparse as sp
from scipy.linalg import block_diag
from myfempy.felib.materset import get_elasticity
from myfempy.felib.crossec import cg_coord
from myfempy.tools.tools import loading_bar_v1

#%%------------------------------------------------------------------------------

class Frame22:
    '''Frame 3D 2-node linear Finite Element'''

    def __init__(self, modelinfo):
                              
        self.dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        self.fulldof = modelinfo["nodedof"][0]*len(modelinfo["coord"])
        self.nodedof = modelinfo['nodedof'][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode =len(modelinfo["coord"])
        self.inci = modelinfo['inci']
        self.coord = modelinfo['coord']
        self.tabmat = modelinfo['tabmat']
        self.tabgeo = modelinfo['tabgeo']
        
        
    def elemset():
        
        dofelem = {'key':'frame22',
                   'id': 141,
                   'def':'struct 1D',
                   'dofs': ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'],
                   'nnodes': ['i', 'j'],
                   'tensor': ['sxx']}
        
        return dofelem

    
    def lockey(self, list_node):
        
        noi = list_node[0]
        noj = list_node[1]
        
        loc = np.array([self.nodedof*noi-6, self.nodedof*noi-5, self.nodedof*noi-4, self.nodedof*noi-3, self.nodedof*noi-2, self.nodedof*noi-1,\
                        self.nodedof*noj-6, self.nodedof*noj-5, self.nodedof*noj-4, self.nodedof*noj-3, self.nodedof*noj-2, self.nodedof*noj-1])
            
        return loc
    
    
    # montagem matriz de rigidez de portico espacial
    def stiff_linear(self, ee):
                
        noi = int(self.inci[ee,4])
        noj = int(self.inci[ee,5]) 
        
        noix = self.coord[noi-1,1]
        noiy = self.coord[noi-1,2]
        noiz = self.coord[noi-1,3]
        nojx = self.coord[noj-1,1]
        nojy = self.coord[noj-1,2]
        nojz = self.coord[noj-1,3]
        
        D = get_elasticity(self.tabmat,self.inci,ee)
        E = D[0]
        G = D[1]
        
        A = self.tabgeo[int(self.inci[ee,3]-1),0]
        Izz = self.tabgeo[int(self.inci[ee,3]-1),1]
        Iyy = self.tabgeo[int(self.inci[ee,3]-1),2]
        Jxx = self.tabgeo[int(self.inci[ee,3]-1),3]
        
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2 + (nojz-noiz)**2)
        
        if (noix == nojx) and (noiy == nojy):
            if nojz > noiz:
                lamb = np.array([[0, 0, 1.0],[0, 1.0, 0],[-1.0, 0, 0]])
            else:
                lamb = np.array([[0, 0, -1.0],[0, 1.0, 0],[1.0, 0, 0]])
        else:
            l = (nojx-noix)/L
            m = (nojy-noiy)/L
            n = (nojz-noiz)/L
            d = np.sqrt(l**2 + m**2)
            lamb = np.array([[l,m,n],[-m/d,l/d,0],[-l*n/d,-m*n/d,d]])
        
        T = block_diag(lamb, lamb, lamb, lamb) 
        
        kef3d2 = np.zeros((self.dofe,self.dofe))
        kef3d2[0,0] = A*E/L
        kef3d2[0,6] = -A*E/L
        kef3d2[6,0] = -A*E/L
        kef3d2[6,6] = A*E/L
        
        kef3d2[1,1] = 12*E*Izz/L**3
        kef3d2[1,7] = -12*E*Izz/L**3
        kef3d2[2,2] = 12*E*Iyy/L**3
        kef3d2[2,8] = -12*E*Iyy/L**3
        kef3d2[7,1] = -12*E*Izz/L**3
        kef3d2[7,7] = 12*E*Izz/L**3
        kef3d2[8,2] = -12*E*Iyy/L**3
        kef3d2[8,8] = 12*E*Iyy/L**3
        
        kef3d2[1,5] = 6*E*Izz/L**2
        kef3d2[1,11] = 6*E*Izz/L**2
        kef3d2[2,4] = -6*E*Iyy/L**2
        kef3d2[2,10] = -6*E*Iyy/L**2
        kef3d2[4,2] = -6*E*Iyy/L**2
        kef3d2[4,8] = 6*E*Iyy/L**2
        kef3d2[5,1] = 6*E*Izz/L**2
        kef3d2[5,7] = -6*E*Izz/L**2
        
        kef3d2[7,5] = -6*E*Izz/L**2
        kef3d2[7,11] = -6*E*Izz/L**2
        kef3d2[8,4] = 6*E*Iyy/L**2
        kef3d2[8,10] = 6*E*Iyy/L**2
        kef3d2[10,2] = -6*E*Iyy/L**2
        kef3d2[10,8] = 6*E*Iyy/L**2
        kef3d2[11,1] = 6*E*Izz/L**2
        kef3d2[11,7] = -6*E*Izz/L**2
     
        kef3d2[4,4] = 4*E*Iyy/L
        kef3d2[4,10] = 2*E*Iyy/L
        kef3d2[5,5] = 4*E*Izz/L
        kef3d2[5,11] = 2*E*Izz/L
        kef3d2[10,4] = 2*E*Iyy/L
        kef3d2[10,10] = 4*E*Iyy/L
        kef3d2[11,5] = 2*E*Izz/L
        kef3d2[11,11] = 4*E*Izz/L
        
        kef3d2[3,3] = G*Jxx/L
        kef3d2[3,9] = -G*Jxx/L
        kef3d2[9,3] = -G*Jxx/L
        kef3d2[9,9] = G*Jxx/L

        kef3d2T = np.dot(np.dot(np.transpose(T),kef3d2),T)
        
            
        list_node = [noi,noj]
        loc = Frame22.lockey(self, list_node)
        
        return kef3d2T, loc
    
    # montagem matriz de rigidez de portico espacial
    def mass(self, ee):
            
        noi = int(self.inci[ee,4])
        noj = int(self.inci[ee,5]) 
        
        noix = self.coord[noi-1,1]
        noiy = self.coord[noi-1,2]
        noiz = self.coord[noi-1,3]
        nojx = self.coord[noj-1,1]
        nojy = self.coord[noj-1,2]
        nojz = self.coord[noj-1,3]
        
        R = self.tabmat[int(self.inci[ee,2]-1),6]
        
        A = self.tabgeo[int(self.inci[ee,3]-1),0]
        Jxx = self.tabgeo[int(self.inci[ee,3]-1),3]
        
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2 + (nojz-noiz)**2)
        
        mef3d = np.zeros((self.dofe,self.dofe))
        mef3d[0,0] = 1/3
        mef3d[6,0] = 1/6
        mef3d[1,1] = 13/35
        mef3d[5,1] = (11*L)/210
        mef3d[7,1] = 9/70
        mef3d[11,1] = (-13*L)/420
        mef3d[2,2] = 13/35
        mef3d[4,2] = (-11*L)/210
        mef3d[8,2] = 9/70
        mef3d[10,2] = (13*L)/420
        mef3d[3,3] = Jxx/(3*A)
        mef3d[9,3] = Jxx/(6*A)
        mef3d[4,4] = (L**2)/105
        mef3d[8,4] = (-13*L)/420
        mef3d[10,4] = (-L**2)/140
        mef3d[5,5] = (L**2)/105
        mef3d[7,5] = (13*L)/420
        mef3d[11,5] = (-L**2)/140
        mef3d[6,6] = 1/3
        mef3d[7,7] = 13/35
        mef3d[11,7] = (-11*L)/210
        mef3d[8,8] = 13/35
        mef3d[10,8] = (11*L)/210
        mef3d[9,9] = Jxx/(3*A)
        mef3d[10,10] = (L**2)/105
        mef3d[11,11] = (L**2)/105
    
        mef3d[0,6]  = mef3d[6,0]
        mef3d[1,5]  = mef3d[5,1]
        mef3d[1,7]  = mef3d[7,1]
        mef3d[1,11] = mef3d[11,1]
        mef3d[2,4]  = mef3d[4,2]
        mef3d[2,8]  = mef3d[8,2]
        mef3d[2,10] = mef3d[10,2]
        mef3d[3,9]  = mef3d[9,3]
        mef3d[4,8]  = mef3d[8,4]
        mef3d[4,10] = mef3d[10,4]
        mef3d[5,7]  = mef3d[7,5]
        mef3d[5,11] = mef3d[11,5]
        mef3d[7,11] = mef3d[11,7]
        mef3d[8,10] = mef3d[10,8]
        
        mef3d = (R*A*L)*mef3d
        
        list_node = [noi,noj]
        loc = Frame22.lockey(self, list_node)
            
        return mef3d, loc
                
    
    
    def intforces(self, U, lines):

        Fint=np.zeros((self.fulldof),dtype=float)
        Nx = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        Vy = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        Vz = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        Tx = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        My = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        Mz = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        domL = np.zeros((len(lines[0][1]),len(lines)),dtype=float)
        for ee in range(self.nelem):
            
            noi = int(self.inci[ee,4])
            noj = int(self.inci[ee,5]) 
            noix = self.coord[noi-1,1]
            noiy = self.coord[noi-1,2]
            noiz = self.coord[noi-1,3]
            nojx = self.coord[noj-1,1]
            nojy = self.coord[noj-1,2]
            nojz = self.coord[noj-1,3]
                        
            L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2 + (nojz-noiz)**2)
                        
            l = (nojx-noix)/L
            m = (nojy-noiy)/L
            n = (nojz-noiz)/L
            d = np.sqrt(l**2 + n**2)
            
            if d == 0.0:
                lamb = np.array([[0,m,0],[-m,0,0],[0,0,1]])
            else:
                lamb = np.array([[l,m,n],[-l*m/d,(l**2 + n**2)/d,-m*n/d],[-n/d,0,l/d]])
            
            T = block_diag(lamb, lamb, lamb, lamb) 
           
            kef3d2T, loc = Frame22.stiff_linear(self, ee)
          
            F = np.dot(kef3d2T,T@U[loc])
            Fint[loc] = [-F[0],F[1],F[2],-F[3],-F[4],-F[5],F[6],-F[7],-F[8],F[9],F[10],F[11]]


        for ed in range(len(lines)):
            nnodes = lines[ed][1]
            
            for nn in range(0,len(nnodes)):
                node = int(nnodes[nn])
                Nx[nn,ed] = Fint[self.nodedof*node-6]
                Vy[nn,ed] = Fint[self.nodedof*node-5]
                Vz[nn,ed] = Fint[self.nodedof*node-4]
                Tx[nn,ed] = Fint[self.nodedof*node-3]
                My[nn,ed] = Fint[self.nodedof*node-2]
                Mz[nn,ed] = Fint[self.nodedof*node-1]
                domL[nn,ed] = self.coord[node-1,1]
                
            domLIdc = np.argsort(domL[:,ed],axis=0)
            domL[:,ed] = domL[domLIdc,ed]
            Nx[:,ed] = Nx[domLIdc,ed]
            Vy[:,ed] = Vy[domLIdc,ed]
            Vz[:,ed] = Vz[domLIdc,ed]
            Tx[:,ed] = Tx[domLIdc,ed]
            My[:,ed] = My[domLIdc,ed]
            Mz[:,ed] = Mz[domLIdc,ed]
        
        ifb = {'le':domL, 'val': [Nx, Vy, Vz, Tx, My, Mz]}
        title =  ["NX","VY","VZ", "TX", "MY", "MZ"]
            
        return ifb, title
    
    # # tensao no elemento de portico 3d
    def matrix_B(self, ee, csc):
        
        y = csc[0]
        z = csc[1]
        r = csc[2]
    
        noi = int(self.inci[ee,4])
        noj = int(self.inci[ee,5])
        # list_node = [noi,noj]
        
        noix = self.coord[noi-1,1]
        noiy = self.coord[noi-1,2]
        noiz = self.coord[noi-1,3]
        nojx = self.coord[noj-1,1]
        nojy = self.coord[noj-1,2]
        nojz = self.coord[noj-1,3]
        
        D = get_elasticity(self.tabmat,self.inci,ee)
        E = D[0]
        G = D[1]
        
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2 + (nojz-noiz)**2)
        
        l = (nojx-noix)/L
        m = (nojy-noiy)/L
        n = (nojz-noiz)/L
        d = np.sqrt(l**2 + n**2)
        
        if d == 0.0:
            lamb = np.array([[0,m,0],[-m,0,0],[0,0,1]])
        else:
            lamb = np.array([[l,m,n],[-l*m/d,(l**2 + n**2)/d,-m*n/d],[-n/d,0,l/d]])
        
        T = block_diag(lamb, lamb, lamb, lamb)
        
        coord_local = np.dot(T[0:6,0:6],np.array([[noix],[noiy],[noiz],[nojx],[nojy],[nojz]]))
        x_mid = (coord_local[3][0]-coord_local[0][0])/2 
                            
        B1 = -E/L
        B2 = (-y*E)*((12*x_mid)/(L**3) - 6/(L**2))
        B3 = (-z*E)*((12*x_mid)/(L**3) - 6/(L**2))
        B4 = -G*r/L
        B5 = (-z*E)*((6*x_mid)/(L**2) - 4/L)
        B6 = (-y*E)*((6*x_mid)/(L**2) - 4/L)
        B7 = E/L
        B8 = (-y*E)*((-12*x_mid)/(L**3) + 6/(L**2))
        B9 = (-z*E)*((-12*x_mid)/(L**3) + 6/(L**2))
        B10 = G*r/L
        B11 = (-z*E)*((6*x_mid)/(L**2) - 2/L)
        B12 = (-y*E)*((6*x_mid)/(L**2) - 2/L)
        
        B = np.array([B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12])
                   
        return B, T