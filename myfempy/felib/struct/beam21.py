# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: v20
_______________________________________________________________________________
 ~~~~~~~~~~        MATRIZ DE RIGIDEZ DOS ELEMENTOS FINITOS           ~~~~~~~~~~

ESTE MODULO DETERMINA A MATRIZ DE RIGIDEZ DO ELEMENTO FINITO E FAZ A MONTAGEM 
DA RIGIDEZ GLOBAL

> FUNCIONALIDADES
--- ENTRADAS: datamesh,inci,coord,tabmat,tabgeo
--- SAIDA: mtKG   

===============================================================================

> ATUALIZACOES DA VERSAO:
--- Matriz TRIAG30 sparce, mais compacta e menos espaco na memoria
--- Matriz de rigidez para elemento FRAME22: Rigidez Axial+Flexao
_______________________________________________________________________________
"""
import sys
import numpy as np
import scipy.sparse as sp
from myfempy.lib.material_behavior import mechanicalMaterial
from myfempy.setup.myfempy_preproc import loading_bar_v1
from scipy.linalg import block_diag

#%%------------------------------------------------------------------------------
# montagem matriz de rigidez de viga plana
def beameul_b2_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
    dofe=2*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    for ee in range(datamesh[2]):
        loading_bar_v1(100*(ee/datamesh[2]))
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ee)
        E = D[0]
        Izz = tabgeo[int(inci[ee,3]-1),1]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        
        keb1 = np.zeros((dofe,dofe))
        keb1[0,0] = 12
        keb1[0,1] = 6*L
        keb1[0,2] = -12
        keb1[0,3] = 6*L
        
        keb1[1,0] = 6*L
        keb1[1,1] = 4*L**2
        keb1[1,2] = -6*L
        keb1[1,3] = 2*L**2
        
        keb1[2,0] = -12
        keb1[2,1] = -6*L
        keb1[2,2] = 12
        keb1[2,3] = -6*L
        
        keb1[3,0] = 6*L
        keb1[3,1] = 2*L**2
        keb1[3,2] = -6*L
        keb1[3,3] = 4*L**2
        keb1 = (E*Izz/L**3)*keb1
        
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = keb1.flatten('F')
            
    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG


# montagem matriz de massa de viga plana
def beameul_b2_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
    dofe=2*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    for ee in range(datamesh[2]):
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        D = tabmat[int(inci[ee,2]-1),6]
        A = tabgeo[int(inci[ee,3]-1),0]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)

        meb2 = np.zeros((dofe,dofe))
        meb2[0,0] = 156
        meb2[0,1] = 22*L
        meb2[0,2] = 54
        meb2[0,3] = -13*L
        
        meb2[1,0] = 22*L
        meb2[1,1] = 4*L**2
        meb2[1,2] = 13*L
        meb2[1,3] = -3*L**2
        
        meb2[2,0] = 54
        meb2[2,1] = 13*L
        meb2[2,2] = 156
        meb2[2,3] = -22*L
        
        meb2[3,0] = -13*L
        meb2[3,1] = -3*L**2
        meb2[3,2] = -22*L
        meb2[3,3] = 4*L**2
        meb2 = (D*A*L/420)*meb2
            
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = meb2.flatten('F')
            
    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return MG


#%% POST-PROCESS
# malha deformada viga
def beameul_b2_deform(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Rdef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,1] = U[datamesh[0]*nn-2,0]
        Rdef[nn-1,2] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = U[datamesh[0]*nn-2,0]
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    meshRotZ = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Rdef)),axis=1)
    return  meshDefU, meshRotZ, Udef, Umag

# esfoco interno elemento viga (internal force - beam element)
def beameul_b2_intforc(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U):
    Fint=np.zeros((datamesh[3],1))
    Vy = np.zeros((datamesh[1],1),dtype=float)
    Mz = np.zeros((datamesh[1],1),dtype=float)
    domL = np.zeros((datamesh[1],1),dtype=float)
    for kk in range(datamesh[2]):
        noi = int(inci[kk,4])
        noj = int(inci[kk,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        D = mechanicalMaterial(typeMechMat,tabmat,inci,kk)
        E = D[0]
        Izz = tabgeo[int(inci[kk,3]-1),1]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        loc = np.array([2*noi-2,2*noi-1,2*noj-2,2*noj-1])
        kemI = (E*Izz/L**3)*np.array([[12,6*L,-12,6*L],\
                                    [6*L,4*L**2,-6*L,2*L**2],\
                                    [-12,-6*L,12,-6*L],\
                                    [6*L,2*L**2,-6*L,4*L**2]])
    
        F = np.dot(kemI,U[loc,0])
        Fint[loc] = [F[0],-F[1],-F[2],F[3]]
        
    
    for i in range(1,datamesh[1]+1):
        Vy[i-1] = Fint[datamesh[0]*i-2]
        Mz[i-1] = Fint[datamesh[0]*i-1]
        domL[i-1] = coord[i-1,1]

    domLOrd = np.sort(domL,axis=0)
    domLIdc = np.argsort(domL,axis=0)
    domL = domL[domLIdc[:,0]]
    Mz = Mz[domLIdc[:,0]]
    Vy = Vy[domLIdc[:,0]]

    return domL,Vy,Mz

# tensao no elemento de viga
def beameul_b2_stress(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max):
    strs_elm = np.zeros((datamesh[2],1),dtype=float)
    for ii in range(datamesh[2]):        
        loading_bar_v1(100*(ii/datamesh[2]))
        noi = int(inci[ii,4])
        noj = int(inci[ii,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        
        x_mid = (nojx-noix)/2
                
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ii)
        E = D[0]
        loc = np.array([2*noi-2,2*noi-1,2*noj-2,2*noj-1])
        
        B = -(y_max*E)*np.array([(12*x_mid)/(L**3) - 6/(L**2),\
                      (6*x_mid)/(L**2) - 4/L,\
                     -(12*x_mid)/(L**3) + 6/(L**2),\
                      (6*x_mid)/(L**2) - 2/L])
        
        strs_elm[ii] = B@U[loc,0]
    
    return strs_elm