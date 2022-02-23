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
# montagem matriz de rigidez trelica 2D
def rod2d_r2_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
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
        A = tabgeo[int(inci[ee,3]-1),0]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        
        T = np.zeros((dofe,dofe))
        T[0,0] = c
        T[0,1] = s
        T[1,0] = -s
        T[1,1] = c
        T[2,2] = c
        T[2,3] = s
        T[3,2] = -s
        T[3,3] = c
                
        ket2 = np.zeros((dofe,dofe))
        ket2[0,0] = 1.0
        ket2[0,2] = -1.0
        ket2[2,0] = -1.0
        ket2[2,2] = 1.0
        ket2 = ((E*A)/L)*ket2

        ket2T = np.dot(np.dot(np.transpose(T),ket2),T)
        
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = ket2T.flatten('F')
            
    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG


#%% POST-PROCESS
# malha deformada trelica
def rod2d_r2_deform(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-2,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-2,0]**2 + U[datamesh[0]*nn-1,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    return   meshDefU, Udef, Umag


# tensao no elemento de barra
def rod2d_r2_stress(datamesh,U,inci,coord,typeMechMat,tabmat):
    strs_elm = np.zeros((datamesh[2],1),dtype=float)
    for ii in range(datamesh[2]):
        loading_bar_v1(100*(ii/datamesh[2]))
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ii)
        E = D[0]
        G = D[1]
        noi = int(inci[ii,4])
        noj = int(inci[ii,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        T = np.array([[c,s,0,0],\
                      [0,0,c,s]])
        loc = np.array([2*noi-2,2*noi-1,2*noj-2,2*noj-1])
        
        B = (E/L)*np.array([-1,1])
        
        strs_elm[ii] = np.dot(B,T@U[loc,0])
    
    return strs_elm