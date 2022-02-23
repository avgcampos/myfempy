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
# montagem matriz de rigidez mola 0D
def spring_k2_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
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
        K = D[0]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        T = np.array([[c,s,0,0],[-s,c,0,0],[0,0,c,s],[0,0,-s,c]])
        
        kes0 = np.zeros((dofe,dofe))
        kes0[0,0] = 1.0
        kes0[0,2] = -1.0
        kes0[2,0] = -1.0
        kes0[2,2] = 1.0
        
        kes0 = K*kes0
        kes0T = np.dot(np.dot(np.transpose(T),kes0),T)
        
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = kes0T.flatten('F')
            
    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG

# montagem matriz de massa LUMBED mola 0D
def spring_k2_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
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
        
        mes0 =  D*np.eye(dofe)
        
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = mes0.flatten('F')
            
    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return MG

#%% POST-PROCESS