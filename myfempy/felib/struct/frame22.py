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
# montagem matriz de rigidez de portico plana
def frame2d_f2_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
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
        Izz = tabgeo[int(inci[ee,3]-1),1]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
            
        T = np.zeros((dofe,dofe))
        T[0,0] = c
        T[0,1] = s
        T[1,0] = -s
        T[1,1] = c
        T[2,2] = 1.0
        T[3,3] = c
        T[3,4] = s
        T[4,3] = -s
        T[4,4] = c
        T[5,5] = 1.0
            
            
        kef2 = np.zeros((dofe,dofe))
        kef2[0,0] = (A*E)/L
        kef2[0,3] = -(A*E)/L
        
        kef2[1,1] = 12*(E*Izz)/L**3
        kef2[1,2] = 6*(E*Izz)/L**2
        kef2[1,4] = -12*(E*Izz)/L**3
        kef2[1,5] = 6*(E*Izz)/L**2
        
        kef2[2,1] = 6*(E*Izz)/L**2
        kef2[2,2] = 4*(E*Izz)/L
        kef2[2,4] = -6*(E*Izz)/L**2
        kef2[2,5] = 2*(E*Izz)/L
        
        kef2[3,0] = -(A*E)/L
        kef2[3,3] = (A*E)/L
        
        kef2[4,1] = -12*(E*Izz)/L**3
        kef2[4,2] = -6*(E*Izz)/L**2
        kef2[4,4] = 12*(E*Izz)/L**3
        kef2[4,5] = -6*(E*Izz)/L**2
        
        kef2[5,1] = 6*(E*Izz)/L**2
        kef2[5,2] = 2*(E*Izz)/L
        kef2[5,4] = -6*(E*Izz)/L**2
        kef2[5,5] = 4*(E*Izz)/L
        kef2 = (E*Izz/L**3)*kef2
            
        kef2T = np.dot(np.dot(np.transpose(T),kef2),T)
        
        loc = np.array([datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = kef2T.flatten('F')
            
    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG


# montagem matriz de massa de portico plana
def frame2d_f2_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
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
            
        mef2d2 = np.zeros((dofe,dofe))
        mef2d2[0,0] = 140
        mef2d2[0,3] = 70
        
        mef2d2[1,1] = 156
        mef2d2[1,2] = 22*L
        mef2d2[1,4] = 54
        mef2d2[1,5] = -13*L
        
        mef2d2[2,1] = 22*L
        mef2d2[2,2] = 4*L**2
        mef2d2[2,4] = 13*L
        mef2d2[2,5] = -3*L**2
        
        mef2d2[3,0] = 70
        mef2d2[3,3] = 140
        
        mef2d2[4,1] = 54
        mef2d2[4,2] = 13*L
        mef2d2[4,4] = 156
        mef2d2[4,5] = -22*L
        
        mef2d2[5,1] = -13*L
        mef2d2[5,2] = -3*L**2
        mef2d2[5,4] = -22*L
        mef2d2[5,5] = 4*L**2
        mef2d2 = (D*A*L/420)*mef2d2
   
        loc = np.array([datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = mef2d2.flatten('F')
            
    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return MG


#%% POST-PROCESS
# malha deformada frame 2D
def frame2d_f2_deform(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Rdef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-3,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-2,0]
        Rdef[nn-1,2] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-3,0]**2 + U[datamesh[0]*nn-2,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    meshRotZ = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Rdef)),axis=1)
    return  meshDefU, meshRotZ, Udef, Umag

def frame2d_f2_intforc(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U):
    Fint=np.zeros((datamesh[3],1))
    Nx = np.zeros((datamesh[1],1),dtype=float)
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
        A = tabgeo[int(inci[kk,3]-1),0]
        Izz = tabgeo[int(inci[kk,3]-1),1]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        T = np.array([[c,s,0,0,0,0],[-s,c,0,0,0,0],\
                      [0,0,1,0,0,0],[0,0,0,c,s,0],\
                      [0,0,0,-s,c,0],[0,0,0,0,0,1]])
        
        loc = np.array([3*noi-3,3*noi-2,3*noi-1,3*noj-3,3*noj-2,3*noj-1])
        kef2 = np.array([[(A*E)/L,0,0,-(A*E)/L,0,0],\
                         [0,12*(E*Izz)/L**3,6*(E*Izz)/L**2,0,-12*(E*Izz)/L**3,6*(E*Izz)/L**2],\
                         [0,6*(E*Izz)/L**2,4*(E*Izz)/L,0,-6*(E*Izz)/L**2,2*(E*Izz)/L],\
                         [-(A*E)/L,0,0,(A*E)/L,0,0],\
                         [0,-12*(E*Izz)/L**3,-6*(E*Izz)/L**2,0,12*(E*Izz)/L**3,-6*(E*Izz)/L**2],\
                         [0,6*(E*Izz)/L**2,2*(E*Izz)/L,0,-6*(E*Izz)/L**2,4*(E*Izz)/L]])
    
        F = np.dot(kef2,T@U[loc,0])
        Fint[loc] = [-F[0],F[1],-F[2],F[3],-F[4],F[5]]
    

    for i in range(1,datamesh[1]+1):
        Nx[i-1] = Fint[datamesh[0]*i-3]
        Vy[i-1] = Fint[datamesh[0]*i-2]
        Mz[i-1] = Fint[datamesh[0]*i-1]
        domL[i-1] = coord[i-1,1]
        
    
    BeamNx = np.zeros((datamesh[1],3),dtype=float)
    for ee in range(datamesh[2]):
        
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
                
        Lx = np.sqrt(((noix-nojx)**2))
        Ly = np.sqrt(((noiy-nojy)**2))
        
        if (Lx>0.0)and(Ly==0.0):
            BeamNx[noi-1,1] = Nx[noi-1]
            BeamNx[noj-1,1] = Nx[noj-1]
        if (Lx==0.0)and(Ly>0.0):
            BeamNx[noi-1,0] = Nx[noi-1]
            BeamNx[noj-1,0] = Nx[noj-1]
    
    meshBeamNx = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],BeamNx)),axis=1)
    
    BeamVy = np.zeros((datamesh[1],3),dtype=float)
    for ee in range(datamesh[2]):
        
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
                
        Lx = np.sqrt(((noix-nojx)**2))
        Ly = np.sqrt(((noiy-nojy)**2))
        
        if (Lx>0.0)and(Ly==0.0):
            BeamVy[noi-1,1] = Vy[noi-1]
            BeamVy[noj-1,1] = Vy[noj-1]
        if (Lx==0.0)and(Ly>0.0):
            BeamVy[noi-1,0] = Vy[noi-1]
            BeamVy[noj-1,0] = Vy[noj-1]
    
    meshBeamVy = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],BeamVy)),axis=1)
    
    BeamMz = np.zeros((datamesh[1],3),dtype=float)
    for ee in range(datamesh[2]):
        
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
                
        Lx = np.sqrt(((noix-nojx)**2))
        Ly = np.sqrt(((noiy-nojy)**2))
        
        if (Lx>0.0)and(Ly==0.0):
            BeamMz[noi-1,1] = Mz[noi-1]
            BeamMz[noj-1,1] = Mz[noj-1]
        if (Lx==0.0)and(Ly>0.0):
            BeamMz[noi-1,0] = Mz[noi-1]
            BeamMz[noj-1,0] = Mz[noj-1]
    
    meshBeamMz = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],BeamMz)),axis=1)
                
    domLOrd = np.sort(domL,axis=0)
    domLIdc = np.argsort(domL,axis=0)
    domL = domL[domLIdc[:,0]]
    Mz = Mz[domLIdc[:,0]]
    Vy = Vy[domLIdc[:,0]]
    Nx = Nx[domLIdc[:,0]]
    
    return domL,Nx,Vy,Mz


# tensao no elemento de portico 2d
def frame2d_f2_stress(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max):
    strs_elm = np.zeros((datamesh[2],1),dtype=float)
    for ii in range(datamesh[2]):
        loading_bar_v1(100*(ii/datamesh[2]))
        noi = int(inci[ii,4])
        noj = int(inci[ii,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ii)
        E = D[0]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2)
        
        s = (nojy-noiy)/L
        c = (nojx-noix)/L
        T = np.array([[c,s,0,0,0,0],\
                      [-s,c,0,0,0,0],\
                      [0,0,1,0,0,0],\
                      [0,0,0,c,s,0],\
                      [0,0,0,-s,c,0],\
                      [0,0,0,0,0,1]])
        
        coord_local = np.dot(T,np.array([[noix],[noiy],[1],[nojx],[nojy],[1]]))    
            
        x_mid = (coord_local[3][0]- coord_local[0][0])/2    
            
        loc = np.array([3*noi-3,3*noi-2,3*noi-1,3*noj-3,3*noj-2,3*noj-1])
        
        B = np.array([-E/L,-(y_max*E)*((12*x_mid)/(L**3) - 6/(L**2)),\
                    -(y_max*E)*((6*x_mid)/(L**2) - 4/L) ,\
                    E/L,-(y_max*E)*(-(12*x_mid)/(L**3) + 6/(L**2)),\
                    -(y_max*E)*((6*x_mid)/(L**2) - 2/L)])
  
        strs_elm[ii] = np.dot(B,T@U[loc,0]) 
        
    return strs_elm