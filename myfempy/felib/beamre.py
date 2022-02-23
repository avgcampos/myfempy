# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: v15
_______________________________________________________________________________
 ~~~~~~~~~~ MODULO DE SIMULACAO PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~~~~

ESTE MODULO CONSTROE RESULTADOS PARA POS PROCESSAMENTO DA ANALISE 
- ENTRADAS: 
- SAIDA:     
===============================================================================

> ATUALIZACOES DA VERSAO:
    --- tensao normal em portico fream22
    --- tensao normal em vigas beam21
    --- atualizacoes de funcoes de acordo com o comportamento mecanico dos elementos
    --- malha deformada mais automatizada para frame
_______________________________________________________________________________
"""
import sys
import numpy as np
import scipy.sparse as sp
from scipy.linalg import block_diag
from myfempy.lib.material_behavior import mechanicalMaterial
from myfempy.lib.integration_points import pointGauss
from myfempy.setup.myfempy_preproc import loading_bar_v1

#------------------------------------------------------------------------------
# malha deformada trelica
def mesh_def_axial_rigid(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-2,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-2,0]**2 + U[datamesh[0]*nn-1,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    return   meshDefU, Udef, Umag


# malha deformada viga
def mesh_def_bending_2d(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Rdef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,1] = U[datamesh[0]*nn-2,0]
        Rdef[nn-1,2] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = U[datamesh[0]*nn-2,0]
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    meshRotZ = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Rdef)),axis=1)
    return  meshDefU,meshRotZ,Udef,Umag


# malha deformada frame 2D
def mesh_def_flexural_2d(datamesh,coord,U,scale_mesh):
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
    return  meshDefU,meshRotZ,Udef,Umag

# malha deformada frame 3D
def mesh_def_flexural_3d(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Rdef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-6,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-5,0]
        Udef[nn-1,2] = U[datamesh[0]*nn-4,0]
        Rdef[nn-1,0] = U[datamesh[0]*nn-3,0]
        Rdef[nn-1,1] = U[datamesh[0]*nn-2,0]
        Rdef[nn-1,2] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-6,0]**2+U[datamesh[0]*nn-5,0]**2 + U[datamesh[0]*nn-4,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    meshRotZ = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Rdef)),axis=1)
    return  meshDefU,meshRotZ,Udef,Umag

#------------------------------------------------------------------------------
# esfoco interno elemento viga (internal force - beam element)
def intn_force_bending_only(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U):
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

def intn_force_flexural_full_2d(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U):
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

def intn_force_flexural_full_3d(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U):
    dofe=2*datamesh[0]
    Fint=np.zeros((datamesh[3],1))
    Nx = np.zeros((datamesh[1],1),dtype=float)
    Vy = np.zeros((datamesh[1],1),dtype=float)
    Vz = np.zeros((datamesh[1],1),dtype=float)
    Tx = np.zeros((datamesh[1],1),dtype=float)
    My = np.zeros((datamesh[1],1),dtype=float)
    Mz = np.zeros((datamesh[1],1),dtype=float)
    domL = np.zeros((datamesh[1],1),dtype=float)
    for ee in range(datamesh[2]):
        noi = int(inci[ee,4])
        noj = int(inci[ee,5]) 
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        noiz = coord[noi-1,3]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        nojz = coord[noj-1,3]
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ee)
        E = D[0]
        G = D[1]
        A = tabgeo[int(inci[ee,3]-1),0]
        Izz = tabgeo[int(inci[ee,3]-1),1]
        Iyy = tabgeo[int(inci[ee,3]-1),2]
        Jxx = tabgeo[int(inci[ee,3]-1),3]
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
       
        kef3d2 = np.zeros((dofe,dofe))
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
                
        loc = np.array([datamesh[0]*noi-6,datamesh[0]*noi-5,datamesh[0]*noi-4,datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,\
                        datamesh[0]*noj-6,datamesh[0]*noj-5,datamesh[0]*noj-4,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1])

    
        F = np.dot(kef3d2T,T@U[loc,0])
        Fint[loc] = [-F[0],F[1],F[2],-F[3],-F[4],-F[5],F[6],-F[7],-F[8],F[9],F[10],F[11]]
   
    for i in range(1,datamesh[1]+1):
        Nx[i-1] = Fint[datamesh[0]*i-6]
        Vy[i-1] = Fint[datamesh[0]*i-5]
        Vz[i-1] = Fint[datamesh[0]*i-4]
        Tx[i-1] = Fint[datamesh[0]*i-3]
        My[i-1] = Fint[datamesh[0]*i-2]
        Mz[i-1] = Fint[datamesh[0]*i-1]
        domL[i-1] = coord[i-1,1]
    
    # domLOrd = np.sort(domL,axis=0)
    domLIdc = np.argsort(domL,axis=0)
    domL = domL[domLIdc[:,0]]
    Mz = Mz[domLIdc[:,0]]
    My = My[domLIdc[:,0]]
    Tx = Tx[domLIdc[:,0]]
    Vz = Vz[domLIdc[:,0]]
    Vy = Vy[domLIdc[:,0]]
    Nx = Nx[domLIdc[:,0]]
        
    return domL,Nx,Vy,Vz,Tx,My,Mz

# Reacoes de suporte/apoio
def react_suport(datamesh,coord,mtKG,U):
    Fr = np.zeros((datamesh[1],1),dtype=float)
    Mr = np.zeros((datamesh[1],1),dtype=float)
    React = np.dot(mtKG,U)
    
    for nn in range(1,datamesh[1]+1):
        Fr[nn-1] = React[datamesh[0]*nn-2]
        Mr[nn-1] = React[datamesh[0]*nn-1]
    
    # meshFr = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],Fr)),axis=1)
    # meshMr = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],Mr)),axis=1)
    return Fr,Mr

#------------------------------------------------------------------------------
# tensao no elemento de barra
def stress_axial_only_elm(datamesh,U,inci,coord,typeMechMat,tabmat):
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

# tensao no elemento de viga
def stress_bending_only_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max):
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


# tensao no elemento de portico 2d
def stress_flexural_full_2d_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max):
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


# tensao no elemento de portico 3d
def stress_flexural_full_3d_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max,z_max,r_max):
    strs_elm = np.zeros((datamesh[2],1),dtype=float)
    for ii in range(datamesh[2]):
        loading_bar_v1(100*(ii/datamesh[2]))
        noi = int(inci[ii,4])
        noj = int(inci[ii,5])
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        noiz = coord[noi-1,3]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        nojz = coord[noj-1,3]
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ii)
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
        
        loc = np.array([datamesh[0]*noi-6,datamesh[0]*noi-5,datamesh[0]*noi-4,datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,\
                        datamesh[0]*noj-6,datamesh[0]*noj-5,datamesh[0]*noj-4,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1])
            
        B1 = -E/L
        B2 = (-y_max*E)*((12*x_mid)/(L**3) - 6/(L**2))
        B3 = (-z_max*E)*((12*x_mid)/(L**3) - 6/(L**2))
        B4 = -G*r_max/L
        B5 = (-z_max*E)*((6*x_mid)/(L**2) - 4/L)
        B6 = (-y_max*E)*((6*x_mid)/(L**2) - 4/L)
        B7 = E/L
        B8 = (-y_max*E)*((-12*x_mid)/(L**3) + 6/(L**2))
        B9 = (-z_max*E)*((-12*x_mid)/(L**3) + 6/(L**2))
        B10 = G*r_max/L
        B11 = (-z_max*E)*((6*x_mid)/(L**2) - 2/L)
        B12 = (-y_max*E)*((6*x_mid)/(L**2) - 2/L)
        
        B = np.array([B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12])
        
        strs_elm[ii] = np.dot(B,T@U[loc,0])
        
    return strs_elm
