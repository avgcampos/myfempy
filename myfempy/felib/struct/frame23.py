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
# montagem matriz de rigidez de portico espacial
def frame3d_f2_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
    dofe=2*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    # mtKG = np.zeros((datamesh[3],datamesh[3]))
    for ee in range(datamesh[2]):
        loading_bar_v1(100*(ee/datamesh[2]))
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
        
        if (noix == nojx) and (noiy == nojy):
            if nojz > noiz:
                lamb = np.array([[0,0,1.0],[0,1.0,0],[-1.0,0,0]])
            else:
                lamb = np.array([[0,0,-1.0],[0,1.0,0],[1.0,0,0]])
        else:
            l = (nojx-noix)/L
            m = (nojy-noiy)/L
            n = (nojz-noiz)/L
            d = np.sqrt(l**2 + m**2)
            lamb = np.array([[l,m,n],[-m/d,l/d,0],[-l*n/d,-m*n/d,d]])
        
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
        kef3d2T = np.dot(np.dot(np.transpose(T),kef3d2),T)
        
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = kef3d2T.flatten('F')
            
    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG

# montagem matriz de rigidez de portico espacial
def frame3d_f2_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo):
    dofe=2*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    # mtKG = np.zeros((datamesh[3],datamesh[3]))
    for ee in range(datamesh[2]):
        noi = int(inci[ee,4])
        noj = int(inci[ee,5]) 
        noix = coord[noi-1,1]
        noiy = coord[noi-1,2]
        noiz = coord[noi-1,3]
        nojx = coord[noj-1,1]
        nojy = coord[noj-1,2]
        nojz = coord[noj-1,3]
        D = tabmat[int(inci[ee,2]-1),6]
        A = tabgeo[int(inci[ee,3]-1),0]
        Jxx = tabgeo[int(inci[ee,3]-1),3]
        L = np.sqrt((nojx-noix)**2 + (nojy-noiy)**2 + (nojz-noiz)**2)
        
        mef3d = np.zeros((dofe,dofe))
        mef3d[0,0] = 1/3
        mef3d[6,0] = 1/6
        mef3d[1,1] = 13/35
        mef3d[5,1] = 11*L/210
        mef3d[7,1] = 9/70
        mef3d[11,1] = -13*L/420
        mef3d[2,2] = 13/35
        mef3d[4,2] = -11*L/210
        mef3d[8,2] = 9/70
        mef3d[10,2] = 13*L/420
        mef3d[3,3] = Jxx/3*A
        mef3d[9,3] = Jxx/6*A
        mef3d[4,4] = (L**2)/105
        mef3d[8,4] = -13*L/420
        mef3d[10,4] = (-L**2)/140
        mef3d[5,5] = (L**2)/105
        mef3d[7,5] = 13*L/420
        mef3d[11,5] = (-L**2)/140
        mef3d[6,6] = 1/3
        mef3d[7,7] = 13/35
        mef3d[11,7] = -11*L/210
        mef3d[8,8] = 13/35
        mef3d[10,8] = 11*L/210
        mef3d[9,9] = Jxx/3*A
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
        
        mef3d = (D*A*L)*mef3d
        
        loc = np.array([datamesh[0]*noi-6,datamesh[0]*noi-5,datamesh[0]*noi-4,datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,\
                        datamesh[0]*noj-6,datamesh[0]*noj-5,datamesh[0]*noj-4,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = mef3d.flatten('F')
            
    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return MG
            

#%% POST-PROCESS
# malha deformada frame 3D
def frame3d_f2_deform(datamesh,coord,U,scale_mesh):
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
    return  meshDefU, meshRotZ, Udef, Umag

def frame3d_f2_intforc(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U):
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

# tensao no elemento de portico 3d
def frame3d_f2_stress(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max,z_max,r_max):
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