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
from ..material import mech_mate
from ..integrat import gauss_point
from ...io.miscel import loading_bar_v1

#%%------------------------------------------------------------------------------
# montagem matriz triangular CST 2D Sparse
def plane_t3_stif(datamesh,inci,coord,tabmat,tabgeo):
    
    dofe=3*datamesh['nodedof']
    nel = datamesh['leninci']
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    
    loading_bar_v1(0,'ASBLK')
    for ee in range(nel):
        loading_bar_v1(100*((ee+1)/nel),'ASBLKG')
        noi=int(inci[ee,4])
        noj=int(inci[ee,5])
        nok=int(inci[ee,6])
        xi=coord[noi-1,1]
        yi=coord[noi-1,2]
        xj=coord[noj-1,1]
        yj=coord[noj-1,2]
        xk=coord[nok-1,1]
        yk=coord[nok-1,2]
        C = np.array([[1,xi,yi],[1,xj,yj],[1,xk,yk]])
        Ae = (1/2)*abs(np.linalg.det(C)) # finite element' area
        # A[ee,0] = Ae
        
        beti=yj-yk
        betj=yk-yi
        betk=yi-yj
        gami=xk-xj
        gamj=xi-xk
        gamk=xj-xi
                
        B = np.zeros((3,dofe))
        B[0,0] = beti
        B[0,2] = betj
        B[0,4] = betk
        B[1,1] = gami
        B[1,3] = gamj
        B[1,5] = gamk
        B[2,0] = gami
        B[2,1] = beti
        B[2,2] = gamj
        B[2,3] = betj
        B[2,4] = gamk
        B[2,5] = betk
        B = 0.5*Ae*B
        
        D = mech_mate(datamesh,tabmat,inci,ee)
        L = tabgeo[int(inci[ee,3]-1),4]
        
        ket3 = L*Ae*np.dot(np.dot(np.transpose(B),D),B)
        # mtKG[np.ix_(loc,loc)] += kem
        
        loc = np.array([datamesh['nodedof']*noi-2,datamesh['nodedof']*noi-1,datamesh['nodedof']*noj-2,datamesh['nodedof']*noj-1,datamesh['nodedof']*nok-2,datamesh['nodedof']*nok-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = ket3.flatten('F')
    
    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh['fulldof'], datamesh['fulldof']))
    # KG = sp.coo_matrix((val, (ith, jth)), shape=[datamesh[3], datamesh[3]])
    # KG = sp.csc_matrix(KG)
    return KG

# montagem matriz triangular CST 2D Sparse
def plane_t3_mass(datamesh,inci,coord,tabmat,tabgeo):
    
    dofe=3*datamesh['nodedof']
    nel = datamesh['leninci']
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    
    for ee in range(datamesh['leninci']):
        noi=int(inci[ee,4])
        noj=int(inci[ee,5])
        nok=int(inci[ee,6])
        xi=coord[noi-1,1]
        yi=coord[noi-1,2]
        xj=coord[noj-1,1]
        yj=coord[noj-1,2]
        xk=coord[nok-1,1]
        yk=coord[nok-1,2]
        C = np.array([[1,xi,yi],[1,xj,yj],[1,xk,yk]])
        Ae = (1/2)*abs(np.linalg.det(C)) # finite element' area
        L = tabgeo[int(inci[ee,3]-1),4]
        D = tabmat[int(inci[ee,2]-1),6]
        # A[ee,0] = Ae
        
        mat_aux = 2*np.eye(dofe)
        mat_aux[0,2] = 1
        mat_aux[0,4] = 1
        mat_aux[1,3] = 1
        mat_aux[1,5] = 1
        mat_aux[2,0] = 1
        mat_aux[2,4] = 1
        mat_aux[3,1] = 1
        mat_aux[3,5] = 1
        mat_aux[4,0] = 1
        mat_aux[4,2] = 1
        mat_aux[5,1] = 1
        mat_aux[5,3] = 1
        
        met3 = (D*L*Ae/12)*mat_aux
        # mtKG[np.ix_(loc,loc)] += kem
                
        loc = np.array([datamesh['nodedof']*noi-2,datamesh['nodedof']*noi-1,datamesh['nodedof']*noj-2,datamesh['nodedof']*noj-1,datamesh['nodedof']*nok-2,datamesh['nodedof']*nok-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = met3.flatten('F')
    
    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh['fulldof'], datamesh['fulldof']))
    return MG


#%% POST-PROCESS
# malha deformada plano 2D
def plane_t3_deform(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh['lencoord'],3),dtype=float)
    Umag = np.zeros((datamesh['lencoord'],1),dtype=float)
    for nn in range(1,datamesh['lencoord']+1):
        Udef[nn-1,0] = U[datamesh['nodedof']*nn-2,0]
        Udef[nn-1,1] = U[datamesh['nodedof']*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh['nodedof']*nn-2,0]**2 + U[datamesh['nodedof']*nn-1,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    return  meshDefU, Udef, Umag

# tensao von-mises no elemento Plane
def plane_t3_stress(U,inci,coord,tabmat,datamesh):
    strs_elm_vm = np.zeros((datamesh['leninci'],1),dtype=float)
    strs_elm_xx = np.zeros((datamesh['leninci'],1),dtype=float)
    strs_elm_yy = np.zeros((datamesh['leninci'],1),dtype=float)
    strs_elm_xy = np.zeros((datamesh['leninci'],1),dtype=float)
    
    strn_elm_vm = np.zeros((datamesh['leninci'],1),dtype=float)
    strn_elm_xx = np.zeros((datamesh['leninci'],1),dtype=float)
    strn_elm_yy = np.zeros((datamesh['leninci'],1),dtype=float)
    strn_elm_xy = np.zeros((datamesh['leninci'],1),dtype=float)
    loading_bar_v1(0,'STRESS')
    for ee in range(datamesh['leninci']):
        loading_bar_v1(100*((ee+1)/datamesh['leninci']),'STRESS')
        noi=int(inci[ee,4])
        noj=int(inci[ee,5])
        nok=int(inci[ee,6])
        
        xi=coord[noi-1,1]
        yi=coord[noi-1,2]
        xj=coord[noj-1,1]
        yj=coord[noj-1,2]
        xk=coord[nok-1,1]
        yk=coord[nok-1,2]
        C = np.array([[1,xi,yi],[1,xj,yj],[1,xk,yk]])
        Ae = (1/2)*abs(np.linalg.det(C)) # finite element' area
        
        beti=yj-yk
        betj=yk-yi
        betk=yi-yj
        gami=xk-xj
        gamj=xi-xk
        gamk=xj-xi
        B=(1/(2*Ae))*np.array([[beti,0,betj,0,betk,0],\
                               [0,gami,0,gamj,0,gamk],\
                               [gami,beti,gamj,betj,gamk,betk]])
        
        D = mech_mate(datamesh,tabmat,inci,ee)    
        
        loc = np.array([2*noi-2,2*noi-1,2*noj-2,2*noj-1,2*nok-2,2*nok-1])
        
        strain =   B@U[loc,0]
        sigma =  np.dot(D,strain)
        
        strn_elm_xx[ee] = strain[0]
        strn_elm_yy[ee] = strain[1]
        strn_elm_xy[ee] = strain[2]
        
        strs_elm_xx[ee] = sigma[0]
        strs_elm_yy[ee] = sigma[1]
        strs_elm_xy[ee] = sigma[2]
        
        strs_elm_vm[ee] = np.sqrt(sigma[0]**2 - sigma[0]*sigma[1] + sigma[1]**2 + 3*sigma[2]**2)
        strn_elm_vm[ee] = np.sqrt(strain[0]**2 - strain[0]*strain[1] + strain[1]**2 + 3*strain[2]**2)
        
        stress_list = np.concatenate((strs_elm_vm,strs_elm_xx,strs_elm_yy,strs_elm_xy),axis=1)
        strain_list = np.concatenate((strn_elm_vm,strn_elm_xx,strn_elm_yy,strn_elm_xy),axis=1)
        
    if datamesh['stresscal'] == 'avr':
        strs_avr_vm = calc_avr_plane(np.reshape(stress_list[:,0],(-1,1)),datamesh,inci)
        strs_avr_xx = calc_avr_plane(np.reshape(stress_list[:,1],(-1,1)),datamesh,inci)
        strs_avr_yy = calc_avr_plane(np.reshape(stress_list[:,2],(-1,1)),datamesh,inci)
        strs_avr_xy = calc_avr_plane(np.reshape(stress_list[:,3],(-1,1)),datamesh,inci)
        
        strn_avr_vm = calc_avr_plane(np.reshape(strain_list[:,0],(-1,1)),datamesh,inci)
        strn_avr_xx = calc_avr_plane(np.reshape(strain_list[:,1],(-1,1)),datamesh,inci)
        strn_avr_yy = calc_avr_plane(np.reshape(strain_list[:,2],(-1,1)),datamesh,inci)
        strn_avr_xy = calc_avr_plane(np.reshape(strain_list[:,3],(-1,1)),datamesh,inci)
        
        stress_list = np.concatenate((strs_avr_vm,strs_avr_xx,strs_avr_yy,strs_avr_xy),axis=1)
        strain_list = np.concatenate((strn_avr_vm,strn_avr_xx,strn_avr_yy,strn_avr_xy),axis=1)
    
    return stress_list, strain_list

def calc_avr_plane(data_elm,datamesh,inci):
    nnode = datamesh['lencoord']
    nel = datamesh['leninci']
    ith = np.zeros((4*nel),dtype=int)
    jth = np.zeros((4*nel),dtype=int)
    val = np.zeros((4*nel),dtype=int)
    
    q0 = 0
    for i in range(datamesh['lencoord']):
        
        elmlist = inci[(np.asarray(np.where(inci[:,4:]==i+1)))[0][:],0]
        q1 = elmlist.size
        ith[q0:q1+q0] = i
        jth[q0:q1+q0] = elmlist-1
        val[q0:q1+q0] = elmlist
        q0=q1+q0 
        
    S = sp.csc_matrix((val, (ith, jth)), shape=(nnode, nel))
    
    data_avr = np.zeros((datamesh['lencoord'],1),dtype=float)
    for mm in range(datamesh['lencoord']):
        data_avr[mm] = np.mean(data_elm[(sp.csr_matrix.nonzero(S[mm,:])[1])])
    return data_avr
