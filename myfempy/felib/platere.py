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
from myfempy.lib.material_behavior import mechanicalMaterial
from myfempy.lib.integration_points import pointGauss
from myfempy.setup.myfempy_preproc import loading_bar_v1

#------------------------------------------------------------------------------
# malha deformada plano 2D
def mesh_def_plane_2d(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-2,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-2,0]**2 + U[datamesh[0]*nn-1,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    
    # Umod = np.zeros((datamesh[1],1),dtype=float)
    # for mm in range(0,datamesh[1]):
    #     Umod[mm,0] = np.sqrt((1/scale_mesh)*meshDefU[mm,1]**2 + (1/scale_mesh)*meshDefU[mm,2]**2)
    
    return  meshDefU,Udef,Umag

#------------------------------------------------------------------------------
# tensao von-mises no elemento Plane
def stress_vm_plane_t3(U,inci,coord,typeMechMat,tabmat,datamesh):
    stress_list=np.zeros((datamesh[2],1),dtype=float)
    strain_list=np.zeros((datamesh[2],1),dtype=float)
    for ee in range(datamesh[2]):
        loading_bar_v1(100*(ee/datamesh[2]))
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
        
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ee)    
        
        loc = np.array([2*noi-2,2*noi-1,2*noj-2,2*noj-1,2*nok-2,2*nok-1])
        
        strain =   B@U[loc,0]
        sigma =  np.dot(D,strain)
        # strs_elm[ee] = np.sqrt(sigma[0]**2 + sigma[1]**2 - sigma[0]*sigma[1] + 3*sigma[2]**2)
        stress_list[ee] = np.sqrt(sigma[0]**2 - sigma[0]*sigma[1] + sigma[1]**2 + 3*sigma[2]**2)
        strain_list[ee] = np.sqrt(strain[0]**2 - strain[0]*strain[1] + strain[1]**2 + 3*strain[2]**2)
    return stress_list, strain_list

#%% 
def stress_vm_plane_q4(U,datamesh,inci,coord,typeMechMat,tabmat,npp):
    strs_elm_eqv = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_xx = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_yy = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_xy = np.zeros((datamesh[2],1),dtype=float)
    
    strn_elm_eqv = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_xx = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_yy = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_xy = np.zeros((datamesh[2],1),dtype=float)
    for ee in range(datamesh[2]):
        loading_bar_v1(100*(ee/datamesh[2]))
        noi=int(inci[ee,4])
        noj=int(inci[ee,5])
        nok=int(inci[ee,6])
        nol=int(inci[ee,7])
        xi=coord[noi-1,1]
        yi=coord[noi-1,2]
        xj=coord[noj-1,1]
        yj=coord[noj-1,2]
        xk=coord[nok-1,1]
        yk=coord[nok-1,2]
        xl=coord[nol-1,1]
        yl=coord[nol-1,2]
        
        matXY = np.array([[xi,yi],[xj,yj],[xk,yk],[xl,yl]])
                
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ee)
        
        xp,wp = pointGauss(npp)
        x = xp[0,:]
        y = xp[1,:]
        
        B = np.zeros((3,8))
        for pp in range(0,npp):
            dN = (1/4)*np.array([[-(1-y[pp]),(1-y[pp]),(1+y[pp]),-(1+y[pp])],\
                                 [-(1-x[pp]),-(1+x[pp]),(1+x[pp]),(1-x[pp])]])
            J = np.dot(dN,matXY)
            
            detJ = J[0,0]*J[1,1] - J[1,0]*J[0,1]
                        
            invJ = (1/detJ)*np.array([[J[1,1],-J[0,1],0.0,0.0],\
                                      [0.0,0.0,-J[1,0],J[0,0]],\
                                      [-J[1,0],J[0,0],J[1,1],-J[0,1]]])
            
            pN = (1/4)*np.array([[-(1-y[pp]),0,(1-y[pp]),0,(1+y[pp]),0,-(1+y[pp]),0],\
                                 [-(1-x[pp]),0,-(1+x[pp]),0,(1+x[pp]),0,(1-x[pp]),0],\
                                 [0,-(1-y[pp]),0,(1-y[pp]),0,(1+y[pp]),0,-(1+y[pp])],\
                                 [0,-(1-x[pp]),0,-(1+x[pp]),0,(1+x[pp]),0,(1-x[pp])]])
                
            B += np.dot(invJ,pN)
        
        loc = np.array([2*noi-2,2*noi-1,2*noj-2,2*noj-1,2*nok-2,2*nok-1,2*nol-2,2*nol-1])
        
        strain =  B@U[loc,0]
        sigma =  np.dot(D,strain)
        
        strn_elm_xx[ee] = strain[0]
        strn_elm_yy[ee] = strain[1]
        strn_elm_xy[ee] = strain[2]
        
        strs_elm_xx[ee] = sigma[0]
        strs_elm_yy[ee] = sigma[1]
        strs_elm_xy[ee] = sigma[2]
        
        strs_elm_eqv[ee] = np.sqrt(sigma[0]**2 - sigma[0]*sigma[1] + sigma[1]**2 + 3*sigma[2]**2)
        strn_elm_eqv[ee] = np.sqrt(strain[0]**2 - strain[0]*strain[1] + strain[1]**2 + 3*strain[2]**2)
        
        stress_list = np.concatenate((strs_elm_eqv,strs_elm_xx,strs_elm_yy,strs_elm_xy),axis=1)
        strain_list = np.concatenate((strn_elm_eqv,strn_elm_xx,strn_elm_yy,strn_elm_xy),axis=1)
        
    return stress_list, strain_list

def stress_avr_plane(strs_elm,datamesh,inci):
    nnode = datamesh[1]
    nel = datamesh[2]
    ith = np.zeros((4*nel),dtype=int)
    jth = np.zeros((4*nel),dtype=int)
    val = np.zeros((4*nel),dtype=int)
    
    q0 = 0
    for i in range(datamesh[1]):
        
        elmlist = inci[(np.asarray(np.where(inci[:,4:]==i+1)))[0][:],0]
        q1 = elmlist.size
        ith[q0:q1+q0] = i
        jth[q0:q1+q0] = elmlist-1
        val[q0:q1+q0] = elmlist
        q0=q1+q0 
        
    S = sp.csc_matrix((val, (ith, jth)), shape=(nnode, nel))
    
    strs_avr = np.zeros((datamesh[1],1),dtype=float)
    for mm in range(datamesh[1]):
        strs_avr[mm] = np.mean(strs_elm[(sp.csr_matrix.nonzero(S[mm,:])[1])])
    return strs_avr
