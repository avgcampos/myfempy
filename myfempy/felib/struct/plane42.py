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
from myfempy.lib.integration_points import pointGauss
from myfempy.setup.myfempy_preproc import loading_bar_v1

#%%------------------------------------------------------------------------------
# Q4 - isoparametric
def plastr_q4_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp):
    
    dofe=4*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    
    # A = np.zeros((datamesh[2],1))
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
        L = tabgeo[int(inci[ee,3]-1),4]
        
        xp,wp = pointGauss(npp)
        x = xp[0,:]
        y = xp[1,:]
        
        keq4 = np.zeros((dofe,dofe))
        for pp in range(0,npp):
            N1x = -(1-y[pp])
            N2x =  (1-y[pp])
            N3x =  (1+y[pp])
            N4x = -(1+y[pp])
            N1y = -(1-x[pp])
            N2y = -(1+x[pp])
            N3y =  (1+x[pp])
            N4y =  (1-x[pp])
            
            # dN = (1/4)*np.array([[-(1-y[pp]),(1-y[pp]),(1+y[pp]),-(1+y[pp])],\
            #                      [-(1-x[pp]),-(1+x[pp]),(1+x[pp]),(1-x[pp])]])
            dN = (1/4)*np.array([[N1x,N2x,N3x,N4x],[N1y,N2y,N3y,N4y]])
            
            J = np.dot(dN,matXY)
            
            detJ = J[0,0]*J[1,1] - J[1,0]*J[0,1]
            # detJ = np.linalg.det(J)
                        
            invJ = (1/detJ)*np.array([[J[1,1],-J[0,1],0.0,0.0],\
                                      [0.0,0.0,-J[1,0],J[0,0]],\
                                      [-J[1,0],J[0,0],J[1,1],-J[0,1]]])
            # invJ = np.linalg.inv(J)
                            
                        
            # pN = (1/4)*np.array([[-(1-y[pp]),0,(1-y[pp]),0,(1+y[pp]),0,-(1+y[pp]),0],\
            #                       [-(1-x[pp]),0,-(1+x[pp]),0,(1+x[pp]),0,(1-x[pp]),0],\
            #                       [0,-(1-y[pp]),0,(1-y[pp]),0,(1+y[pp]),0,-(1+y[pp])],\
            #                       [0,-(1-x[pp]),0,-(1+x[pp]),0,(1+x[pp]),0,(1-x[pp])]])
                                            
            
            
            pN = np.zeros((4,dofe))
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
        
            B = np.dot(invJ,pN)
            
            keq4 += L*np.dot(np.dot(np.transpose(B),D),B)*detJ*wp[pp]*wp[pp]
            
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1,datamesh[0]*nok-2,datamesh[0]*nok-1,datamesh[0]*nol-2,datamesh[0]*nol-1])
            
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = keq4.flatten('F')

    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG
    
# Q4 - isoparametric
def plastr_q4_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp):
    
    dofe=4*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    
    # A = np.zeros((datamesh[2],1))
    for ee in range(datamesh[2]): 
       
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
                
        L = tabgeo[int(inci[ee,3]-1),4]
        D = tabmat[int(inci[ee,2]-1),6]
        
        xp,wp = pointGauss(npp)
        x = xp[0,:]
        y = xp[1,:]
        
        meq4 = np.zeros((dofe,dofe))
        for pp in range(0,npp):
            
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
            
            detJ = J[0,0]*J[1,1] - J[1,0]*J[0,1]
            # detJ = np.linalg.det(J)
            
            meq4 += np.dot(np.transpose(N),N)*D*L*detJ*wp[pp]*wp[pp]
            
        loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-2,datamesh[0]*noj-1,datamesh[0]*nok-2,datamesh[0]*nok-1,datamesh[0]*nol-2,datamesh[0]*nol-1])
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = meq4.flatten('F')

    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return MG


#%% POST-PROCESS
# malha deformada plano 2D
def plastr_q4_deform(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-2,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-2,0]**2 + U[datamesh[0]*nn-1,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)    
    return  meshDefU,Udef,Umag

# tensao von-mises no elemento Plane
def plastr_q4_streqv(U,datamesh,inci,coord,typeMechMat,tabmat,npp):
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

