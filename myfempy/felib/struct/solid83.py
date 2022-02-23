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
# LINEAR HEXAHEDRAL ELEMENT - isoparametric
def solid_hx8_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp):
    
    dofe = 8*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    
    # A = np.zeros((datamesh[2],1))
    for ee in range(datamesh[2]): 
        loading_bar_v1(100*(ee/datamesh[2]))
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        nok = int(inci[ee,6])
        nol = int(inci[ee,7])
        nom = int(inci[ee,8])
        non = int(inci[ee,9])
        noo = int(inci[ee,10])
        nop = int(inci[ee,11])
        
        xi = coord[noi-1,1]
        yi = coord[noi-1,2]
        zi = coord[noi-1,3]        
        xj = coord[noj-1,1]
        yj = coord[noj-1,2]
        zj = coord[noj-1,3]        
        xk = coord[nok-1,1]
        yk = coord[nok-1,2]
        zk = coord[nok-1,3]        
        xl = coord[nol-1,1]
        yl = coord[nol-1,2]
        zl = coord[nol-1,3]        
        xm = coord[nom-1,1]
        ym = coord[nom-1,2]
        zm = coord[nom-1,3]        
        xn = coord[non-1,1]
        yn = coord[non-1,2]
        zn = coord[non-1,3]        
        xo = coord[noo-1,1]
        yo = coord[noo-1,2]
        zo = coord[noo-1,3]        
        xp = coord[nop-1,1]
        yp = coord[nop-1,2]
        zp = coord[nop-1,3]        
        
        matXY = np.array([[xi,yi,zi],[xj,yj,zj],[xk,yk,zk],[xl,yl,zl],\
                          [xm,ym,zm],[xn,yn,zn],[xo,yo,zo],[xp,yp,zp]])
                
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ee)
        
        xp,wp = pointGauss(npp)
        x = xp[0,:]
        y = xp[1,:]
        z = xp[2,:]
        
        keh8 = np.zeros((dofe,dofe))
        for pp in range(0,npp):
            N1x = -(1-y[pp])*(1-z[pp])
            N2x =  (1-y[pp])*(1-z[pp])
            N3x =  (1+y[pp])*(1-z[pp])
            N4x = -(1+y[pp])*(1-z[pp])
            N5x = -(1-y[pp])*(1+z[pp])
            N6x =  (1-y[pp])*(1+z[pp])
            N7x =  (1+y[pp])*(1+z[pp])
            N8x = -(1+y[pp])*(1+z[pp]) 
            N1y = -(1-x[pp])*(1-z[pp])
            N2y = -(1+x[pp])*(1-z[pp])
            N3y =  (1+x[pp])*(1-z[pp])
            N4y =  (1-x[pp])*(1-z[pp])
            N5y = -(1-x[pp])*(1+z[pp])
            N6y = -(1+x[pp])*(1+z[pp])
            N7y =  (1+x[pp])*(1+z[pp])
            N8y =  (1-x[pp])*(1+z[pp])
            N1z = -(1-x[pp])*(1-y[pp])
            N2z = -(1+x[pp])*(1-y[pp])
            N3z = -(1+x[pp])*(1+y[pp])
            N4z = -(1-x[pp])*(1+y[pp])
            N5z =  (1-x[pp])*(1-y[pp])
            N6z =  (1+x[pp])*(1-y[pp])
            N7z =  (1+x[pp])*(1+y[pp])
            N8z =  (1-x[pp])*(1+y[pp])
            
            dN = (1/8)*np.array([[N1x,N2x,N3x,N4x,N5x,N6x,N7x,N8x],\
                                 [N1y,N2y,N3y,N4y,N5y,N6y,N7y,N8y],\
                                 [N1z,N2z,N3z,N4z,N5z,N6z,N7z,N8z]])
            
            J = np.dot(dN,matXY)
            
            detJ = np.linalg.det(J)
            
            invJ = np.linalg.inv(J)
            
            dN1 = np.array([[N1x],[N1y],[N1z]])
            dN2 = np.array([[N2x],[N2y],[N2z]])
            dN3 = np.array([[N3x],[N3y],[N3z]])
            dN4 = np.array([[N4x],[N4y],[N4z]])
            dN5 = np.array([[N5x],[N5y],[N5z]])
            dN6 = np.array([[N6x],[N6y],[N6z]])
            dN7 = np.array([[N7x],[N7y],[N7z]])
            dN8 = np.array([[N8x],[N8y],[N8z]])
            
            pN1 = np.dot(invJ,dN1)
            pN2 = np.dot(invJ,dN2)
            pN3 = np.dot(invJ,dN3)
            pN4 = np.dot(invJ,dN4)
            pN5 = np.dot(invJ,dN5)
            pN6 = np.dot(invJ,dN6)
            pN7 = np.dot(invJ,dN7)
            pN8 = np.dot(invJ,dN8)
            
            B1 = np.array([[pN1[0,0],0,0],\
                          [0,pN1[1,0],0],\
                          [0,0,pN1[2,0]],\
                          [pN1[1,0],pN1[0,0],0],\
                          [0,pN1[2,0],pN1[1,0]],\
                          [pN1[2,0],0,pN1[0,0]]])
           

            B2 = np.array([[pN2[0,0],0,0],\
                          [0,pN2[1,0],0],\
                          [0,0,pN2[2,0]],\
                          [pN2[1,0],pN2[0,0],0],\
                          [0,pN2[2,0],pN2[1,0]],\
                          [pN2[2,0],0,pN2[0,0]]])
                
            B3 = np.array([[pN3[0,0],0,0],\
                          [0,pN3[1,0],0],\
                          [0,0,pN3[2,0]],\
                          [pN3[1,0],pN3[0,0],0],\
                          [0,pN3[2,0],pN3[1,0]],\
                          [pN3[2,0],0,pN3[0,0]]])
                
            B4 = np.array([[pN4[0,0],0,0],\
                          [0,pN4[1,0],0],\
                          [0,0,pN4[2,0]],\
                          [pN4[1,0],pN4[0,0],0],\
                          [0,pN4[2,0],pN4[1,0]],\
                          [pN4[2,0],0,pN4[0,0]]])

            B5 = np.array([[pN5[0,0],0,0],\
                          [0,pN5[1,0],0],\
                          [0,0,pN5[2,0]],\
                          [pN5[1,0],pN5[0,0],0],\
                          [0,pN5[2,0],pN5[1,0]],\
                          [pN5[2,0],0,pN5[0,0]]])
           

            B6 = np.array([[pN6[0,0],0,0],\
                          [0,pN6[1,0],0],\
                          [0,0,pN6[2,0]],\
                          [pN6[1,0],pN6[0,0],0],\
                          [0,pN6[2,0],pN6[1,0]],\
                          [pN6[2,0],0,pN6[0,0]]])
                
            B7 = np.array([[pN7[0,0],0,0],\
                          [0,pN7[1,0],0],\
                          [0,0,pN7[2,0]],\
                          [pN7[1,0],pN7[0,0],0],\
                          [0,pN7[2,0],pN7[1,0]],\
                          [pN7[2,0],0,pN7[0,0]]])
                
            B8 = np.array([[pN8[0,0],0,0],\
                          [0,pN8[1,0],0],\
                          [0,0,pN8[2,0]],\
                          [pN8[1,0],pN8[0,0],0],\
                          [0,pN8[2,0],pN8[1,0]],\
                          [pN8[2,0],0,pN8[0,0]]])
            
            
            B = (1/8)*np.concatenate((B1,B2,B3,B4,B5,B6,B7,B8),axis=1)
           
            keh8 += np.dot(np.dot(np.transpose(B),D),B)*detJ*wp[pp]*wp[pp]
            
        loc = np.array([datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1,\
                        datamesh[0]*nok-3,datamesh[0]*nok-2,datamesh[0]*nok-1,datamesh[0]*nol-3,datamesh[0]*nol-2,datamesh[0]*nol-1,\
                        datamesh[0]*nom-3,datamesh[0]*nom-2,datamesh[0]*nom-1,datamesh[0]*non-3,datamesh[0]*non-2,datamesh[0]*non-1,\
                        datamesh[0]*noo-3,datamesh[0]*noo-2,datamesh[0]*noo-1,datamesh[0]*nop-3,datamesh[0]*nop-2,datamesh[0]*nop-1])
            
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = keh8.flatten('F')

    KG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return KG

# LINEAR HEXAHEDRAL ELEMENT - isoparametric
def solid_hx8_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp):
    
    dofe = 8*datamesh[0]
    nel = datamesh[2]
    ith = np.zeros((nel*(dofe*dofe)))
    jth = np.zeros((nel*(dofe*dofe)))
    val = np.zeros((nel*(dofe*dofe)))
    
    # A = np.zeros((datamesh[2],1))
    for ee in range(datamesh[2]): 
       
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        nok = int(inci[ee,6])
        nol = int(inci[ee,7])
        nom = int(inci[ee,8])
        non = int(inci[ee,9])
        noo = int(inci[ee,10])
        nop = int(inci[ee,11])
        
        xi = coord[noi-1,1]
        yi = coord[noi-1,2]
        zi = coord[noi-1,3]        
        xj = coord[noj-1,1]
        yj = coord[noj-1,2]
        zj = coord[noj-1,3]        
        xk = coord[nok-1,1]
        yk = coord[nok-1,2]
        zk = coord[nok-1,3]        
        xl = coord[nol-1,1]
        yl = coord[nol-1,2]
        zl = coord[nol-1,3]        
        xm = coord[nom-1,1]
        ym = coord[nom-1,2]
        zm = coord[nom-1,3]        
        xn = coord[non-1,1]
        yn = coord[non-1,2]
        zn = coord[non-1,3]        
        xo = coord[noo-1,1]
        yo = coord[noo-1,2]
        zo = coord[noo-1,3]        
        xp = coord[nop-1,1]
        yp = coord[nop-1,2]
        zp = coord[nop-1,3] 
        
        matXY = np.array([[xi,yi,zi],[xj,yj,zj],[xk,yk,zk],[xl,yl,zl],\
                          [xm,ym,zm],[xn,yn,zn],[xo,yo,zo],[xp,yp,zp]])
                
        D = tabmat[int(inci[ee,2]-1),6]
        
        xp,wp = pointGauss(npp)
        x = xp[0,:]
        y = xp[1,:]
        z = xp[2,:]
        
        meh8 = np.zeros((dofe,dofe))
        for pp in range(0,npp):
            
            N1 = (1-x[pp])*(1-y[pp])*(1-z[pp])
            N2 = (1+x[pp])*(1-y[pp])*(1-z[pp])
            N3 = (1+x[pp])*(1+y[pp])*(1-z[pp])
            N4 = (1-x[pp])*(1+y[pp])*(1-z[pp])
            N5 = (1-x[pp])*(1-y[pp])*(1+z[pp])
            N6 = (1+x[pp])*(1-y[pp])*(1+z[pp])
            N7 = (1+x[pp])*(1+y[pp])*(1+z[pp])
            N8 = (1-x[pp])*(1+y[pp])*(1+z[pp]) 
    
            
            N = (1/8)*np.array([[N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0,0],\
                                 [0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0],\
                                 [0,0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8]])
            
            N1x = -(1-y[pp])*(1-z[pp])
            N2x =  (1-y[pp])*(1-z[pp])
            N3x =  (1+y[pp])*(1-z[pp])
            N4x = -(1+y[pp])*(1-z[pp])
            N5x = -(1-y[pp])*(1+z[pp])
            N6x =  (1-y[pp])*(1+z[pp])
            N7x =  (1+y[pp])*(1+z[pp])
            N8x = -(1+y[pp])*(1+z[pp]) 
            N1y = -(1-x[pp])*(1-z[pp])
            N2y = -(1+x[pp])*(1-z[pp])
            N3y =  (1+x[pp])*(1-z[pp])
            N4y =  (1-x[pp])*(1-z[pp])
            N5y = -(1-x[pp])*(1+z[pp])
            N6y = -(1+x[pp])*(1+z[pp])
            N7y =  (1+x[pp])*(1+z[pp])
            N8y =  (1-x[pp])*(1+z[pp])
            N1z = -(1-x[pp])*(1-y[pp])
            N2z = -(1+x[pp])*(1-y[pp])
            N3z = -(1+x[pp])*(1+y[pp])
            N4z = -(1-x[pp])*(1+y[pp])
            N5z =  (1-x[pp])*(1-y[pp])
            N6z =  (1+x[pp])*(1-y[pp])
            N7z =  (1+x[pp])*(1+y[pp])
            N8z =  (1-x[pp])*(1+y[pp])
            
            dN = (1/8)*np.array([[N1x,N2x,N3x,N4x,N5x,N6x,N7x,N8x],\
                                 [N1y,N2y,N3y,N4y,N5y,N6y,N7y,N8y],\
                                 [N1z,N2z,N3z,N4z,N5z,N6z,N7z,N8z]])
                
            J = np.dot(dN,matXY)
            
            detJ = np.linalg.det(J)
               
            meh8 += np.dot(np.transpose(N),N)*D*detJ*wp[pp]*wp[pp]
            
        loc = np.array([datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1,\
                        datamesh[0]*nok-3,datamesh[0]*nok-2,datamesh[0]*nok-1,datamesh[0]*nol-3,datamesh[0]*nol-2,datamesh[0]*nol-1,\
                        datamesh[0]*nom-3,datamesh[0]*nom-2,datamesh[0]*nom-1,datamesh[0]*non-3,datamesh[0]*non-2,datamesh[0]*non-1,\
                        datamesh[0]*noo-3,datamesh[0]*noo-2,datamesh[0]*noo-1,datamesh[0]*nop-3,datamesh[0]*nop-2,datamesh[0]*nop-1])
            
        loc_T = loc.reshape(1,dofe).T
        Ycopy_loc = np.tile(loc_T,(1,dofe))
        Ycopy_loc_T = np.transpose(Ycopy_loc)
        ith[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc.flatten('F')
        jth[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = Ycopy_loc_T.flatten('F')
        val[(dofe*dofe)*ee:(dofe*dofe)*(ee+1)] = meh8.flatten('F')

    MG = sp.csc_matrix((val, (ith, jth)), shape=(datamesh[3], datamesh[3]))
    return MG


#%% POST-PROCESS
# malha deformada solido 3D
def solid_hx8_deform(datamesh,coord,U,scale_mesh):
    Udef = np.zeros((datamesh[1],3),dtype=float)
    Umag = np.zeros((datamesh[1],1),dtype=float)
    for nn in range(1,datamesh[1]+1):
        Udef[nn-1,0] = U[datamesh[0]*nn-3,0]
        Udef[nn-1,1] = U[datamesh[0]*nn-2,0]
        Udef[nn-1,2] = U[datamesh[0]*nn-1,0]
        Umag[nn-1,0] = np.sqrt(U[datamesh[0]*nn-3,0]**2 + U[datamesh[0]*nn-2,0]**2 + U[datamesh[0]*nn-1,0]**2)
    
    meshDefU = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh*Udef)),axis=1)
    return  meshDefU,Udef,Umag

# tensao von-mises no elemento Solid
def solid_hx8_streqv(U,datamesh,inci,coord,typeMechMat,tabmat,npp):  
    strs_elm_eqv = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_xx = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_yy = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_zz = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_xy = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_yz = np.zeros((datamesh[2],1),dtype=float)
    strs_elm_zx = np.zeros((datamesh[2],1),dtype=float)
    
    strn_elm_eqv = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_xx = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_yy = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_zz = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_xy = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_yz = np.zeros((datamesh[2],1),dtype=float)
    strn_elm_zx = np.zeros((datamesh[2],1),dtype=float)
    
    for ee in range(datamesh[2]): 
        loading_bar_v1(100*(ee/datamesh[2]))
        noi = int(inci[ee,4])
        noj = int(inci[ee,5])
        nok = int(inci[ee,6])
        nol = int(inci[ee,7])
        nom = int(inci[ee,8])
        non = int(inci[ee,9])
        noo = int(inci[ee,10])
        nop = int(inci[ee,11])
        
        xi = coord[noi-1,1]
        yi = coord[noi-1,2]
        zi = coord[noi-1,3]        
        xj = coord[noj-1,1]
        yj = coord[noj-1,2]
        zj = coord[noj-1,3]        
        xk = coord[nok-1,1]
        yk = coord[nok-1,2]
        zk = coord[nok-1,3]        
        xl = coord[nol-1,1]
        yl = coord[nol-1,2]
        zl = coord[nol-1,3]        
        xm = coord[nom-1,1]
        ym = coord[nom-1,2]
        zm = coord[nom-1,3]        
        xn = coord[non-1,1]
        yn = coord[non-1,2]
        zn = coord[non-1,3]        
        xo = coord[noo-1,1]
        yo = coord[noo-1,2]
        zo = coord[noo-1,3]        
        xp = coord[nop-1,1]
        yp = coord[nop-1,2]
        zp = coord[nop-1,3]        
        
        matXY = np.array([[xi,yi,zi],[xj,yj,zj],[xk,yk,zk],[xl,yl,zl],\
                          [xm,ym,zm],[xn,yn,zn],[xo,yo,zo],[xp,yp,zp]])
                
        D = mechanicalMaterial(typeMechMat,tabmat,inci,ee)
        
        xp,wp = pointGauss(npp)
        x = xp[0,:]
        y = xp[1,:]
        z = xp[2,:]
        
        B = np.zeros((6,24))
        for pp in range(0,npp):
            N1x = -(1-y[pp])*(1-z[pp])
            N2x =  (1-y[pp])*(1-z[pp])
            N3x =  (1+y[pp])*(1-z[pp])
            N4x = -(1+y[pp])*(1-z[pp])
            N5x = -(1-y[pp])*(1+z[pp])
            N6x =  (1-y[pp])*(1+z[pp])
            N7x =  (1+y[pp])*(1+z[pp])
            N8x = -(1+y[pp])*(1+z[pp]) 
            N1y = -(1-x[pp])*(1-z[pp])
            N2y = -(1+x[pp])*(1-z[pp])
            N3y =  (1+x[pp])*(1-z[pp])
            N4y =  (1-x[pp])*(1-z[pp])
            N5y = -(1-x[pp])*(1+z[pp])
            N6y = -(1+x[pp])*(1+z[pp])
            N7y =  (1+x[pp])*(1+z[pp])
            N8y =  (1-x[pp])*(1+z[pp])
            N1z = -(1-x[pp])*(1-y[pp])
            N2z = -(1+x[pp])*(1-y[pp])
            N3z = -(1+x[pp])*(1+y[pp])
            N4z = -(1-x[pp])*(1+y[pp])
            N5z =  (1-x[pp])*(1-y[pp])
            N6z =  (1+x[pp])*(1-y[pp])
            N7z =  (1+x[pp])*(1+y[pp])
            N8z =  (1-x[pp])*(1+y[pp])
            
            dN = (1/8)*np.array([[N1x,N2x,N3x,N4x,N5x,N6x,N7x,N8x],\
                                 [N1y,N2y,N3y,N4y,N5y,N6y,N7y,N8y],\
                                 [N1z,N2z,N3z,N4z,N5z,N6z,N7z,N8z]])
            
            J = np.dot(dN,matXY)
                      
            invJ = np.linalg.inv(J)
            
            dN1 = np.array([[N1x],[N1y],[N1z]])
            dN2 = np.array([[N2x],[N2y],[N2z]])
            dN3 = np.array([[N3x],[N3y],[N3z]])
            dN4 = np.array([[N4x],[N4y],[N4z]])
            dN5 = np.array([[N5x],[N5y],[N5z]])
            dN6 = np.array([[N6x],[N6y],[N6z]])
            dN7 = np.array([[N7x],[N7y],[N7z]])
            dN8 = np.array([[N8x],[N8y],[N8z]])
            
            pN1 = np.dot(invJ,dN1)
            pN2 = np.dot(invJ,dN2)
            pN3 = np.dot(invJ,dN3)
            pN4 = np.dot(invJ,dN4)
            pN5 = np.dot(invJ,dN5)
            pN6 = np.dot(invJ,dN6)
            pN7 = np.dot(invJ,dN7)
            pN8 = np.dot(invJ,dN8)
            
            B1 = np.array([[pN1[0,0],0,0],\
                          [0,pN1[1,0],0],\
                          [0,0,pN1[2,0]],\
                          [pN1[1,0],pN1[0,0],0],\
                          [0,pN1[2,0],pN1[1,0]],\
                          [pN1[2,0],0,pN1[0,0]]])
           

            B2 = np.array([[pN2[0,0],0,0],\
                          [0,pN2[1,0],0],\
                          [0,0,pN2[2,0]],\
                          [pN2[1,0],pN2[0,0],0],\
                          [0,pN2[2,0],pN2[1,0]],\
                          [pN2[2,0],0,pN2[0,0]]])
                
            B3 = np.array([[pN3[0,0],0,0],\
                          [0,pN3[1,0],0],\
                          [0,0,pN3[2,0]],\
                          [pN3[1,0],pN3[0,0],0],\
                          [0,pN3[2,0],pN3[1,0]],\
                          [pN3[2,0],0,pN3[0,0]]])
                
            B4 = np.array([[pN4[0,0],0,0],\
                          [0,pN4[1,0],0],\
                          [0,0,pN4[2,0]],\
                          [pN4[1,0],pN4[0,0],0],\
                          [0,pN4[2,0],pN4[1,0]],\
                          [pN4[2,0],0,pN4[0,0]]])

            B5 = np.array([[pN5[0,0],0,0],\
                          [0,pN5[1,0],0],\
                          [0,0,pN5[2,0]],\
                          [pN5[1,0],pN5[0,0],0],\
                          [0,pN5[2,0],pN5[1,0]],\
                          [pN5[2,0],0,pN5[0,0]]])
           

            B6 = np.array([[pN6[0,0],0,0],\
                          [0,pN6[1,0],0],\
                          [0,0,pN6[2,0]],\
                          [pN6[1,0],pN6[0,0],0],\
                          [0,pN6[2,0],pN6[1,0]],\
                          [pN6[2,0],0,pN6[0,0]]])
                
            B7 = np.array([[pN7[0,0],0,0],\
                          [0,pN7[1,0],0],\
                          [0,0,pN7[2,0]],\
                          [pN7[1,0],pN7[0,0],0],\
                          [0,pN7[2,0],pN7[1,0]],\
                          [pN7[2,0],0,pN7[0,0]]])
                
            B8 = np.array([[pN8[0,0],0,0],\
                          [0,pN8[1,0],0],\
                          [0,0,pN8[2,0]],\
                          [pN8[1,0],pN8[0,0],0],\
                          [0,pN8[2,0],pN8[1,0]],\
                          [pN8[2,0],0,pN8[0,0]]])
            
            
            B += (1/8)*np.concatenate((B1,B2,B3,B4,B5,B6,B7,B8),axis=1)
             
        loc = np.array([datamesh[0]*noi-3,datamesh[0]*noi-2,datamesh[0]*noi-1,datamesh[0]*noj-3,datamesh[0]*noj-2,datamesh[0]*noj-1,\
                        datamesh[0]*nok-3,datamesh[0]*nok-2,datamesh[0]*nok-1,datamesh[0]*nol-3,datamesh[0]*nol-2,datamesh[0]*nol-1,\
                        datamesh[0]*nom-3,datamesh[0]*nom-2,datamesh[0]*nom-1,datamesh[0]*non-3,datamesh[0]*non-2,datamesh[0]*non-1,\
                        datamesh[0]*noo-3,datamesh[0]*noo-2,datamesh[0]*noo-1,datamesh[0]*nop-3,datamesh[0]*nop-2,datamesh[0]*nop-1])
      
        
        strain =  B@U[loc,0]
        sigma =  np.dot(D,strain)
        
        strn_elm_xx[ee] = strain[0]
        strn_elm_yy[ee] = strain[1]
        strn_elm_zz[ee] = strain[2]
        strn_elm_xy[ee] = strain[3]
        strn_elm_yz[ee] = strain[4]
        strn_elm_zx[ee] = strain[5]
        
        
        strs_elm_xx[ee] = sigma[0]
        strs_elm_yy[ee] = sigma[1]
        strs_elm_zz[ee] = sigma[2]
        strs_elm_xy[ee] = sigma[3]
        strs_elm_yz[ee] = sigma[4]
        strs_elm_zx[ee] = sigma[5]
        
        strs_elm_eqv[ee] = np.sqrt(0.5*((sigma[0]-sigma[1])**2 + (sigma[1]-sigma[2])**2 + (sigma[2]-sigma[0])**2 + 6*(sigma[3]**2 + sigma[4]**2 + sigma[5]**2))) 
        strn_elm_eqv[ee] = np.sqrt(0.5*((strain[0]-strain[1])**2 + (strain[1]-strain[2])**2 + (strain[2]-strain[0])**2 + 6*(strain[3]**2 + strain[4]**2 + strain[5]**2))) 
        
        stress_list = np.concatenate((strs_elm_eqv,strs_elm_xx,strs_elm_yy,strs_elm_zz,strs_elm_xy,strs_elm_yz,strs_elm_zx),axis=1)
        strain_list = np.concatenate((strn_elm_eqv,strn_elm_xx,strn_elm_yy,strn_elm_zz,strn_elm_xy,strn_elm_yz,strn_elm_zx),axis=1)
        
    return stress_list, strain_list