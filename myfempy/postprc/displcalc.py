# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __                                
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/ 

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""

import numpy as np
# import scipy.sparse as sp
# from myfempy.felib.material import mech_mate
# from myfempy.felib.intpoint import gauss_point
from myfempy.bin.tools import loading_bar_v1

#%%------------------------------------------------------------------------------

class Deformation:
    
    def __init__(self, modelinfo):
        
        self.nodedof = modelinfo['nodedof'][0]
        self.nnode = len(modelinfo['coord'])
        self.coord = modelinfo['coord']
        
    def ux(self, U, scale):
    
        Udef = np.zeros((self.nnode,3),dtype=float)
        Umag = np.zeros((self.nnode,1),dtype=float)
        loading_bar_v1(0,'POST-PROCESSING')
        for nn in range(1,self.nnode+1):
            loading_bar_v1(100*((nn)/self.nnode),'POST-PROCESSING')
            Udef[nn-1,0] = U[self.nodedof*nn-2]
            Udef[nn-1,1] = U[self.nodedof*nn-1]
            Umag[nn-1,0] = np.sqrt(U[self.nodedof*nn-2]**2 + U[self.nodedof*nn-1]**2)
    
        meshDefU = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Udef)),axis=1)
        
        result = np.concatenate((Umag, Udef),axis=1) 
        title =  ['DISPL_X','DISPL_Y','DISPL_Z']
        
        return   meshDefU, result
    
    def uy_rz(self, U, scale):
    
        Udef = np.zeros((self.nnode,3),dtype=float)
        Umag = np.zeros((self.nnode,1),dtype=float)
        loading_bar_v1(0,'POST-PROCESSING')
        for nn in range(1,self.nnode+1):
            loading_bar_v1(100*((nn)/self.nnode),'POST-PROCESSING')
            Udef[nn-1,1] = U[self.nodedof*nn-2]
            # Rdef[nn-1,2] = U[self.nodedof*nn-1,0]
            Umag[nn-1,0] = Udef[nn-1,1]
        
        meshDefU = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Udef)),axis=1)
        # meshRotZ = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Rdef)),axis=1)
        
        result = np.concatenate((Umag, Udef),axis=1)
        title =  ['DISPL_X','DISPL_Y','DISPL_Z']
        
        return  meshDefU, result
      
    def ux_uy_rz(self, U, scale):
                
        Udef = np.zeros((self.nnode,3),dtype=float)
        Umag = np.zeros((self.nnode,1),dtype=float)
        loading_bar_v1(0,'POST-PROCESSING')
        for nn in range(1,self.nnode+1):
            loading_bar_v1(100*((nn)/self.nnode),'POST-PROCESSING')
            Udef[nn-1,0] = U[self.nodedof*nn-3]
            Udef[nn-1,1] = U[self.nodedof*nn-2]
            # Rdef[nn-1,2] = U[self.nodedof*nn-1,0]
            Umag[nn-1,0] = np.sqrt(U[self.nodedof*nn-3]**2 + U[self.nodedof*nn-2]**2)
        
        meshDefU = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Udef)),axis=1)
        # meshRotZ = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Rdef)),axis=1)
        
        result = np.concatenate((Umag, Udef),axis=1) 
        title =  ['DISPL_X','DISPL_Y','DISPL_Z']
        
        return  meshDefU, result
      
              
      # # malha deformada frame 3D
    def ux_uy_uz_rx_ry_rz(self, U, scale):
        
        Udef = np.zeros((self.nnode,3),dtype=float)
        Umag = np.zeros((self.nnode,1),dtype=float)
        for nn in range(1,self.nnode+1):
            Udef[nn-1,0] = U[self.nodedof*nn-6]
            Udef[nn-1,1] = U[self.nodedof*nn-5]
            Udef[nn-1,2] = U[self.nodedof*nn-4]
            # Rdef[nn-1,0] = U[self.nodedof*nn-3,0]
            # Rdef[nn-1,1] = U[self.nodedof*nn-2,0]
            # Rdef[nn-1,2] = U[self.nodedof*nn-1,0]
            Umag[nn-1,0] = np.sqrt(U[self.nodedof*nn-6]**2 + U[self.nodedof*nn-5]**2 + U[self.nodedof*nn-4]**2)
        
        meshDefU = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Udef)),axis=1)
        # meshRotZ = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Rdef)),axis=1)
        
        result = np.concatenate((Umag, Udef),axis=1) 
        title =  ['DISPL_X','DISPL_Y','DISPL_Z']
        
        return  meshDefU, result
    
   
    def ux_uy(self, U, scale):
        
        Udef = np.zeros((self.nnode,3),dtype=float)
        Umag = np.zeros((self.nnode,1),dtype=float)
        loading_bar_v1(0,'POST-PROCESSING')
        for nn in range(1,self.nnode+1):
            loading_bar_v1(100*((nn)/self.nnode),'POST-PROCESSING')
            Udef[nn-1,0] = U[self.nodedof*nn-2]
            Udef[nn-1,1] = U[self.nodedof*nn-1]
            Umag[nn-1,0] = np.sqrt(U[self.nodedof*nn-2]**2 + U[self.nodedof*nn-1]**2)
        
        meshDefU = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Udef)),axis=1)
        result = np.concatenate((Umag, Udef),axis=1)
        title =  ['DISPL_X','DISPL_Y','DISPL_Z']
        
        return  meshDefU, result
    
    
    def ux_uy_uz(self, U, scale):
    
        Udef = np.zeros((self.nnode,3),dtype=float)
        Umag = np.zeros((self.nnode,1),dtype=float)
        loading_bar_v1(0,'POST-PROCESSING')
        for nn in range(1,self.nnode+1):
            loading_bar_v1(100*((nn)/self.nnode),'POST-PROCESSING')
            Udef[nn-1,0] = U[self.nodedof*nn-3]
            Udef[nn-1,1] = U[self.nodedof*nn-2]
            Udef[nn-1,2] = U[self.nodedof*nn-1]
            Umag[nn-1,0] = np.sqrt(U[self.nodedof*nn-3]**2 + U[self.nodedof*nn-2]**2 + U[self.nodedof*nn-1]**2)
        
        meshDefU = np.concatenate((self.coord[:,0][:, np.newaxis],np.add(self.coord[:,1:],scale*Udef)),axis=1)
           
        result = np.concatenate((Umag, Udef),axis=1)    
        title =  ['DISPL_X','DISPL_Y','DISPL_Z']
        
        return  meshDefU, result