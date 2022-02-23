# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
_______________________________________________________________________________
 ~~~~~~~~~~              COMPORTAMENTO DO MATERIAL                   ~~~~~~~~~~

> FUNCIONALIDADES
--- ENTRADAS: 
--- SAIDA: 

===============================================================================

> ATUALIZACOES DA VERSAO:

_______________________________________________________________________________
"""

import numpy as np

#%% Material
def mech_mate(datamesh,tabmat,inci,num_elm):
    if datamesh['matdefi'] == 'lumped':
        if datamesh['matbeha'] == 'springlinear':
            E = tabmat[int(inci[num_elm,2])-1,0]  
            C = tabmat[int(inci[num_elm,2])-1,1]  
            D = [E, C]
            
        elif datamesh['matbeha'] == 'springnonlin':
            print('no implement yet')
   
    if datamesh['matdefi'] == 'isotropic':
        if datamesh['matbeha'] == 'uniaxial':
            E = tabmat[int(inci[num_elm,2])-1,0]  
            G = tabmat[int(inci[num_elm,2])-1,2]  
            D = [E, G]
            
        elif datamesh['matbeha'] == 'planestress':   
            E = tabmat[int(inci[num_elm,2])-1,0]  
            v = tabmat[int(inci[num_elm,2])-1,1]  
            D = E/(1-v**2)*(np.array([[1,v,0],\
                                      [v,1,0],\
                                      [0,0,(1-v)/2]]))
            
        elif datamesh['matbeha'] == 'planestrain':  
            E = tabmat[int(inci[num_elm,2])-1,0]  
            v = tabmat[int(inci[num_elm,2])-1,1]
            D = E/((1+v)*(1-2*v))*(np.array([[1-v,v,0],\
                                             [v,1-v,0],\
                                             [0,0,(1-2*v)/2]]))
        
        elif datamesh['matbeha'] == 'threeaxial':
            E = tabmat[int(inci[num_elm,2])-1,0]  
            v = tabmat[int(inci[num_elm,2])-1,1]
            D = E/((1+v)*(1-2*v))*np.array([[1-v,v,v,0,0,0],\
                                            [v,1-v,v,0,0,0],\
                                            [v,v,1-v,0,0,0],\
                                            [0,0,0,(1-2*v)/2,0,0],\
                                            [0,0,0,0,(1-2*v)/2,0],\
                                            [0,0,0,0,0,(1-2*v)/2]])
                    
        else:
             print("input erro: mat_behavior don't defined")
    
    elif datamesh['matdefi'] == 'orthotropic':
        if datamesh['matbeha'] == 'uniaxial':
            D = tabmat[int(inci[num_elm,2])-1,0]
        elif datamesh['matbeha'] == 'planestress':   
            E = tabmat[int(inci[num_elm,2])-1,0]  
            v = tabmat[int(inci[num_elm,2])-1,1]
            D = E/(1-v**2)*(np.array([[1,v,0],\
                                      [v,1,0],\
                                      [0,0,(1-v)/2]]))
            
        elif datamesh['matbeha'] == 'planestrain':  
            E = tabmat[int(inci[num_elm,2])-1,0]  
            v = tabmat[int(inci[num_elm,2])-1,1]
            D = E/((1+v)*(1-2*v))*(np.array([[1-v,v,0],\
                                             [v,1-v,0],\
                                             [0,0,(1-2*v)/2]]))
        else:
             print("input erro: mat_behavior don't defined")    
    
    return D
