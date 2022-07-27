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

# Material
#-----------------------------------------------------------------------------#
class Elasticity:
    '''Isotropic material'''
    
    def __init__(self, tabmat, inci, num_elm):
        
        self.E = tabmat[int(inci[num_elm,2])-1,0] # material elasticity
        self.v = tabmat[int(inci[num_elm,2])-1,1] # material poisson ratio    
    
    def isotropic(self):
                       
        D = np.zeros((3,3))
        D[0,0] = self.E*(1.0-self.v)/((1+self.v)*(1.0-2.0*self.v))
        D[0,1] = D[0,0]*self.v/(1.0-self.v)
        D[1,0] = D[0,1]
        D[1,1] = D[0,0]
        D[2,2] = D[0,0]*0.5*(1.0-2.0*self.v)/(1.0-self.v)
        
        
        return D
               
               
# class Ortotropic: