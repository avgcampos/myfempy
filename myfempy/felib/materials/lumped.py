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
    '''SpringLinear material'''
    
    def __init__(self, tabmat, inci, num_elm):
        
        self.S = tabmat[int(inci[num_elm,2])-1,7] # spring stiffness
        self.D = tabmat[int(inci[num_elm,2])-1,8] # spring dampe
        
    def springlinear(self):
        
        S = self.S 
        D = self.D 
        E = [S, D]
                
        return E
        
# class SpringNonlin: