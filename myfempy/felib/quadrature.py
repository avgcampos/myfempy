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

sqrt3 = np.sqrt(3)

class Quadrature:
        
    def no_interpol(npp):
        if npp == 1:
           xp = np.array([[1.0], [1.0], [1.0]])
           wp = np.array([1])
           return xp, wp
           
    #-----------------------------------------------------------------------------#
    #%% Pontos de integração Gauss
    def gaussian(npp):

        if npp == 1:
            xp = np.array([[1/3],[1/3]])
            wp = np.array([1])
            return xp, wp
        
        elif npp == 2:
            xp = np.array([[-0.577350269],[0.577350269]])
            wp = np.array([1])
            return xp, wp
        
        elif npp == 3:
           xp = np.array([[2/3, 1/6, 1/6],\
                          [1/6, 1/6, 2/3]])
           wp = np.array([1/3,1/3,1/3])
           return xp, wp
        
        elif npp == 4:
           xp = np.array([[-1/sqrt3 , 1/sqrt3 , 1/sqrt3 , -1/sqrt3 ],\
                          [-1/sqrt3 , -1/sqrt3 , 1/sqrt3 , 1/sqrt3 ]])
           wp = np.array([1,1,1,1])
           return xp, wp
          
        elif npp==8:
            xp = np.array([[-1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3), -1/np.sqrt(3)],\
                          [-1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)],
                          [-1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]])
            wp = np.array([1,1,1,1,1,1,1,1])
            return xp, wp
               
        else:
            print('erro pontos de integração Gauss')
        
        
