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


import sys
import os
import imp
import numpy as np
# from colorama import Fore, Back, Style


#-----------------------------------------------------------------------------#
def gmsh_elm_type(elemid):
    # 'vtkkey' >>' Type "elem" in myfempy':CellType --> https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    celltype = {
        '110':1,
        '120':1,
        '130':1,
        '131':8,
        '140':1,
        '141':1,
        '142':8,
        '210':2,
        '211':9,
        '220':3,
        '221':16,
        '230':16,
        '231':16,
        '240':16,
        '310':4,
        '320':5,
        '321':17,
        }
    
    gmshtyp = celltype[elemid]
    return gmshtyp


def convert_from_msh1(filename):
    
    file_imp = (filename+'.msh1')
    with open(file_imp,'r') as file_object:
        # line = file_object.readline()
        # lineaux = line.split()
        
        file_object.readline()
        NNOD = int(file_object.readline())
        nodelist = [[None]*4]
        for ii in range(0,NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            nodelist.append([float(contstr[0]),float(contstr[1]),float(contstr[2]),float(contstr[3])])
        nodelist = nodelist[1::][::]
        
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = []
        for kk in range(0,NELM):
            line = file_object.readline()
            lineaux = line.split()
            conec_elm.append(list(map(float, lineaux[:])))     
        # conec_elm = conec_elm[1::][::]
            
    return conec_elm, nodelist