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

# import sys
import os

# sys.path.insert(0, os.getcwd())

def create_user_path(file_dir):
    
    path = os.getcwd()
    
    if not os.path.exists(path+file_dir):
        # print('nova pasta criada no diretório "myfempy/user/"')
        os.makedirs(path+file_dir)
    
    # print('esta pasta já existe no diretório  "myfempy/user/"')
    path_user  = str(path+file_dir)
    return path_user 

