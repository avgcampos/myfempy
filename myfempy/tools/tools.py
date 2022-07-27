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
# sys.path.insert(0, os.getcwd() )


def get_version():
    try:
        VERSIONFILE = (os.getcwd()+"/myfempy/bin/version.py.py")
        verstrline = open(VERSIONFILE,"rt").read()
        version = verstrline.split("=")[1].replace("\n", "").replace("'", "")
    except:
        version = "unknown"
        
    return version
    


def get_logo():
    f = open((os.getcwd()+'/myfempy/bin/logo.txt'),'r', encoding="utf8")
    file_contents = f.read()
    print(file_contents)
    print('myfempy version:',get_version())
    f.close()
    
    
def clear_console():
    if os.name == 'posix': # linux/mac
        _ = os.system('clear')
        
    else: # windows
        _ = os.system('cls')
     

def gen_user_dir(path, file_dir):
    
    # sys.path.insert(0, os.getcwd() )
    
    if not os.path.exists(path+file_dir):
        # print('nova pasta criada no diretório "myfempy/user/"')
        os.makedirs(path+file_dir)
    
    # print('esta pasta já existe no diretório  "myfempy/user/"')
    path_user  = str(path+file_dir)
    return path_user 



def loading_bar_v1(pct, name):
    # sys.stdout.write(Fore.CYAN + Style.BRIGHT+"\r|%-50s"%('#'*int(pct*0.5)+'|'+str(round(pct))+'%'))
    sys.stdout.write("\r"+name+": "+"|%-50s"%('%'*int(pct*0.5)+'|'+str(round(pct))+'%'))
    sys.stdout.flush()

      
 
def print_console(sc):
    
    if sc == 'pre':
        
        clear_console()
        
        print('=================================================================================')
        
        get_logo()
        
        print('')
        print('\r******************   P R E - P R O C E S S   L O A D I N G   ******************')
    
    elif sc == 'mesh':
        print('\r**********************   G E N E R A T I N G   M E S H   **********************')
        
    elif sc == 'solver':
        print('\r********************   S O L V I N G   E Q U A T I O N S   ********************')
                                                                      
        
    elif sc == 'post':
        print('\r***************   P O S T - P R O C E S S   C O M P U T I N G   ***************')
        
    elif sc == 'thank':
        print('\n')
        print('***************       A N A L Y S I S   S U C C E S S F U L       ***************')
        print('***************             M Y F E M P Y @ 2 0 2 2               ***************')
        print('***************        T H A N K   Y O U   F O R   U S E !        ***************')
        print('=================================================================================')
    
    