
__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""

import importlib.metadata
import os
import sys
from art import tprint

def get_about():
    try:
        with open(os.getcwd()+'/myfempy/utils/about.txt', 'r', encoding='utf-8') as file:
            conteudo = file.read()
            print(conteudo)
    except FileNotFoundError:
        print("About file not found")

def get_version():
    try:
        __version__ = importlib.metadata.version("myfempy")
    except:
        __version__ = "dev"
    return __version__


def get_logo():
    print(
        "================================================================================="
    )
    tprint("             myfempy", font="ogre")
    print("version: ", get_version())
    print(
        "================================================================================="
    )


def newDir(file_dir):
    path = os.getcwd()
    if not os.path.exists(path + "/" + file_dir):
        os.makedirs(path + "/" + file_dir)
    path_user = str(path + "/" + file_dir)
    return path_user


def clear_console():
    if os.name == "posix":  # linux/mac
        _ = os.system("clear")
    else:  # windows
        _ = os.system("cls")


def loading_bar_v1(pct, name):
    sys.stdout.write(
        "\r"
        + name
        + ">> "
        + "|%-50s" % ("%" * int(pct * 0.7) + "|" + str(round(pct)) + "%")
    )
    sys.stdout.flush()


def print_console(sc):
    if sc == "mesh":
        print(
            "\r[1 / 5]   G E N E R A T I N G   M O D E L"
        )
    elif sc == "phy":
        print(
            "\r[2 / 5]   P H Y S I C S ' S   L O A D I N G"
        )
    elif sc == "solver":
        print(
            "\r[3 / 5]   S O L V I N G   E Q U A T I O N S"
        )
    elif sc == 'succ':
        print(
            "\r[4 / 5]   A N A L Y S I S   S U C C E S S F U L"
        )
    elif sc == "post":
        print(
            "\r[5 / 5]   P O S T - P R O C E S S   C O M P U T I N G"
        )
    elif sc == "thank":
        # print(
        #     "\r***************                   M Y F E M P Y                   ***************"
        # )
        print(
            "\r***************        T H A N K   Y O U   F O R   U S E !        ***************"
        )
