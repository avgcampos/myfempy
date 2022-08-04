#!/usr/bin/env python
__doc__ ="""
Tools and Misc.
"""
import sys
import os
import art


def get_version():
    from myfempy import version
    return version.__version__


def get_logo():
    art.tprint('myfempy', font="ogre")
    print('myfempy version:', get_version())


def clear_console():
    if os.name == 'posix':  # linux/mac
        _ = os.system('clear')
    else:  # windows
        _ = os.system('cls')


def gen_user_dir(path, file_dir):
    if not os.path.exists(path+file_dir):
        os.makedirs(path+file_dir)
    path_user = str(path+file_dir)
    return path_user


def loading_bar_v1(pct, name):
    sys.stdout.write("\r"+name+": "+"|%-50s" %
                     ('%'*int(pct*0.5)+'|'+str(round(pct))+'%'))
    sys.stdout.flush()


def print_console(sc):
    if sc == 'pre':
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
