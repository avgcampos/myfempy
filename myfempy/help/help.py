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

from colorama import Fore, Style
import os
import version
import webbrowser as wb

def about():
    f = open('menu/about.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
     
    
def copyr():
    f = open('menu/copyr.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)

    
def comm():
    f = open('menu/comm.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
    

def features():   
    f = open('menu/feat.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
        

def release():
    f = open('menu/relea.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
    print(" ")
    print(Fore.WHITE + Style.BRIGHT+'> LAST VERSION')
    currversion()

def license():
    f = open('menu/license.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
    
def paths():
    f = open('menu/paths.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
    

def head():
    f = open('menu/head.txt','r',encoding="utf8")
    file_contents = f.read()
    print(Fore.WHITE + Style.BRIGHT+file_contents)
    f.close()
    

def help_webpage():
    wb.open_new_tab('file:///'+os.getcwd()+'/' + 'help.html')


def currversion():
    project_name = 'myfempy'
    version_launch = 'mar 2022'
    version_current = version._version
    print(Fore.WHITE + Style.BRIGHT+project_name+' CURRENT VERSION: '+' '+version_current)
    print(Fore.WHITE + Style.BRIGHT+project_name+' LAUNCH DATA: '+' '+version_launch)
    
    
def clear_console():
    if os.name == 'posix': # mac/linux
        _ = os.system('clear')
        
    else: # windows
        _ = os.system('cls')
    
    
#%%
def help_menu():
    print(" ")
    print("MYFEMPY HELP MENU")
    print("ENTER WITH THE FOLLOWING OPTIONS")
    print("   a) About")
    print("   f) Features")
    print("   v) Version Release")
    print("   c) Citation")
    print("   l) License")
    print("   p) Paths")
    print("   w) WebPage Help")
    print("   q) Quit")
    choice = input("Choice: ")    
    if choice.lower() in ['a','c', 'f','v', 'p', 'q', 'w', 'l']:
        return choice.lower()
    else:
        print(choice +"?")
        print("Invalid option")
        return None    
    
def menu_loop():
    while True:
        choice = help_menu()
        if choice == None:
            continue
        if choice == 'q':
            print( "Exiting...")
            break     # jump out of while loop
        elif choice == 'a':
            head()
            currversion()
            about()
        elif choice == 'c':
            copyr()
        elif choice == 'f':
            features()
        elif choice == 'v':
            release()
        elif choice == 'p':
            paths()
        elif choice == 'w':
            help_webpage()
        elif choice == 'l':
            license()
        else:
            print("Invalid choice.")
    

def mfpy_help():
    """
    MYFEMPY HELP MENU

    Returns
    -------
    None.

    """
    clear_console()
    head()
    currversion()
    menu_loop()
    
# The following makes this program start running at main_loop() 
# when executed as a stand-alone program.    
if __name__ == '__main__':
   mfpy_help()
    