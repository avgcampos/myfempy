from colorama import Fore, Back, Style

def myfempy_logo():
    import os
    if os.name == 'posix': # mac/linux
        _ = os.system('clear')
        print("                                                                               ")
        print("                              ___   ___   __  __                               ")
        print("               _ __    _  _  | __| | __| |  \/  |  _ __   _  _                 ")
        print("              | '  \  | || | | _|  | _|  | |\/| | | '_ \ | || |                ")
        print("              |_|_|_|  \_, | |_|   |___| |_|  |_| | .__/  \_, |                ")
        print("                       |__/                       |_|     |__/                 ")
        print('~~~~~~      Mechanical studY with Finite Element Method in PYthon       ~~~~~~~')
        print('~~~~~~                  copyright all rights reserved                   ~~~~~~~')
        print('~~~~~~                       VERSION BETA 0.0.2                         ~~~~~~~')
        print("                                                                               ")
    else: # windows
        _ = os.system('cls')
        print("                                                                               ")
        print(Fore.RED + Style.BRIGHT+"                              ___   ___   __  __                               ")
        print(Fore.RED + Style.BRIGHT+"               _ __    _  _  | __| | __| |  \/  |  _ __   _  _                 ")
        print(Fore.YELLOW + Style.BRIGHT+"              | '  \  | || | | _|  | _|  | |\/| | | '_ \ | || |                ")
        print(Fore.GREEN + Style.BRIGHT+"              |_|_|_|  \_, | |_|   |___| |_|  |_| | .__/  \_, |                ")
        print(Fore.CYAN + Style.BRIGHT+"                       |__/                       |_|     |__/                 ")
        print(Fore.WHITE + Style.BRIGHT+'~~~~~~      Mechanical studY with Finite Element Method in PYthon       ~~~~~~~')
        print(Fore.WHITE + Style.BRIGHT+'~~~~~~                  copyright all rights reserved                   ~~~~~~~')
        print(Fore.WHITE + Style.BRIGHT+'~~~~~~                       VERSION BETA 0.0.2                         ~~~~~~~')
        print(Fore.WHITE + Style.BRIGHT+"                                                                               ")

#%% TEXT MYFEMPY ABOUT
def myfempy_about(welcome_true):        
    if welcome_true == 'true':
        f = open('myfempy_about.txt','r',encoding="utf8")
        file_contents = f.read()
        print(file_contents)
        f.close()
        input('PRESS ENTER TO CONTINUE ...')
    else:
        print("===============================================================================")
        
def myfempy_version():
    version_current = 'vBeta-1'
    myfempy_main_name = 'MYFEMPY'
    myfempy_version = myfempy_main_name+' - '+version_current
    return myfempy_version