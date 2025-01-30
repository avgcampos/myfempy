#!/usr/bin/env python
__doc__ = """
Utils
"""
import os
import sys
from art import tprint
import importlib.metadata

def get_version():
    """_summary_

    Returns:
        _description_
    """
    try:
        __version__ = importlib.metadata.version("myfempy")
    except:
        __version__ = "dev"
    return __version__


def get_logo():
    """_summary_"""
    print("=================================================================================")
    tprint("             myfempy", font="ogre")
    print("version: ", get_version())
    print("=================================================================================")


def newDir(file_dir):
    """_summary_

    Arguments:
        file_dir -- _description_

    Returns:
        _description_
    """
    path = os.getcwd()
    if not os.path.exists(path + "/" + file_dir):
        os.makedirs(path + "/" + file_dir)
    path_user = str(path + "/" + file_dir)
    return path_user


def clear_console():
    """_summary_"""
    if os.name == "posix":  # linux/mac
        _ = os.system("clear")
    else:  # windows
        _ = os.system("cls")


def loading_bar_v1(pct, name):
    """_summary_

    Arguments:
        pct -- _description_
        name -- _description_
    """
    sys.stdout.write(
        "\r"
        # + name
        + ">> "
        + "|%-50s" % ("%" * int(pct * 0.7) + "|" + str(round(pct)) + "%")
    )
    sys.stdout.flush()


def print_console(sc):
    """_summary_

    Arguments:
        sc -- _description_
    """
    if sc == "pre":
        print(
            "\r******************   P R E - P R O C E S S   L O A D I N G   ******************"
        )
    elif sc == "mesh":
        print(
            "\r**********************   G E N E R A T I N G   M E S H   **********************"
        )
    elif sc == "solver":
        print(
            "\r********************   S O L V I N G   E Q U A T I O N S   ********************\n"
        )
    elif sc == "post":
        print(
            "\r***************   P O S T - P R O C E S S   C O M P U T I N G   ***************"
        )
    elif sc == "thank":
        print("\n")
        print(
            "***************       A N A L Y S I S   S U C C E S S F U L       ***************"
        )
        print(
            "***************                   M Y F E M P Y                   ***************"
        )
        print(
            "***************        T H A N K   Y O U   F O R   U S E !        ***************"
        )
