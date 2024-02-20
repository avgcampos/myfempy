#!/usr/bin/env python
__doc__ = """
Utils
"""
import os
import sys


def get_version():
    """_summary_

    Returns:
        _description_
    """
    from myfempy import __version__

    return __version__


def get_logo():
    """_summary_"""
    import art

    art.tprint("myfempy", font="ogre")
    print("myfempy version:", get_version())


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
        + name
        + ": "
        + "|%-50s" % ("%" * int(pct * 0.5) + "|" + str(round(pct)) + "%")
    )
    sys.stdout.flush()


def print_console(sc):
    """_summary_

    Arguments:
        sc -- _description_
    """
    if sc == "pre":
        print(
            "================================================================================="
        )
        get_logo()
        print("")
        print(
            "\r******************   P R E - P R O C E S S   L O A D I N G   ******************"
        )
    elif sc == "mesh":
        print(
            "\r**********************   G E N E R A T I N G   M E S H   **********************"
        )
    elif sc == "solver":
        print(
            "\r********************   S O L V I N G   E Q U A T I O N S   ********************"
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
            "***************             M Y F E M P Y @ 2 0 2 2               ***************"
        )
        print(
            "***************        T H A N K   Y O U   F O R   U S E !        ***************"
        )
        print(
            "================================================================================="
        )
