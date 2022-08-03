#!/usr/bin/env python
"""
Path Setting
"""
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"

# import sys
import os


def create_user_path(file_dir):
    path = os.getcwd()
    if not os.path.exists(path+file_dir):
        os.makedirs(path+file_dir)
    path_user = str(path+file_dir)
    return path_user
