#!/usr/bin/env python
__doc__ = """
Path Setting
"""
import os

def newDir(file_dir):
    """_summary_

    Arguments:
        file_dir -- _description_

    Returns:
        _description_
    """
    path = os.getcwd()
    if not os.path.exists(path+'/'+file_dir):
        os.makedirs(path+'/'+file_dir)
    path_user = str(path+'/'+file_dir)
    return path_user



# def gen_user_dir(path, file_dir):
#     """_summary_

#     Arguments:
#         path -- _description_
#         file_dir -- _description_

#     Returns:
#         _description_
#     """
#     if not os.path.exists(path + file_dir):
#         os.makedirs(path + file_dir)
#     path_user = str(path + file_dir)
#     return path_user