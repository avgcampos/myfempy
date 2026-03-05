'''

Simulação de vibração modal

'''
import sys
# setting path
sys.path.append('../myfempy')
import sys
import numpy as np
import matplotlib.pyplot as plt

from myfempy import newAnalysis
from myfempy import DynamicEigenLinear

import time

#===============================================================================
#                                   FEA
#===============================================================================
fea = newAnalysis(DynamicEigenLinear)

mat1 = {
    "NAME": "mat1",
    "VXY": 0.3,
    "EXX": 210E6,
    "RHO": 7850,
    }

mat2 = {
    "NAME": "mat2",
    "VXX": 0.3,
    "EXX": 1E6,
    "RHO": 1,
    }

geo = {
    "NAME": "espessura1",
    "THICKN": 0.1,
    # "DIM": [b, h, t, d],
    }

# MODEL SET

# nelx = 10
# nely = 5

# nodes = [[1, 0, 0, 0],
#          [2, nelx, 0, 0],
#          [3, nelx, nely, 0],
#          [4, 0, nely, 0],
#         #  [5, 1, 1, 0],
#         #  [6, 0, 1, 0],
#          ]

# conec = [[1, 1, 1, 1, 2, 3, 4],
#         #  [2, 2, 1, 2, 3, 4, 5],
#          ]

# MODEL SET
LX = 8 # 80
LY = 8 # 60

nelx = 4 # 40 # 80 # 160
nely = 4 # 30 # 60 # 120

# esize = 0.25

# # gmsh config
# points = [
#     [0, 0, 0],
#     [LX, 0, 0],
#     [LX, LY, 0],
#     [0, LY, 0]
# ]
         
# lines = [[1, 2],
#          [2, 3],
#          [3, 4],
#          [4, 1],
#          ]

# plane = [[1, 2, 3, 4]]

modeldata = {
    # "MESH": {
    #     'TYPE': 'add',
    #     'COORD': nodes,
    #     'INCI': conec,
    # },
        
    "MESH": {
        'TYPE': 'legacy',
        'LX': LX,
        'LY': LY,
        'NX': nelx,
        'NY': nely,
        },
        
    # "MESH": {
    #         'TYPE': 'gmsh',
    #         'filename': 'test_mfoop_gmsh_01',
    #         # "meshimport": 'object_dir',
    #         'pointlist': points,
    #         'linelist': lines,
    #         'planelist': plane,
    #         # 'arc': arcs,
    #         'meshconfig': {
    #             'mesh': 'quad4',   #quad4
    #             'sizeelement': 2 * esize,
    #             # 'extrude': 20,
    #             'meshmap': {'on': True,
    #                         'edge': 'all', #'all'
    #                     #  "numbernodes": 10,
    #             }}
    #     },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat2],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

forcespring_left = {
    'TYPE': 'forcenode',
    'DOF': 'spring2ground',
    'DIR': 'point',
    'TAG': 1,
    'VAL': [1000.0],
    }

forcespring_right = {
    'TYPE': 'forcenode',
    'DOF': 'spring2ground',
    'DIR': 'point',
    'TAG': 2,
    'VAL': [1000.0],
    }

massadd = {
    'TYPE': 'forcenode',
    'DOF': 'masspoint',
    'DIR': 'node', #'point', #'node',
    'LOC': {'x': LX/2, 'y': LY/2, 'z': 0},
    'VAL': [1.4],
    }

bc_node_left = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'node',                 # node
    'LOC': {'x': 0, 'y': LY/2, 'z': 0},
    }

bc_node_right = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'node',
    'LOC': {'x': LX, 'y': LY/2, 'z': 0},
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
                "LOAD": [],
                "BOUNDCOND": [],
    },
    }
fea.Physic(physicdata)

# loadaply = fea.getLoadApply()
# # bcaply = fea.getBCApply()
# print(loadaply)
# print(bcaply)
# print(fea.getLoadArray(loadaply))
# print(fea.getDirichletNH(bcaply))
# print('forca total [N]', np.sum(loadaply[:,2]))

previewset = {'RENDER': {'filename': 'squad', 'show': True, 'scale': 10, 'savepng': True, 'lines': False,
                        # 'plottags': {
                        # # 'line': True
                        # 'point': True
                        #       }
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# sys.exit()

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'mode',  # mode, freq, time ...
                        'start': 0,
                        'end': 10,
                        'step': 1},
             'SYMM':False,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

print(solverdata['solution']['FREQ'])

postprocdata = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    'structural': {'modes': True},
                    # 'structural': {'frf': True}
                    },
                "PLOTSET": {'filename': 'test_dynamic', 'savepng': True},
                # "TRACKER": {'frf': {'x': 16, 'y': 1, 'z': 0, 'dof':0}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }

postporc_result = fea.PostProcess(postprocdata)