
import sys
import numpy as np
import matplotlib.pyplot as plt

from myfempy import newAnalysis
from myfempy import DynamicHarmonicResponseLinear

import time

#===============================================================================
#                                   FEA
#===============================================================================

fea = newAnalysis(DynamicHarmonicResponseLinear)

mat1 = {
    "NAME": "mat1",
    "VXX": 0.25,
    "EXX": 250E9,
    "RHO": 7800
    }

mat2 = {
    "NAME": "mat2",
    "VXX": 0.25,
    "EXX": 1E6,
    "RHO": 1
    }

geo = {
    "NAME": "espessura1",
    "THICKN": 1.0,
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


# gmsh config
# points = [
#     [0, 0, 0],
#     [10, 0, 0],
#     [20, 0, 0],
#     [20, 10, 0],
#     [10, 10, 0],
#     [0, 10, 0],
# ]
         
# lines = [[1, 2],
#          [2, 5],
#          [5, 6],
#          [6, 1],
#          [2, 3],
#          [3, 4],
#          [4, 5],
#          ]

# plane = [[1, 2, 3, 4],
#          [5, 6, 7, 2]]

modeldata = {
    # "MESH": {
    #     'TYPE': 'add',
    #     'COORD': nodes,
    #     'INCI': conec,
    # },
        
    "MESH": {
        'TYPE': 'legacy',
        'LX': 16,
        'LY': 2,
        'NX': 4*16,
        'NY': 4*2,
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
    #             'sizeelement': 1,
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
        "PROPMAT": [mat1],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

s1 = {
    'TYPE': 'forcenode',
    'DOF': 'masspoint',
    'DIR': 'node',#'point', #'node',
    'LOC': {'x': 20, 'y': 5, 'z': 0},
    # 'TAG': 4,
    'VAL': [0],
    }

f2 = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 16, 'y': 1, 'z': 0},
    'VAL': [-1000.0],
    }

# bc2 = {
#     'TYPE': 'fixed',
#     'DOF': 'all',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }

bc1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }

bc2 = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'node',
    'LOC': {'x': 10000, 'y': 0, 'z': 0},
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [f2],
               "BOUNDCOND": [bc1],
    },
}

fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'michell_2ss', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# sys.exit()
# #-------------------------------- SOLVER -------------------------------------#
solverset = {
            # "SOLVER": 'SLD', #SLD
            # 'TOL': 1E-10,
            "STEPSET": {
                # 'type': 'table',  # mode, freq, time ...
                'start': 0,
                'end': 200,
                'step': 0.5},
             }
solverdata = fea.Solve(solverset)

# print(solverdata['solution']['U'].shape)
# print(solverdata['solution']['FREQ'])

postprocdata = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    # 'structural': {'modes': True},
                    'structural': {'frf': True}
                    },
                "PLOTSET": {'filename': 'test_dynamic', 'savepng': True},
                "TRACKER": {'frf': {'x': 16, 'y': 1, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }

postporc_result = fea.PostProcess(postprocdata)