
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

# MODEL SET
LX = 16 # 80
LY = 2 # 60
nelx = 4 # 40 # 80 # 160
nely = 3 # 30 # 60 # 120

# gmsh config
points = [
    [0, 0, 0],
    [LX, 0, 0],
    [LX, LY, 0],
    [0, LY, 0]
]
         
lines = [[1, 2],
         [2, 3],
         [3, 4],
         [4, 1],
         ]

plane = [[1, 2, 3, 4]]

modeldata = {
    # "MESH": {
    #     'TYPE': 'add',
    #     'COORD': nodes,
    #     'INCI': conec,
    # },
        
    # "MESH": {
    #     'TYPE': 'legacy',
    #     'LX': LX,
    #     'LY': LY,
    #     'NX': nelx,
    #     'NY': nely,
    #     },
        
    "MESH": {
            'TYPE': 'gmsh',
            'filename': 'test_mfoop_gmsh_01',
            # "meshimport": 'object_dir',
            'pointlist': points,
            'linelist': lines,
            'planelist': plane,
            # 'arc': arcs,
            'meshconfig': {
                'mesh': 'quad4',   #quad4
                'sizeelement': 1,
                # 'extrude': 20,
                'meshmap': {'on': True,
                            'edge': 'all', #'all'
                        #  "numbernodes": 10,
                }}
        },

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

forcespring = {
    'TYPE': 'forcenode',
    'DOF': 'spring2ground',
    'DIR': 'line',
    'TAG': 3,
    'VAL': [100000.0],
    }

forcesmass = {
    'TYPE': 'forcenode',
    'DOF': 'masspoint',
    'DIR': 'point', #'point', #'node',
    'TAG': 3,
    'VAL': [1000000],
    }

forcesmass2 = {
    'TYPE': 'forcenode',
    'DOF': 'masspoint',
    'DIR': 'point', #'point', #'node',
    'TAG': 2,
    'VAL': [1000000],
    }

f2 = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'node',
    'LOC': {'x': 16, 'y': 1, 'z': 0},
    'VAL': [-1.0],
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
               "LOAD": [forcesmass],
               "BOUNDCOND": [bc1],
    },
}
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'michell_2ss', 'show': True, 'scale': 10, 'savepng': True, 'lines': True,
                        # 'plottags': {
                        # # 'line': True
                        # 'point': True
                        #      }
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# sys.exit()

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'mode',  # mode, freq, time ...
                        'start': 0,
                        'end': 6,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

postprocdata = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    'structural': {'modes': True},
                    # 'structural': {'frf': True}
                    },
                "PLOTSET": {'filename': 'test_dynamic', 'savepng': True},
                "TRACKER": {'frf': {'x': 16, 'y': 1, 'z': 0, 'dof':0}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }

postporc_result = fea.PostProcess(postprocdata)