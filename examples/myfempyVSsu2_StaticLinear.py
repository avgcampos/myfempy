
import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import SteadyStateLinear

import numpy as np

# ===============================================================================
#                                   FEA
# ===============================================================================

fea = newAnalysis(SteadyStateLinear, 'su2_test')

mat = {
    "NAME": "material",
    "VXX": 0.35,
    "EXX": 5.0E9,
    }

geo = {"NAME": "geo1",
       "THICKN": 1.0}

# MODEL SET
LX = 0.5
LY = 10 

# gmsh config
points = [
        [0, 0, 0],
        [LX, 0, 0],
        [LX, LY, 0],
        [0, LY, 0]
        ]
         
lines = [
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 1],
        ]

plane = [
        [1, 2, 3, 4]
        ]

modeldata = {

   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'mesh_gmsh',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'quad4',
            'sizeelement': 0.2,
            'meshmap': {'on': True,
                        'edge': 'all',
                    #  "numbernodes": 10,
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

f1 = {
    'TYPE': 'forceedge',
    'DOF': 'fx',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'VAL': [1.0],
    }

bc1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [f1],
               "BOUNDCOND": [bc1],
    },
    }

fea.Physic(physicdata)


previewset = {'RENDER': {'filename': 'model_preview', 'show': True, 'scale': 2, 'savepng': True, 'lines': False,
                         'plottags': {'line': True}
                         },
              }
fea.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
             }
solverdata = fea.Solve(solverset)

print(np.max(np.abs(solverdata['solution']['U'])))

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'solution', 'savepng': True},
                "REPORT": {'log': True,
                           'get': {
                                'nelem': True,
                                'nnode': True
                                }
                                },
            }
postprocdata = fea.PostProcess(postprocset)
