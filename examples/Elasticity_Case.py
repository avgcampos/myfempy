from myfempy import newAnalysis
from myfempy import StaticLinear

from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================

fea = newAnalysis(StaticLinear)

mat = {
    "NAME": "material",
    "VXX": 0.35,
    "EXX": 5.0E9,
    "MAT": 'isotropic',
    "DEF": 'planestress'
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
        'filename': 'data_mesh',
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
        'TYPE': 'plane',
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
    "DOMAIN": 'structural',
    "LOAD": [f1],
    "BOUNDCOND": [bc1],
    }

fea.Physic(physicdata)


previewset = {'RENDER': {'filename': 'model_preview', 'show': True, 'scale': 2, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
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

# print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "PLOTSET": {'show': True, 'filename': 'solution', 'savepng': True},
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "TRACKER": {'point': {'x': 0.5*LX, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True, 'inci': True, 'coord': True}},
            }
postprocdata = fea.PostProcess(postprocset)
