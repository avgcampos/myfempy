import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import DynamicEigenLinear

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(DynamicEigenLinear)

# MODEL SET
mat = {
    "NAME": "Aluminum_Alloy",
    "VXY": 0.33,
    "EXX": 71E6,        # MPa (N/mm2)
    "RHO": 2.77E-06,    # kg/mm3
    }

geo = {"NAME": "Solid"}

# MODEL SET
LX = 100
LY = 100
esize = 20

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
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'mesh_gmsh',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'hexa8',   #quad4 tria3
            'sizeelement': 2*esize,
            'extrude': 100,
            'meshmap': {'on': True,
                        'edge': 'all', #[[1,3], [2,4]], #'all'
                        # "numbernodes": [20, 10],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structsolid',
        'SHAPE': 'hexa8',
        'INTGAUSS': 9,
    },

    "MATERIAL": {
        "MAT": 'solidelastic',
        "TYPE": 'isotropic',
        "PROPMAT": [mat],
    },

    "GEOMETRY": {
        "GEO": 'solid',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [],
               "BOUNDCOND": [],
    },
}
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'preview', 'show': True, 'scale': 8, 'savepng': True, 'lines': True,
                        #  'plottags': {'point': True},
                         },
              }
fea.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 12,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

print(solverdata['solution']['FREQ'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    'structural': {'modes': True},
                    },
                "PLOTSET": {'filename': 'output', 'savepng': True},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }
postprocdata = fea.PostProcess(postprocset)