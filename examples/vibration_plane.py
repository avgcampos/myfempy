from myfempy import newAnalysis
from myfempy import DynamicEigenLinear

#===============================================================================
#                                   FEA
#===============================================================================
fea = newAnalysis(DynamicEigenLinear)

mat = {
    "NAME": "mat1",
    "VXY": 0.25,
    "EXX": 10000,
    "RHO": 10,
    }

geo = {
    "NAME": "espessura",
    "THICKN": 2,
    }

# MODEL SET
LX = 10
LY = 5

# nelx = 40
# nely = 30

esize = 0.25

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
    #     'TYPE': 'legacy',
    #     'LX': LX,
    #     'LY': LY,
    #     'NX': nelx,
    #     'NY': nely,
    #     },
        
    "MESH": {
            'TYPE': 'gmsh',
            'filename': 'bench_vibra',
            'pointlist': points,
            'linelist': lines,
            'planelist': plane,
            'meshconfig': {
                'mesh': 'quad4',   #quad4
                'sizeelement': 2 * esize,
                'meshmap': {'on': True,
                            'edge': 'all', #'all'
                        #  "numbernodes": 10,
                }}
        },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        # 'INTGAUSS': 4,
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
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
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
                "BOUNDCOND": [bc_node_left],
    },
    }
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'bench_vibra_preview', 'show': True, 'scale': 10, 'savepng': True, 'lines': False,
                        'plottags': {
                        # 'line': True
                        'point': True
                              }
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'mode',  # mode, freq, time ...
                        'start': 0,
                        'end': 3,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

print(solverdata['solution']['FREQ'])

postprocdata = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    'structural': {'modes': True},
                    },
                "PLOTSET": {'filename': 'bench_vibra', 'savepng': True},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }

postporc_result = fea.PostProcess(postprocdata)