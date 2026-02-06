from myfempy import newAnalysis
from myfempy import SteadyStateLinearIterative

from time import time

fea = newAnalysis(SteadyStateLinearIterative)
# MODEL SET
mat = {
    "NAME": "mat",
    "EXX": 30000,       # N/mm^2 --> MPa
    "GXY": 10000,
    }

geo = {
    "NAME": "secao1",
    "AREACS": 10,
    "INERYY": 100,
    "INERZZ": 100,
    "INERXX": 50,
    # "CG":{
    #     "y_max":1,
    #     "y_min":-1,
    #     "z_max":1,
    #     "z_min":-1,
    #     "r_max":1
    # },
    }

# gmsh config
points = [
    [100, 0, 100],
    [100, 100, 100],
    [0, 100, 100],
    [100, 100, 0],
]
         
lines = [[2, 1],
         [3, 2],
         [2, 4],
         ]

modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'mesh_gmsh',
        'pointlist': points,
        'linelist': lines,
        'meshconfig': {
            'mesh': 'line2',
            "numbernodes": 2,
            }
    },

    "ELEMENT": {
        'TYPE': 'structbeam',
        'SHAPE': 'line2',
        # 'INTGAUSS': 2,
    },

    "MATERIAL": {
        "MAT": 'uniaxialstress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat, mat, mat],
    },
    
    "GEOMETRY": {
        "GEO": 'frame',
        "SECTION": 'userdefined',
        "PROPGEO": [geo, geo, geo],
    },
}
fea.Model(modeldata)

bc1 = {
    'TYPE': 'fixed',  
    'DOF': 'full',
    'DIR': 'point',
    'TAG': 1,
    }

bc2 = {
    'TYPE': 'fixed',  
    'DOF': 'full',
    'DIR': 'point',
    'TAG': 3,
    }

bc3 = {
    'TYPE': 'fixed',  
    'DOF': 'full',
    'DIR': 'point',
    'TAG': 4,
    }

fy = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'point',
    'TAG': 2,
    'VAL': [-50],
    }

tx = {
    'TYPE': 'forcenode',
    'DOF': 'tx',
    'DIR': 'point',
    'TAG': 2,
    'VAL': [-1000],
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [fy, tx],
               "BOUNDCOND": [bc1, bc2, bc3],
    },
}
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'preview', 'show': True, 'scale': 15, 'savepng': True, 'lines': True,
                         'plottags': {'point': True},
                        # 'cs': True,
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)
print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'output', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "REPORT": {'log': True, 'get':{
                        'nelem': True,
                        'nnode': True,
                        'inci': True,
                        'coord':True,
                        'tabmat':True,
                        'tabgeo':True,
                        'boundcond_list':True,
                        'forces_list':True,
                    },
            }}
postprocdata = fea.PostProcess(postprocset)