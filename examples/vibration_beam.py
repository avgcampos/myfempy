from myfempy import newAnalysis
from myfempy import DynamicEigenLinear

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(DynamicEigenLinear)
LX = 10 
nelx = 10
# MODEL SET
mat = {
    "NAME": "aco",
    "EXX": 10000,       # N/mm^2 --> MPa
    "GXY": 4000,
    "RHO": 10,    # kg/mm^2
    }

geo = {
    "NAME": "secao",
    "AREACS": 10,
    "INERYY": 100,
    "INERZZ": 100,
    "INERXX": 100,
    }


# MODEL SET
# gmsh config
points = [
    [0, 0, 0],
    [LX, 0, 0],
]
         
lines = [
    [1, 2],
         ]

modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'test_line_force_x',
        'pointlist': points,
        'linelist': lines,
        'meshconfig': {
            'mesh': 'line3',
            "numbernodes": nelx + 1,
            }
    },

    "ELEMENT": {
        'TYPE': 'structbeam',
        'SHAPE': 'line3',
        # 'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'uniaxialstress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat],
    },
    
    "GEOMETRY": {
        "GEO": 'frame',
        "SECTION": 'userdefined',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': LX, 'y': 0, 'z': 0},
    'VAL': [-10],
    }

f1b = {
    'TYPE': 'forcebeam',
    'DOF': 'fy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'VAL': [-100],
    }

f2 = {
    'TYPE': 'forcenode',
    'DOF': 'tz',
    'DIR': 'node',
    'LOC': {'x': 10, 'y': 10, 'z': 0},
    'VAL': [5000],
    }

bc1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 0, 'z': 0},
    }

bc2 = {
    'TYPE': 'fixed',
    'DOF': 'uz',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [],
               "BOUNDCOND": [bc1],
    },
}
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'beam', 'show': True, 'scale': 20, 'savepng': True, 'lines': False,
                         },
              }
fea.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'mode',  # mode, freq, time ...
                        'start': 0,
                        'end': 30,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

print(solverdata['solution']['FREQ'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'modes': True}},
                "PLOTSET": {'show': True, 'filename': 'test_shakedown', 'savepng': True},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)
