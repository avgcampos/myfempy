import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import DynamicHarmonicResponseLinear

#===============================================================================
#                                   FEA
#===============================================================================
fea = newAnalysis(DynamicHarmonicResponseLinear)

mat1 = {
    "NAME": "mat1",
    "VXY": 0.25,
    "EXX": 250E9,
    "RHO": 7800
    }


geo = {
    "NAME": "espessura",
    "THICKN": 1.0,
    }

# MODEL SET
modeldata = {
        
    "MESH": {
        'TYPE': 'legacy',
        'LX': 16,
        'LY': 2,
        'NX': 4*16,
        'NY': 4*2,
        },
        
    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        # 'INTGAUSS': 4,
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
    'DIR': 'node',
    'LOC': {'x': 20, 'y': 5, 'z': 0},
    # 'TAG': 4,
    'VAL': [0],
    }

f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 16, 'y': 1, 'z': 0},
    'VAL': [-1.0],
    }
bc1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [f1],
               "BOUNDCOND": [bc1],
    },
}

fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'vibration_beam', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
              }
fea.PreviewAnalysis(previewset)
# #-------------------------------- SOLVER -------------------------------------#
solverset = {
            "STEPSET": {
                'type': 'freq',
                'start': 0,
                'end': 200,
                'step': 0.5},
             }
solverdata = fea.Solve(solverset)

postprocdata = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    'structural': {'frf': True}
                    },
                "PLOTSET": {'filename': 'test_dynamic', 'savepng': True},
                "PLOT": {'freq': {'x': 16, 'y': 1, 'z': 0, 'dof':2}},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }
postporc_result = fea.PostProcess(postprocdata)