import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import SteadyStateLinear
import numpy as np

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(SteadyStateLinear, 'sim')

mat = {
    "NAME": "aco",
    "VXY": 0.3,
    "EXX": 210E3,       # N/mm^2 --> MPa
    }

geo = {
    "NAME": "espessura",
    "THICKN": 100,
    }

modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'mesh_gmsh',
        "meshimport": {'object': 'nafems_le1_mesh_quad4_fine'}, 
        'meshconfig': {
            'mesh': 'quad4',
            }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 2,
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

press = {
    'TYPE': 'forceedge',
    'DOF': 'pressure',
    'DIR': 'line',
    'TAG': 2,
    'VAL': [10],
    }

bc1 = {
    'TYPE': 'fixed',  
    'DOF': 'ux',
    'DIR': 'line',
    'TAG': 1,
    }


bc2 = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'line',
    'TAG': 3,
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [press],
               "BOUNDCOND": [bc1, bc2],
    },
}
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
                         'plottags': {'line': True}
                         },
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

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'output', 'savepng': True},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

STR_YY_D = max(postprocdata['STRESS_YY'])
print('STRESS_YY_D ',STR_YY_D)

STR_YY_D_REF = 92.7
erro = (abs(STR_YY_D_REF - STR_YY_D)/abs(STR_YY_D_REF))*100
print('ERRO [%] ',erro)