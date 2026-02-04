
"""
==============================================================================
MYFEMPY example runscript
==============================================================================
@File    :   run_fist_test.py
@Date    :   2026/02/03
@Author  :   Antonio V. G. Campos
@Description : This example can be used as an initial analysis to test the complete
installation of the myfempy code. A report (PT-BR) is available to verify the program's
compatibility. See the user guide or _help_ for more information about commands.
"""

# environ['OMP_NUM_THREADS'] = '1'        #win

# ===============================================================================
# imports, set SteadyStateLinear as a Solver or SteadyStateLinearIterative, see help
# ===============================================================================
from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative
import numpy as np

# ===============================================================================
# FEA API
# ===============================================================================
fea = newAnalysis(SteadyStateLinear)

# ===============================================================================
# set material prop.
# ===============================================================================
mat = {
    "NAME": "aco",
    "VXY": 0.30,
    "EXX": 200E3,       # N/mm^2 --> MPa
    "GXY": 0.5*200E3,
    "RHO": 7.85E-06,    # kg/mm^2
    }

# ===============================================================================
# set geometry prop.
# ===============================================================================
geo = {
    "NAME": "sec",
    "THICKN": 40,
    # "DIM":{
    #     "b":203,
    #     "h":203,
    #     "t":25.4,
    #     "d":0,
    # },
    }

# ===============================================================================
# set gmsh config.
# ===============================================================================
# MODEL SET
LX = 1200
LY = 80
esize = 5
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


# ===============================================================================
# set modeldata
# ===============================================================================
modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'tests',
        # "meshimport": 'object_dir',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        # 'arc': arcs,
        'meshconfig': {
            'mesh': 'hexa8',   #quad4 tria3
            "numbernodes": 2,
            'sizeelement': 2*esize,
            'extrude': 40,
            'meshmap': {'on': True,
                        'edge': 'all', #[[1,3], [2,4]], #'all'
                        # "numbernodes": [20, 10],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structsolid',
        'SHAPE': 'hexa8',
        # 'INTGAUSS': 8,
    },

    "MATERIAL": {
        "MAT": 'uniaxialstress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat],
    },
    
    "GEOMETRY": {
        "GEO": 'solid',
        # "SECTION": 'rectangle',
        "PROPGEO": [geo],
    },
}

# ===============================================================================
# pass modeldata to Model API
# ===============================================================================
fea.Model(modeldata)

# ===============================================================================
# set bound. cond.
# ===============================================================================
bc = {
    'TYPE': 'fixed',  
    'DOF': 'full',
    'DIR': 'point',
    'TAG': 1,
    }

# ===============================================================================
# set load
# ===============================================================================
fr = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'point',
    'TAG': 2,
    'VAL': [1.5625],
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [fr],
               "BOUNDCOND": [bc],
    },
}

# ===============================================================================
# pass physicdata to Physic API
# ===============================================================================
fea.Physic(physicdata)

# ===============================================================================
# preview config.
# ===============================================================================
previewset = {'RENDER': {'filename': 'plane_stress', 'show': True, 'scale': 5, 'savepng': True, 'lines': False,
                         'plottags': {'point': True},
                        #  'cs': True,
                         },
              }
fea.PreviewAnalysis(previewset)

# ===============================================================================
# SOLVER
# ===============================================================================
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':False,
             'MP':False,
            }
solverdata = fea.Solve(solverset)

# ===============================================================================
# post process
# ===============================================================================
postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'test', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}}, #'coord': True
            }
postprocdata = fea.PostProcess(postprocset)

print(np.max(np.abs(solverdata['solution']['U'])))