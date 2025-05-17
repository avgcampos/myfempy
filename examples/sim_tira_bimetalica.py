
from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative

from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================

fea = newAnalysis(SteadyStateLinear) # StaticLinearIterative FASTER THEN StaticLinear

steel = {
    "NAME": "steel",
    "KXX": 0.0605,	    # W/mm·°C  
    "KYY": 0.0605,
    "VXX": 0.30,
    "EXX": 200E3,       # N/mm^2 --> MPa
    "RHO": 7.85E-06,    # kg/mm^3
    "CTE": 12E-6        # (mm/mm) /°C
    }

copper = {
    "NAME": "copper",
    "KXX": 0.4230,	    # W/mm·°C  
    "KYY": 0.4230,
    "VXX": 0.33,
    "EXX": 110E3,       # N/mm^2 --> MPa
    "RHO": 7.85E-06,    # kg/mm^3
    "CTE": 17E-6        # (mm/mm) /°C
    }

geo = {
    "NAME": "espessura",
    "THICKN": 5, # mm
    }

points = [
    [0, 0, 0],
    [200, 0, 0],
    [0, 10, 0],
    [200, 10, 0],
    [0, 20, 0],
    [200, 20, 0]
]
         
lines = [[1, 2],
         [2, 4],
         [3, 4],
         [3, 1],
         [6, 4],
         [5, 6],
         [5, 3],
         ]

plane = [[1, 2, 3, 4], [3, 5, 6, 7]]

modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'test_line_force_x',
        # "meshimport": 'object_dir',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        # 'arc': arcs,
        'meshconfig': {
            'mesh': 'quad4',   #quad4
            'sizeelement': 1,
            'meshmap': {'on': True,
                        'edge': [[1,3,6], [2,4,5,7]], #'all'
                        "numbernodes": [100, 8],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'heatplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'heatplane',
        "TYPE": 'isotropic',
        "PROPMAT": [copper, steel],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo, geo],
    },
}
fea.Model(modeldata)

# # sys.exit()
q = -0.1 #1      # W/mm2
# Q = 0.2         # W/mm3
h = 0.000005    # W/mm2.ºC 0.000005 (air)
Tinf = 20       # ºC

heatflux = {
    'TYPE': 'heatfluxedge',
    'DOF': 'heatflux',
    'DIR': 'edgex',
    'LOC': {'x': 200, 'y': 999, 'z': 0},
    'VAL': [0*q],
    }

# heatgen = {
#     'TYPE': 'heatgeneration',
#     'DOF': 'heatflux',
#     'DIR': 'node',
#     # 'LOC': {'x': 0, 'y': 999, 'z': 0},
#     'VAL': [Q],
#     }

convc_s1 = {
    'TYPE': 'convectionedge',
    'DOF': 'convection',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'VAL': [-h*Tinf],
    }

convc_s2 = {
    'TYPE': 'convectionedge',
    'DOF': 'convection',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 20, 'z': 0},
    'VAL': [-h*Tinf],
    }

# convc_s3 = {
#     'TYPE': 'convectionedge',
#     'DOF': 'convection',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': LY, 'z': 0},
#     'VAL': [h*Tinf],
#     }

# bc1 = {
#     'TYPE': 'insulated',  
#     'DOF': 'full',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }

# bc2 = {
#     'TYPE': 'insulated',
#     'DOF': 'full',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 0, 'z': 0},
#     }

# bc3 = {
#     'TYPE': 'insulated',
#     'DOF': 'full',
#     'DIR': 'edgex',
#     'LOC': {'x': LX, 'y': 999, 'z': 0},
#     }

bcNH1 = {
    'TYPE': 'temperature',
    'DOF': 't',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'VAL': [200.0],
    }

bcNH2 = {
    'TYPE': 'temperature',
    'DOF': 't',
    'DIR': 'edgex',
    'LOC': {'x': 200, 'y': 999, 'z': 0},
    'VAL': [200.0],
    }

# bcNH3 = {
#     'TYPE': 'temperature',
#     'DOF': 't',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 0, 'z': 0},
#     'VAL': [100.0],
#     }

physicdata = {
    "PHYSIC": {"DOMAIN": "thermal",           # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [heatflux],
               "BOUNDCOND": [bcNH1],
    }
}
fea.Physic(physicdata)


previewset = {'RENDER': {'filename': 'heat_transfer', 'show': True, 'scale': 5, 'savepng': True, 'lines': False,
                         'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# sys.exit()

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

# print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'thermal': {'temp': True, 'heatflux': True}},
                "PLOTSET": {'show': True, 'filename': 'SIM_TIRA_BIMETALICA_HEAT', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

### new sim

# sys.exit()

modeldata = {
    "MESH": {
    'TYPE': 'gmsh',
    'filename': 'test_line_force_x',
    # "meshimport": 'object_dir',
    'pointlist': points,
    'linelist': lines,
    'planelist': plane,
    # 'arc': arcs,
    'meshconfig': {
        'mesh': 'quad4',   #quad4
        'sizeelement': 2,
        'meshmap': {'on': True,
                    'edge': [[1,3,6], [2,4,5,7]], #'all'
                    "numbernodes": [100, 8],
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
        "PROPMAT": [copper, steel],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo, geo],
    },
}
fea.Model(modeldata)

# fc = {
#     'TYPE': 'forcenode',
#     'DOF': 'fy',
#     'DIR': 'node',
#     'LOC': {'x': 200, 'y': LY, 'z': 0},
#     'VAL': [0.0],
#     }

fc = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'edgex',
    'LOC': {'x': 200, 'y': 999, 'z': 0},
    'VAL': [0.0],
    }

# bcb = {
#     'TYPE': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 0, 'z': 0},
#     }

bcl = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }

# bcr = {
#     'TYPE': 'fixed',
#     'DOF': 'full',
#     'DIR': 'edgex',
#     'LOC': {'x': LX, 'y': 999, 'z': 0},
#     }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [fc],
               "BOUNDCOND": [bcl],
    },
    "COUPLING": {"TYPE": 'thermalstress', # fsi asi
                 "POST": [postprocdata]
    },
    }
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'structural', 'show': True, 'scale': 10, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
# fea.PreviewAnalysis(previewset)

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
                "PLOTSET": {'show': True, 'filename': 'SIM_TIRA_BIMETALICA_DISP', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)