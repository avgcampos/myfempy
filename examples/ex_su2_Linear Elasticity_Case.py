'''
SU2 LINEAR ELASTICITY CASE - TEST BENCH
'''

from myfempy.mesh.genmesh import ModelGen
from myfempy.core.solver import Solver
from myfempy.postprc.postcomp import PostProcess
from myfempy.plots.postplot import postproc_plot
from myfempy.plots.prevplot import preview_plot
import time


# Case: https://su2code.github.io/tutorials/Linear_Elasticity/
#Sanchez, R. Et Al. (2016), Towards A Fluid-Structure Interaction Solver For Problems With Large Deformations Within The Open-Source Su2 Suite, 57th Aiaa/Asce/Ahs/Asc Structures, Structural Dynamics, And Materials Conference, 4-8 January 2016, San Diego, California, Usa. Doi: 10.2514/6.2016-0205


#----------------------------- PRE-PROCESS -----------------------------------#
mat = {
    "NAME": "material",
    "VXX": 0.35,
    "EXX": 5.0E9,
    "MAT": 'isotropic',
    "DEF": 'planestress'
}

geo = {"NAME": "geo1", "THICKN": 1.0}

force = {'DEF': 'forceedge',
         'DOF': 'fx',
         'DIR': 'edgex',
         'LOC': {'x': 0, 'y': 999, 'z': 0},
         'VAL': [1.0],
         }

bondcond = {'DEF': 'fixed',
            'DOF': 'all',
            'DIR': 'edgey',
            'LOC': {'x': 999, 'y': 0, 'z': 0},
            }

#-------------------------- GEN WITH LEGACY MESH ----------------------------#
meshdata = {"LEGACY": {'lx': 0.5, 'ly': 10, 'mesh': 'quad4', 'elem': 'plane41', 'nx': 10, 'ny': 100},
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            "FORCES": [force],
            "BOUNDCOND": [bondcond],
            "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
            "DOMAIN":'structural'
            }


# points = [[0, 0, 0],
#           [4000, 0, 0],
#           [4000, 1000, 0],
#           [0, 1000, 0]]

# lines = [[1, 2],
#          [2, 3],
#          [3, 4],
#          [4, 1]]

# plane = [[1, 2, 3, 4]]

# meshdata = {"GMSH": {'filename': 'malha',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'quad4',
#                          'elem': 'plane41',
#                          'sizeelement': 100,
#                          'meshmap': {'on': True,
#                                      'edge': 'all'}}},
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [force],
#             "BOUNDCOND": [bondcond],
#             "QUADRATURE": {'meth': 'gaussian', 'npp': 4}
#             }

#-------------------------------- GEN MESH -----------------------------------#
modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'malha_triagular_pre', 'show': True, 'scale': 4, 'savepng': False, 'lines': False, 'plottags': {'point': False}},
              'QUALITY': {'show': False, 'method': 1, 'scale': 0.1},
              'LABELS': {'show': False, 'lines': True, 'scale': 1},
              }

preview_plot(previewset, modelinfo)
#-------------------------------- SOLVER -------------------------------------#
solverset = {"SOLVER": 'SLI',
             'TOL': 1E-8,
             "STEPSET": {'type': 'table',  # mode, freq, time ...
                         'start': 0,
                         'end': 1,
                         'step': 1},
             #"TRACKER": {'show': False, 'result2plot': 'displ', 'max': []}
             }

solution = Solver.get_static_solve(solverset, modelinfo)
#----------------------------- POST-PROCESS ----------------------------------#
postprocset = {"SOLUTION": solution,
            "COMPUTER": {'elasticity': {'displ': True, 'stress': True, 'average': True}},
                # 'step':2
            "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'su2_test', 'savepng': False},
            # "TRACKER": {'show': True, 'result2plot':'stress', 'point': {'x':6,'y':1.5,'z':0}}
            }

postporc_result = PostProcess(modelinfo).compute(postprocset)
#----------------------------- VIEW SOLUTION ---------------------------------#
postproc_plot(postprocset, postporc_result, modelinfo)
