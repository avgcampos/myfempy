# -*- coding: utf-8 -*-

import sys

# from myfempy.io.ioctrl import postproc2vtk
from myfempy.mesh.genmesh import ModelGen
# from myfempy.solver.assembler import stiffness, loads
# from myfempy.solver.solverset import get_constrains_dofs
from myfempy.core.solver import gen_static_solution
from myfempy.postprc.postcomp import PostComputer, postproc_plot
from myfempy.plots.prevplot import preview_plot
from myfempy.bin.tools import print_console

#----------------------------- PRE-PROCESS -----------------------------------#
# print_console('pre')


mat1 = {
    "NAME": "steel",
    "VXX": 0.33,
    "EXX": 200E3,
    "MAT": 'isotropic',
    "DEF": 'planestress'
}


mat2 = {
    "NAME": "mola",
    "DAMP": 0,
    "STIF": 1000,
    "MAT": 'springlinear',
    "DEF": 'lumped'
}


geo1 = {"NAME": "geo1", "THICKN": 100}

geo2 = {"NAME": "geo2"}

f1 = {'DEF': 'forceedge',
      'DOF': 'fx',
      'DIR': 'edge',
      'TAG': 2,
      # 'LOC': {'x':800,'y':999,'z':0},
      'VAL': [100],  # table: [val_setp1, val_step2,...]
      }


f2 = {'DEF': 'forceedge',
      'DOF': 'fy',
      'DIR': 'edge',
      'TAG': 3,
      # 'LOC': {'x':800,'y':999,'z':0},
      'VAL': [-100],  # table: [val_setp1, val_step2,...]
      }


s1 = {'DEF': 'forcenode',
      'DOF': 'spring2ground',
      'DIR': 'edgey',
      # 'TAG': 9,
      'LOC': {'x': 999, 'y': 3, 'z': 0},
      'VAL': [100],  # table: [val_setp1, val_step2,...]
      }


# f2 = {'DEF':'forcenode',
#           'DIR': 'point',
#           'DOF': 'fy',
#           'LOC': {'x':0,'y':3,'z':0},
#           'VAL': [-1000, 0], # table: [val_setp1, val_step2,...]
#           }


# f3 = {'DEF':'forcenode',
#           'DIR': 'point',
#           'DOF': 'fy',
#           'LOC': {'x':6,'y':3,'z':0},
#           'VAL': [-1000,-2000,-3000,-4000,-5000], # table: [val_setp1, val_step2,...]
#           }


# bc1 = {'DEF': 'fixed',
#        'DOF': 'all',
#        'DIR': 'pointx',
#        # 'TAG':4,
#        'LOC': {'x': 0, 'y': 999, 'z': 0},
#        }

bc1 = {'DEF': 'fixed',
       'DOF': 'all',
       'DIR': 'edge',
       'TAG': 2,
       # 'LOC': {'x':0,'y':999,'z':0},
       }


bc2 = {'DEF': 'fixed',
       'DOF': 'all',
       'DIR': 'edge',
       'TAG': 4,
       # 'LOC': {'x':0,'y':999,'z':0},
       }

# bc1 = {'DEF': 'fixed',
#        'DOF': 'all',
#        'DIR': 'point',
#        'TAG': 1,
#        # 'LOC': {'x':0,'y':999,'z':0},
#        }

# bc2 = {'DEF': 'fixed',
#        'DOF': 'all',
#        'DIR': 'point',
#        'TAG': 2,
#        # 'LOC': {'x':0,'y':999,'z':0},
#        }

#-------------------------------- BY ADD ELEMENTS -----------------------------------#

# elementos = [[1,"plane31","steel","geo1",[1, 6, 5]],
#             [2,"plane31","steel","geo1",[1, 2, 6]],
#             [3,"plane31","steel","geo1",[2, 7, 6]],
#             [4,"plane31","steel","geo1",[2, 3, 7]],
#             [5,"plane31","steel","geo1",[3, 8, 7]],
#             [6,"plane31","steel","geo1",[3, 4, 8]],
#             [7,"plane31","steel","geo1",[5, 10, 9]],
#             [8,"plane31","steel","geo1",[5, 6, 10]],
#             [9,"plane31","steel","geo1",[6, 11, 10]],
#             [10,"plane31","steel","geo1",[6, 7, 11]],
#             [11,"plane31","steel","geo1",[7, 12, 11]],
#             [12,"plane31","steel","geo1",[7, 8, 12]]]


# coordenadas = [[1, 0, 0, 0],
#                 [2, 2, 0, 0],
#                 [3, 4, 0, 0],
#                 [4, 6, 0, 0],
#                 [5, 0, 1.5, 0],
#                 [6, 3, 1.5, 0],
#                 [7, 4, 1.5, 0],
#                 [8, 6, 1.5, 0],
#                 [9, 0, 3, 0],
#                 [10, 2, 3, 0],
#                 [11, 4, 3, 0],
#                 [12, 6, 3, 0]]


# meshdata = {"ELEMLIST": elementos,
#             "NODELIST": coordenadas,
#             "PROPMAT": [mat1],
#             "PROPGEO": [geo1],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1, bc2],
#         }


#-------------------------------- BT DEF GEOMETRY -----------------------------------#
# elementos = [[1,'spring21','mola','geo1',[1,2]]]

# coordenadas = [[1, 6, 0, 0],
#                [2, 7, 0, 0]]


# meshdata = {"LEGACY": {'lx':800,'ly':400,'mesh':'tria3', 'elem': 'plane31','nx':32, 'ny':16},
#             "PROPMAT": [mat1],
#             "PROPGEO": [geo1],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1],
#             "QUADRATURE":{'meth':'gaussian','npp':3}
#             }


#-------------------------------- BY GMSH -----------------------------------#

points = [[0, 0, 0],
          [4000, 0, 0],
          [4000, 1000, 0],
          [0, 1000, 0]]
# [5,0,0],
# [5,1,0],
# ]


lines = [[1, 2],
         [2, 3],
         [3, 4],
         [4, 1]]
# [2,5],
# [5,6],
# [6,3],
# ]


plane = [[1, 2, 3, 4]]
# [5,6,7,2]]

# # hole = [planeN, dia, [center coord], [arc lenght]]
# # hole = [[1, 0.2, [2, 0.5, 0], ['0', '2*Pi']],
# #         [1, 0.1, [1, 0.5, 0], ['0', '2*Pi']],
# #         [1, 0.1, [3, 0.5, 0], ['0', '2*Pi']]]
# #         # [2, 0.1, [4.5, 0.5, 0], ['0', '2*Pi']]]

# hole1 =   [['hole', 1, 0.1, [4.5, 0.5, 0]]]


# points = [[0,0,0],
#           [0.5,0,0],
#           [0.6464466094067263,0.3535533905932738,0],
#           [0,1,0],
#           [1,0.5,0],
#           [1,1,0],
#           [1,0,0],
#           ]


# lines = [[1,2],
#           [3,4],
#           [4,1],
#           [5,6],
#           [6,4],
#           ]


# hole2 =   [['arc', [2, 7, 3]],
#             ['arc', [3, 7, 5]]]


# plane = [[1,6,2,3],
#           [7,4,5,2]]

meshdata = {"GMSH": {'filename': 'test_gmsh',
                     'pointlist': points,
                     'linelist': lines,
                     'planelist': plane,
                     # 'curve': hole2,
                     'meshconfig': {
                         'mesh': 'quad4',
                         'elem': 'plane41',
                         'sizeelement': 100,
                         'meshmap': {'on': True,
                                     'edge': 'all'}}},  # 'numbernodes':50
            "PROPMAT": [mat1],
            "PROPGEO": [geo1],
            "FORCES": [f1],
            "BOUNDCOND": [bc2],
            "QUADRATURE": {'meth': 'gaussian', 'npp': 4}
            }


#-------------------------------- GEN MESH -----------------------------------#
print_console('mesh')

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'test_gmsh_pre', 'show': True, 'scale': 4, 'savepng': True, 'lines': False, 'plottags': {'point': True}},
              'QUALITY': {'show': False, 'method': 1, 'scale': 0.1},
              'LABELS': {'show': False, 'lines': True, 'scale': 1},
              }

preview_plot(previewset, modelinfo)

# sys.exit()
# #
#-------------------------------- SOLVER -------------------------------------#
print_console('solver')

solverset = {"SOLVER": 'SLD',
             'TOL': 2E-3,
             "STEPSET": {'type': 'table',  # mode, freq, time ...
                         'start': 0,
                         'end': 1,
                         'step': 1},
             "TRACKER": {'show': False, 'result2plot': 'displ', 'max': []}
             }


solution = gen_static_solution(solverset, modelinfo)

#----------------------------- POST-PROCESS ----------------------------------#

print_console('post')
postprocset = {
    "SOLUTION": solution['U'],
    "COMPUTER": {'displ': True, 'stress': True, 'average': True},
    "SCALE": 20,
    # 'step':2
    "PLOTSET": {'show': True, 'result2plot': {'displ': []}, 'filename': 'quadrilatero_load_edge', 'savepng': True},
    "TRACKER": {'show': False, 'result2plot': 'stress', 'point': {'x': 6, 'y': 1.5, 'z': 0}}
}


postporc_result = PostComputer(modelinfo).main(postprocset)

#----------------------------- VIEW SOLUTION ---------------------------------#
postproc_plot(postprocset, postporc_result, modelinfo)

print_console('thank')
