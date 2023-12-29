'''
myfempy Tutorial 03

Geração da malha com gmsh

Necessario instalação prévia do gmsh (não nativo do myfempy)

'''

# Imports

from myfempy import ModelGen
from myfempy import Solver
from myfempy import PostProcess
from myfempy import postproc_plot
from myfempy import preview_plot

# Definição do material, geometria 

mat = {
    "NAME": "material",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'planestress'
    }

geo = {
    "NAME": "geometria",
    "THICKN": 1.0
    }

# ----------------------------------------------------------

# Placa L estado plano (malha quad 4) 100 x 100 mm

# f1 = {
#     'DEF': 'forceedge',
#     'DOF': 'fy',
#     'DIR': 'edge',
#     'TAG': 3,
#     'VAL': [-1000.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'edge',
#     'TAG': 5,
#     }


# points = [[0, 0, 0],
#           [100, 0, 0],
#           [100, 50, 0],
#           [50, 50, 0],
#           [50, 100, 0],
#           [0, 100, 0],
#           ]

# lines = [[1, 2],
#          [2, 3],
#          [3, 4],
#          [4, 5],
#          [5, 6],
#          [6, 1],
#          ]

# plane = [[1, 2, 3, 4, 5, 6]]

# meshdata = {"GMSH": {'filename': 'tutorial_03a',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'tria3',   #tria3
#                          'elem': 'plane31', #plane31
#                          'sizeelement': 10,
#                          'meshmap': {'on': True,
#                                      'edge': 'all',
#                                     #  "numbernodes": 10,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1],
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_03a',
#                          'show': True,
#                          'scale': 5,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {'point': True,
#                                       'edge': True}},
#               }

# preview_plot(previewset, modelinfo)


# --------------------------------------------------

# # Placa H estado plano (malha quad 4) 100 x 100 mm

f1 = {
    'DEF': 'forceedge',
    'DOF': 'fx',
    'DIR': 'edge',
    'TAG': 8,
    'VAL': [100.0],
    }

f2 = {
    'DEF': 'forceedge',
    'DOF': 'fx',
    'DIR': 'edge',
    'TAG': 10,
    'VAL': [-100.0],
    }

f3 = {
    'DEF': 'forceedge',
    'DOF': 'fx',
    'DIR': 'edge',
    'TAG': 4,
    'VAL': [100.0],
    }

f4 = {
    'DEF': 'forceedge',
    'DOF': 'fx',
    'DIR': 'edge',
    'TAG': 2,
    'VAL': [-100.0],
    }

bc1 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'edge',
    'TAG': 3,
    }

bc2 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'edge',
    'TAG': 9,
    }

points = [
    [0, 0, 0],
    [20, 0, 0],
    [20, 30, 0],
    [80, 30, 0],
    [80, 0, 0],
    [100, 0, 0],
    [100, 100, 0],
    [80, 100, 0],
    [80, 70, 0],
    [20, 70, 0],
    [20, 100, 0],
    [0, 100, 0],
    ]

lines = [[1, 2],
         [2, 3],
         [3, 4],
         [4, 5],
         [5, 6],
         [6, 7],
         [7, 8],
         [8, 9],
         [9, 10],
         [10, 11],
         [11, 12],
         [12, 1],
         ]

plane = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]

meshdata = {"GMSH": {'filename': 'tutorial_03b',
                     'pointlist': points,
                     'linelist': lines,
                     'planelist': plane,
                     'meshconfig': {
                         'mesh': 'quad4',   #tria3
                         'elem': 'plane41', #plane31
                         'sizeelement': 5,
                         'meshmap': {'on': True,
                                     'edge': 'all',
                                    #  "numbernodes": 10,
                                     }
                         }
                     },
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            "FORCES": [f1, f2, f3, f4],
            "BOUNDCOND": [bc1, bc2],
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_03b',
                         'show': True,
                         'scale': 4,
                         'savepng': True,
                         'lines': True,
                         'plottags': {
                            # 'point': True,
                            'edge': True,
                            # 'surf': True
                            }
                        },
              }

preview_plot(previewset, modelinfo)