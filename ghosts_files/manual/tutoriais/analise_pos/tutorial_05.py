'''
myfempy Tutorial 05

Geração da malha com gmsh/ solidos

Necessario instalação prévia do gmsh (não nativo do myfempy)

'''

# Imports

from myfempy import ModelGen
from myfempy import Solver
from myfempy import PostProcess
from myfempy import postproc_plot
from myfempy import preview_plot

import sys

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


# #-----------------------------------------------------

# # 1/4 de placa com furo central (dia. 100 mm) estado plano (malha quad 4) 200 x 200 mm, espessura da placa de 20 mm

# points = [
#     [100, 0, 0],
#     [100, 100, 0],
#     [0, 100, 0],
# ]
         
# lines = [[1, 2],
#          [2, 3],
#          [3, 5],
#          [4, 1],
#          ]

# arcs = [[50, [0, 0, 0], ['0', 'Pi/2']]]

# plane = [[1, 2, 3, 4, 5]]

# meshdata = {"GMSH": {'filename': 'tutorial_05a',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'arc': arcs,
#                      'meshconfig': {
#                          'mesh': 'tetr4',   #tetr4 hexa8
#                          'elem': 'solid41', #solid41 solid81
#                          'sizeelement': 5,
#                          'extrude': 20,
#                          'meshmap': {'on': True,
#                                      'edge': 'all', #'all'
#                                     #  "numbernodes": 10,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_05a',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {'point': True,
#                                       'edge': True}},
#               }

# preview_plot(previewset, modelinfo)


#------------------------------------------


# # Placa completa multiplos furos (dia. 10, 10, 30 mm) estado plano (malha quad 4) 200 x 200 mm, espessura da placa de 20 mm

# points = [
#     [0, 0, 0],
#     [200, 0, 0],
#     [200, 200, 0],
#     [0, 200, 0],
# ]
         

# lines = [[1, 2], # line 1
#          [2, 3], # line 2
#          [3, 4], # line 3
#          [4, 1], # line 4
#          ]

# arcs = [[30, [150, 100, 0], ['0', '2*Pi']], # line 5
#         [10, [30, 150, 0], ['0', '2*Pi']],  # line 6
#         [10, [30, 50, 0], ['0', '2*Pi']],   # line 7
#         ]
         
# plane = [[1, 2, 3, 4],
#         [-5],
#         [-6],
#         [-7],
#         # [5],
#         # [6],
#         # [7]
#         ]

# meshdata = {"GMSH": {'filename': 'tutorial_05b',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'arc': arcs,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'tetr4',   #tetr4 hexa8
#                          'elem': 'solid41', #solid41 solid81
#                          'sizeelement': 5,
#                          'extrude': 20,
#                          'meshmap': {'on': True,
#                                      'edge': 'all', #'all'
#                                     #  "numbernodes": 10,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }


# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_05b',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {
#                             # 'point': True,
#                             # 'edge': True,
#                             'surf': True
#                                       }},
#               }

# preview_plot(previewset, modelinfo)




# # ------------------------------------------------------

# # Placa completa furo quadrado central (dia. 50 mm) estado plano (malha quad 4) 200 x 200 mm, espessura da placa de 20 mm

# points = [
#     [0, 0, 0],
#     [200, 0, 0],
#     [200, 100, 0],
#     [0, 100, 0],
#     [75, 25, 0],
#     [125, 25, 0],
#     [125, 75, 0],
#     [75, 75, 0],
#     ]
         

# lines = [[1, 2], # line 1
#          [2, 3], # line 2
#          [3, 4], # line 3
#          [4, 1], # line 4
#          [5, 6], # line 5
#          [6, 7], # line 6
#          [7, 8], # line 7
#          [8, 5], # line 8
#          ]

         
# plane = [[1, 2, 3, 4],
#         [-5, -6, -7, -8],
#         # [5, 6, 7, 8]
#         ]

# meshdata = {"GMSH": {'filename': 'tutorial_05c',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'hexa8',   #tetr4 hexa8
#                          'elem': 'solid81', #solid41 solid81
#                          'sizeelement': 5,
#                          'extrude': 20,
#                          'meshmap': {'on': True,
#                                      'edge': 'all', #'all'
#                                     #  "numbernodes": 10,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }


# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_05c',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {
#                             'point': True,
#                             # 'edge': True,
#                             # 'surf': True
#                             }},
#               }

# preview_plot(previewset, modelinfo)



# ------------------------------------------------------


# # Duas Placa composta estado plano (malha quad 4) 200 x 200 mm, espessura da placa de 20 mm

# points = [
#     [0, 0, 0],     # point 1
#     [100, 0, 0],   # point 2
#     [100, 100, 0], # point 3
#     [0, 100, 0],   # point 4
#     [200, 0, 0],   # point 5
#     [200, 100, 0], # point 6
#     ]
         

# lines = [[1, 2], # line 1
#          [2, 3], # line 2
#          [3, 4], # line 3
#          [4, 1], # line 4
#          [2, 5], # line 5
#          [5, 6], # line 6
#          [6, 3], # line 7
#          ]

         
# plane = [[1, 2, 3, 4], # plane 1
#         [5, 6, 7, 2],  # plane 2
#         ]

# meshdata = {"GMSH": {'filename': 'tutorial_05d',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'tetr4',   #tetr4 hexa8
#                          'elem': 'solid41', #solid41 solid81
#                          'sizeelement': 5,
#                          'extrude': 20,
#                          'meshmap': {'on': True,
#                                      'edge': [2, 6], #'all'
#                                      "numbernodes": 10,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }


# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_05d',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {
#                             # 'point': True,
#                             'edge': True,
#                             # 'surf': True
#                             }},
#               }

# preview_plot(previewset, modelinfo)


# # ------------------------------------------------------

# # Duas Placa composta com furo em 0,0,0 (dia. 50 mm) estado plano (malha quad 4) 200 x 200 mm, espessura da placa de 20 mm


f1 = {
    'DEF': 'forcesurf',
    'DOF': 'fx',
    'DIR': 'surf',
    'TAG': 10,
    'VAL': [100.0],
    }

bc1 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'surf',
    'TAG': 5,
    }

bc2 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'surf',
    'TAG': 11,
    }


bc3 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'surf',
    'TAG': 3,
    }


bc4 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'surf',
    'TAG': 9,
    }


bc5 = {
    'DEF': 'fixed',
    'DOF': 'ux',
    'DIR': 'surf',
    'TAG': 6,
    }


mat1 = {
    "NAME": "material1",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'solid'
    }


mat2 = {
    "NAME": "material2",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'solid'
    }


points = [
    [100, 0, 0],   # point 1
    [100, 100, 0], # point 2
    [0, 100, 0],   # point 3
    [200, 0, 0],   # point 4
    [200, 100, 0], # point 5
    ]
         

lines = [
    [6, 1], # line 1
    [1, 2], # line 2
    [2, 3], # line 3
    [3, 7], # line 4
    [1, 4], # line 5
    [4, 5], # line 6
    [5, 2], # line 7
    ]   

arcs = [[50, [0, 0, 0], ['0', 'Pi/2']], # line 8
        ]
         
plane = [[1, 2, 3, 4, 8], # plane 1
        [5, 6, 7, 2],  # plane 2
        ]

meshdata = {"GMSH": {'filename': 'tutorial_05e',
                     'pointlist': points,
                     'linelist': lines,
                     'arc': arcs,
                     'planelist': plane,
                     'meshconfig': {
                         'mesh': 'tetr4',   #tetr4 hexa8
                         'elem': 'solid41', #solid41 solid81
                         'sizeelement': 5,
                         'extrude': 20,
                        #  'meshmap': {'on': True,
                        #              'edge': 'all',
                        #              "numbernodes": 20,
                        #              }
                         }
                     },
            "PROPMAT": [mat1,mat2],
            "PROPGEO": [geo,geo],
            "FORCES": [f1],
            "BOUNDCOND": [bc1, bc2, bc3, bc4, bc5],
            # "QUADRATURE":{'meth':'gaussian','npp':4},
            "DOMAIN":'structural'
            }


modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_05e',
                         'show': True, 'scale': 10,
                         'savepng': True,
                         'lines': True,
                         'plottags': {
                            # 'point': True,
                            # 'edge': True,
                            'surf': True
                            }},
              }

preview_plot(previewset, modelinfo)


#-------------------------------- SOLVER -------------------------------------#
solverset = {"SOLVER": 'SLIPRE', #SLI
             'TOL': 1E-8,
             "STEPSET": {'type': 'table',  # mode, freq, time ...
                         'start': 0,
                         'end': 1,
                         'step': 1},
             }

solution = Solver.get_static_solve(solverset, modelinfo)

#----------------------------- POST-PROCESS ----------------------------------#
postprocset = {"SOLUTION": solution,
                "COMPUTER": {'elasticity': {'displ': True, 'stress': True, 'average': True}},
                "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'tutorial_05e_sim', 'savepng': True},
            }

postporc_result = PostProcess(modelinfo).compute(postprocset)

#----------------------------- VIEW SOLUTION ---------------------------------#
postproc_plot(postprocset, postporc_result, modelinfo)
