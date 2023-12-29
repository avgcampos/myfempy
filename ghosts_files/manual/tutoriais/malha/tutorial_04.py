'''
myfempy Tutorial 04

Geração da malha com gmsh/ criação de curvas e furos e áreas compostas

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

# # 1/4 de placa com furo central (dia. 100 mm) estado plano (malha quad 4) 200 x 200 mm

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

# meshdata = {"GMSH": {'filename': 'tutorial_04a',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'arc': arcs,
#                      'meshconfig': {
#                          'mesh': 'quad4',   #quad4
#                          'elem': 'plane41', #plane41
#                          'sizeelement': 5,
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

# previewset = {'RENDER': {'filename': 'tutorial_04a',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {'point': True,
#                                       'edge': True}},
#               }

# preview_plot(previewset, modelinfo)


# sys.exit()


#------------------------------------------


# # Placa completa multiplos furos (dia. 10, 10, 30 mm) estado plano (malha quad 4) 200 x 200 mm

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
#         [5],
#         [6],
#         [7]
#         ]

# meshdata = {"GMSH": {'filename': 'tutorial_04b',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'arc': arcs,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'quad4',   #quad4
#                          'elem': 'plane41', #plane41
#                          'sizeelement': 10,
#                         #  'meshmap': {'on': True,
#                         #              'edge': [5], #'all'
#                         #              "numbernodes": 60,
#                         #              }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }


# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_04b',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {
#                             # 'point': True,
#                             # 'edge': True,
#                             'surf': True
#                                       }},
#             'LABELS': {'show': True, 'lines': True, 'scale': 0.2},
#               }

# preview_plot(previewset, modelinfo)




# # ------------------------------------------------------

# # Placa completa furo quadrado central (dia. 50 mm) estado plano (malha quad 4) 200 x 200 mm

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

# meshdata = {"GMSH": {'filename': 'tutorial_04c',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'quad4',   #quad4
#                          'elem': 'plane41', #plane41
#                          'sizeelement': 2,
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

# previewset = {'RENDER': {'filename': 'tutorial_04c',
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


# # Duas Placa composta estado plano (malha quad 4) 200 x 200 mm

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

# meshdata = {"GMSH": {'filename': 'tutorial_04d',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'quad4',   #quad4
#                          'elem': 'plane41', #plane41
#                          'sizeelement': 5,
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

# previewset = {'RENDER': {'filename': 'tutorial_04d',
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


# ------------------------------------------------------

# # Duas Placa composta com furo em 0,0,0 (dia. 50 mm) estado plano (malha quad 4) 200 x 200 mm

# points = [
#     [100, 0, 0],   # point 1
#     [100, 100, 0], # point 2
#     [0, 100, 0],   # point 3
#     [200, 0, 0],   # point 4
#     [200, 100, 0], # point 5
#     ]
         

# lines = [
#     [6, 1], # line 1
#     [1, 2], # line 2
#     [2, 3], # line 3
#     [3, 7], # line 4
#     [1, 4], # line 5
#     [4, 5], # line 6
#     [5, 2], # line 7
#     ]   

# arcs = [[50, [0, 0, 0], ['0', 'Pi/2']], # line 8
#         ]
         
# plane = [[1, 2, 3, 4, 8], # plane 1
#         [5, 6, 7, 2],  # plane 2
#         ]

# meshdata = {"GMSH": {'filename': 'tutorial_04e',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'arc': arcs,
#                      'planelist': plane,
#                      'meshconfig': {
#                          'mesh': 'tria3',   #quad4
#                          'elem': 'plane31', #plane41
#                          'sizeelement': 10,
#                          'meshmap': {'on': True,
#                                      'edge': [8], #'all'
#                                      "numbernodes": 20,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }


# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_04e',
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

# -----------------------------------

# # disco (dia. 200 mm) com furo central (dia. 100 mm) estado plano (malha quad 4)

# points = [
#     # [100, 0, 0],
#     # [100, 100, 0],
#     # [0, 100, 0],
# ]
         
# lines = [
#     # [1, 2],
#     # [2, 3],
#     # [3, 5],
#     # [4, 1],
#          ]

# arcs = [[200, [0, 0, 0], ['0', '2*Pi']],[100, [0, 0, 0], ['0', '2*Pi']]]

# plane = [[1],[-2]]

# meshdata = {"GMSH": {'filename': 'tutorial_04f',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'arc': arcs,
#                      'meshconfig': {
#                          'mesh': 'tria3',   #quad4
#                          'elem': 'plane31', #plane41
#                          'sizeelement': 40,
#                          'meshmap': {'on': True,
#                                      'edge': 'all', #'all'
#                                     #  "numbernodes": 80,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_04f',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {
#                             'edge': True,
#                             # 'surf': True
#                             }},
#             # 'LABELS': {'show': True, 'lines': True, 'scale': 10},
#               }

# preview_plot(previewset, modelinfo)


# -----------------------------------


# # 1/4 de placa com bara no furo central (dia. 100 mm) estado plano (malha quad 4) 200 x 200 mm

# points = [
#     [100, 0, 0],
#     [100, 100, 0],
#     [0, 100, 0],
#     [0, 0, 0],
# ]
         
# lines = [[1, 2],
#          [2, 3],
#          [3, 6],
#          [5, 1],
#          [4, 5],
#          [6, 4],
#          ]

# arcs = [[50, [0, 0, 0], ['0', 'Pi/2']]]

# plane = [[1, 2, 3, 4, 7],[5, 6, 7]]

# meshdata = {"GMSH": {'filename': 'tutorial_04g',
#                      'pointlist': points,
#                      'linelist': lines,
#                      'planelist': plane,
#                      'arc': arcs,
#                      'meshconfig': {
#                          'mesh': 'quad4',   #quad4
#                          'elem': 'plane41', #plane41
#                          'sizeelement': 10,
#                          'meshmap': {'on': True,
#                                      'edge': [7], #'all'
#                                     #  "numbernodes": 10,
#                                      }
#                          }
#                      },
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_04g',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {
#                                 # 'point': True,
#                                 # 'edge': True,
#                                 'surf': True,
#                                 }},
#               }

# preview_plot(previewset, modelinfo)


# #-----------------------------------------------------

# 1/4 de placa com furo central (dia. 100 mm) estado plano (malha quad 4) 200 x 200 mm malha mapeada

points = [
    [100, 0, 0],
    [100, 100, 0],
    [0, 100, 0],
]
         
lines = [[1, 2], # 1
         [2, 3], # 2
         [3, 7], # 3
         [4, 1], # 4
         [2, 5], # 5
         ]

arcs = [[50, [0, 0, 0], ['0', 'Pi/4']],    # 6
        [50, [0, 0, 0], ['Pi/4', 'Pi/2']]] # 7

plane = [[1, 5, 6, 4],
         [5, 2, 3, 7]]

meshdata = {"GMSH": {'filename': 'tutorial_04h',
                     'pointlist': points,
                     'linelist': lines,
                     'planelist': plane,
                     'arc': arcs,
                     'meshconfig': {
                         'mesh': 'quad4',   #quad4
                         'elem': 'plane41', #plane41
                         'sizeelement': 10,
                         'meshmap': {'on': True,
                                     'edge': 'all', #'all'
                                     "numbernodes": 10,
                                     }
                         }
                     },
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_04h',
                         'show': True, 'scale': 10,
                         'savepng': True,
                         'lines': True,
                         'plottags': {
                            'point': True,
                            # 'edge': True
                            }},
            # 'LABELS': {'show': True, 'lines': True, 'scale': 0.5},   
              }

preview_plot(previewset, modelinfo)