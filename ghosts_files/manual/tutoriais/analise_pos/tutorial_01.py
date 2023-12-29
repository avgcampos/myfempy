'''
myfempy Tutorial 01

Geração da malha manual

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

# # Retangulo estado plano (malha manual quad 4) 100 x 50 mm com forcas e bc. nodais

# f1 = {
#     'DEF': 'forcenode',
#     'DOF': 'fy',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     'VAL': [-1000.0],
#     }

# f2 = {
#     'DEF': 'forcenode',
#     'DOF': 'fy',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 50, 'z': 0},
#     'VAL': [-1000.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 0, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 50, 'z': 0},
#     }


# elementos = [[1,"plane41","material","geometria",[1, 2, 3, 4]]]


# coordenadas = [[1, 0, 0, 0],
#                [2, 100, 0, 0],
#                [3, 100, 50, 0],
#                [4, 0, 50, 0]]


# meshdata = {"ADD121": elementos,
#             "NODELIST": coordenadas,
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1,f2],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE":{'meth':'gaussian','npp':4},
#             "DOMAIN":'structural'
#             }

# modelinfo = ModelGen.get_model(meshdata)

# # previewset = {'RENDER': {'filename': 'tutorial_01a', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
# #               'LABELS': {'show': True, 'lines': True, 'scale': 10},
# #               }

# # preview_plot(previewset, modelinfo)

# #-------------------------------- SOLVER -------------------------------------#
# solverset = {"SOLVER": 'SLD', #SLI
#              'TOL': 1E-8,
#              "STEPSET": {'type': 'table',  # mode, freq, time ...
#                          'start': 0,
#                          'end': 1,
#                          'step': 1},
#              #"TRACKER": {'show': False, 'result2plot': 'displ', 'max': []}
#              }

# solution = Solver.get_static_solve(solverset, modelinfo)

# #----------------------------- POST-PROCESS ----------------------------------#
# postprocset = {"SOLUTION": solution,
#                 "COMPUTER": {'elasticity': {'displ': True, 'stress': True, 'average': True}},
#                     # 'step':2
#                 "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'tutorial_01a_sim', 'savepng': True},
#                 # "TRACKER": {'show': True, 'result2plot':'stress', 'point': {'x':6,'y':1.5,'z':0}}
#             }

# postporc_result = PostProcess(modelinfo).compute(postprocset)

# #----------------------------- VIEW SOLUTION ---------------------------------#
# postproc_plot(postprocset, postporc_result, modelinfo)



#----------------------------------------------

# # Retangulo estado plano (malha legacy quad 4 10x5) 100 x 50 mm com forcas e bc distribuida

# f1 = {
#     'DEF': 'forceedge',
#     'DOF': 'fy',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 50, 'z': 0},
#     'VAL': [-1000.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     }

# meshdata = {"LEGACY": {'lx': 100, 'ly': 50, 'mesh': 'quad4', 'elem': 'plane41', 'nx': 10, 'ny': 5},
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
#             "DOMAIN":'structural'
#             }

# modelinfo = ModelGen.get_model(meshdata)

# # previewset = {'RENDER': {'filename': 'tutorial_01c', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
# #               'LABELS': {'show': True, 'lines': True, 'scale': 10},
# #               }

# # preview_plot(previewset, modelinfo)


# #-------------------------------- SOLVER -------------------------------------#
# solverset = {"SOLVER": 'SLD', #SLI
#              'TOL': 1E-8,
#              "STEPSET": {'type': 'table',  # mode, freq, time ...
#                          'start': 0,
#                          'end': 1,
#                          'step': 1},
#              #"TRACKER": {'show': False, 'result2plot': 'displ', 'max': []}
#              }

# solution = Solver.get_static_solve(solverset, modelinfo)

# #----------------------------- POST-PROCESS ----------------------------------#
# postprocset = {"SOLUTION": solution,
#                 "COMPUTER": {'elasticity': {'displ': True, 'stress': True, 'average': True}},
#                     # 'step':2
#                 "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'tutorial_01b_sim', 'savepng': True},
#                 # "TRACKER": {'show': True, 'result2plot':'stress', 'point': {'x':6,'y':1.5,'z':0}}
#             }

# postporc_result = PostProcess(modelinfo).compute(postprocset)

# #----------------------------- VIEW SOLUTION ---------------------------------#
# postproc_plot(postprocset, postporc_result, modelinfo)






# --------------------------------



# # Retangulo estado plano (malha legacy quad 4 10x5) 100 x 50 mm com forcas e bc distribuida

# f1 = {
#     'DEF': 'forceedge',
#     'DOF': 'fy',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 50, 'z': 0},
#     'VAL': [0, -10, -50, -100],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     }

# meshdata = {"LEGACY": {'lx': 100, 'ly': 50, 'mesh': 'quad4', 'elem': 'plane41', 'nx': 10, 'ny': 5},
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
#             "DOMAIN":'structural'
#             }

# modelinfo = ModelGen.get_model(meshdata)

# # previewset = {'RENDER': {'filename': 'tutorial_01c', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
# #               'LABELS': {'show': True, 'lines': True, 'scale': 10},
# #               }

# # preview_plot(previewset, modelinfo)


# #-------------------------------- SOLVER -------------------------------------#
# solverset = {"SOLVER": 'SLI', #SLI
#              'TOL': 1E-8,
#              "STEPSET": {'type': 'table',  # mode, freq, time ...
#                          'start': 0,
#                          'end': 4,
#                          'step': 1},
#              }

# solution = Solver.get_static_solve(solverset, modelinfo)

# #----------------------------- POST-PROCESS ----------------------------------#
# postprocset = {"SOLUTION": solution,
#                 "COMPUTER": {'elasticity': {'displ': True, 'stress': True, 'average': True}},
#                 "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'tutorial_01c_sim', 'savepng': True},
#             }

# postporc_result = PostProcess(modelinfo).compute(postprocset)

# #----------------------------- VIEW SOLUTION ---------------------------------#
# postproc_plot(postprocset, postporc_result, modelinfo)









# ----------------------------








# # # Retangulo estado plano dois materiais (malha manual quad 4) 200 x 100 mm com forcas e bc. nodais

# mat2 = {
#     "NAME": "material2",
#     "VXX": 0.25,
#     "EXX": 100E6,
#     "MAT": 'isotropic',
#     "DEF": 'planestress'
#     }

# geo2 = {
#     "NAME": "geometria2",
#     "THICKN": 1.0
#     }

# f1 = {
#     'DEF': 'forcenode',
#     'DOF': 'fx',
#     'DIR': 'node',
#     'LOC': {'x': 200, 'y': 0, 'z': 0},
#     'VAL': [100.0],
#     }

# f2 = {
#     'DEF': 'forcenode',
#     'DOF': 'fx',
#     'DIR': 'node',
#     'LOC': {'x': 200, 'y': 100, 'z': 0},
#     'VAL': [100.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 0, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     }

# bc3 = {
#     'DEF': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'node',
#     'LOC': {'x': 200, 'y': 0, 'z': 0},
#     }

# bc4 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 100, 'z': 0},
#     }

# bc5 = {
#     'DEF': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 100, 'z': 0},
#     }

# bc6 = {
#     'DEF': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'node',
#     'LOC': {'x': 200, 'y': 100, 'z': 0},
#     }


# elementos = [[1,"plane41","material2","geometria",[1, 2, 3, 4]],
#              [2,"plane41","material","geometria2",[2, 5, 6, 3]]]


# coordenadas = [[1, 0, 0, 0],
#                [2, 100, 0, 0],
#                [3, 100, 100, 0],
#                [4, 0, 100, 0],
#                [5, 200, 0, 0],
#                [6, 200, 100, 0]]


# meshdata = {"ADD121": elementos,
#             "NODELIST": coordenadas,
#             "PROPMAT": [mat2,mat],
#             "PROPGEO": [geo,geo2],
#             "FORCES": [f1,f2],
#             "BOUNDCOND": [bc1,bc2,bc3,bc4,bc5,bc6],
#             "QUADRATURE":{'meth':'gaussian','npp':4},
#             "DOMAIN":'structural'
#             }

# modelinfo = ModelGen.get_model(meshdata)

# print(modelinfo['inci'])

# previewset = {'RENDER': {'filename': 'tutorial_01d', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
#               'LABELS': {'show': True, 'lines': True, 'scale': 10},
#               }

# preview_plot(previewset, modelinfo)



# #-------------------------------- SOLVER -------------------------------------#
# solverset = {"SOLVER": 'SLI', #SLI
#              'TOL': 1E-8,
#              "STEPSET": {'type': 'table',  # mode, freq, time ...
#                          'start': 0,
#                          'end': 1,
#                          'step': 1},
#              }

# solution = Solver.get_static_solve(solverset, modelinfo)

# #----------------------------- POST-PROCESS ----------------------------------#
# postprocset = {"SOLUTION": solution,
#                 "COMPUTER": {'elasticity': {'displ': True, 'stress': True, 'average': True}},
#                 "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'tutorial_01d_sim', 'savepng': True},
#             }

# postporc_result = PostProcess(modelinfo).compute(postprocset)

# #----------------------------- VIEW SOLUTION ---------------------------------#
# postproc_plot(postprocset, postporc_result, modelinfo)




# -----------------------------------

# Dois Solidos 100 x 100 x 100 mm diferentes materiais

mat1 = {
    "NAME": "material1",
    "VXX": 0.25,
    "EXX": 100E6,
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

geo2 = {
    "NAME": "geometria2",
    "THICKN": 0.0
    }


coordenadas = [[1, 0, 0, 0],
               [2, 100, 0, 0],
               [3, 100, 100, 0],
               [4, 0, 100, 0],
               [5, 0, 0, 100],
               [6, 100, 0, 100],
               [7, 100, 100, 100],
               [8, 0, 100, 100],
               [9, 200, 0, 0],
               [10, 200, 100, 0],
               [11, 200, 0, 100],
               [12, 200, 100, 100],
               ]

elementos = [[1,"solid81","material1","geometria2",[1, 2, 3, 4, 5, 6, 7, 8]],
             [2,"solid81","material1","geometria2",[2, 9, 10, 3, 6, 11, 12, 7]]
             ]


# f1 = {
#     'DEF': 'forcesurf', #forcesurf
#     'DOF': 'fx',
#     'DIR': 'surfyz',
#     'LOC': {'x': 200, 'y': 999, 'z': 999},
#     'VAL': [100.0],
#     }

f1 = {
    'DEF': 'forcenode', #forcesurf
    'DOF': 'fz',
    'DIR': 'node',
    'LOC': {'x': 200, 'y': 100, 'z': 100},
    'VAL': [100.0],
    }

f2 = {
    'DEF': 'forcenode', #forcesurf
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 200, 'y': 0, 'z': 100},
    'VAL': [-100.0],
    }

f3 = {
    'DEF': 'forcenode', #forcesurf
    'DOF': 'fz',
    'DIR': 'node',
    'LOC': {'x': 200, 'y': 0, 'z': 0},
    'VAL': [-100.0],
    }

f4 = {
    'DEF': 'forcenode', #forcesurf
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 200, 'y': 100, 'z': 0},
    'VAL': [100.0],
    }


bc1 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'surfzx',
    'LOC': {'x': 999, 'y': 100, 'z': 999},
    }

bc2 = {
    'DEF': 'fixed',
    'DOF': 'uy',
    'DIR': 'surfzx',
    'LOC': {'x': 999, 'y': 0, 'z': 999},
    }

bc3 = {
    'DEF': 'fixed',
    'DOF': 'uz',
    'DIR': 'surfzx',
    'LOC': {'x': 999, 'y': 100, 'z': 999},
    }

bc4 = {
    'DEF': 'fixed',
    'DOF': 'uz',
    'DIR': 'surfzx',
    'LOC': {'x': 999, 'y': 0, 'z': 999},
    }

bc5 = {
    'DEF': 'fixed',
    'DOF': 'all',
    'DIR': 'surfyz',
    'LOC': {'x': 0, 'y': 999, 'z': 999},
    }


meshdata = {"ADD121": elementos,
            "NODELIST": coordenadas,
            "PROPMAT": [mat1,mat2],
            "PROPGEO": [geo2,geo2],
            "FORCES": [f1,f2,f3,f4],
            "BOUNDCOND": [bc5],
            "QUADRATURE":{'meth':'gaussian','npp':8},
            "DOMAIN":'structural'
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_01e', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
              'LABELS': {'show': True, 'lines': True, 'scale': 10},
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
                "PLOTSET": {'show': True, 'data': {'displ': []}, 'filename': 'tutorial_01e_sim', 'savepng': True},
            }

postporc_result = PostProcess(modelinfo).compute(postprocset)

#----------------------------- VIEW SOLUTION ---------------------------------#
postproc_plot(postprocset, postporc_result, modelinfo)


