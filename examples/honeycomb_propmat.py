import sys
# setting path
sys.path.append('../myfempy')

import numpy as np
from myfempy import newAnalysis
from myfempy import HomogenizationPlane

from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(HomogenizationPlane)

E_solid = 71E3 # #200E3 #0.91
v = 0.33
r_solid = 2.77E-9 #7.850E-9 #1

mat1 = {
    "NAME": "solid",
    "VXY": v,
    "EXX": E_solid,       # N/mm^2 --> MPa
    "RHO": r_solid,       # kg/mm^2
    }

xmin = 1E-3                # [] Xmin
penal = 3                   # [] Penal Material Interpolation Factor
mat2 = {
    "NAME": "void",
    "VXY": v,
    "EXX": E_solid*xmin**penal,
    "RHO": r_solid*xmin**penal,
    }

geo = {
    "NAME": "espessura",
    "THICKN": 1,
    }


# # HONEYCOMB MODEL SET
a = 20 #3**(-3/4)           # Comprimento do lado do hexágono
t = 4 #a * (np.sqrt(3)/6)  # Espessura da parede

# DImensoes da celula periodica
Lax = 3 * a
Lay = np.sqrt(3) * a    # Altura exata para um hexágono equilátero

# Centro da célula
cx = Lax / 2
cy = Lay / 2

# Fatores de correção para espessura constante t em ângulos de 60°
# dx é o deslocamento horizontal nas quinas para manter espessura t
ty_bench = t / 2
dx = (t/2) * np.sqrt(3) 

points = [
    # --- CONTORNO EXTERNO (10 pontos) ---
    [0, cy + t/2, 0],                # 1: Início braço esquerdo
    [a/2, cy + t/2, 0],              # 2: Encontro braço/hexágono
    [a, Lay, 0],                     # 3: Quina superior esquerda
    [2*a, Lay, 0],                   # 4: Quina superior direita
    [2.5*a, cy + t/2, 0],            # 5: Encontro hexágono/braço direito
    [3*a, cy + t/2, 0],              # 6: Fim braço direito
    [3*a, cy - t/2, 0],              # 7: Fim braço direito baixo
    [2.5*a, cy - t/2, 0],            # 8: Encontro hexágono/braço direito baixo
    [2*a, 0, 0],                     # 9: Quina inferior direita
    [a, 0, 0],                       # 10: Quina inferior esquerda
    [a/2, cy - t/2, 0],              # 11: Encontro braço/hexágono esquerdo baixo
    [0, cy - t/2, 0],                # 12: Início braço esquerdo baixo

    # --- CONTORNO INTERNO (6 pontos - Hexágono Equilátero) ---
    [a/2 + dx, cy, 0],               # 13: Ponta interna esquerda
    [a + dx, Lay - ty_bench, 0],     # 14: Topo interno esquerdo
    [2*a - dx, Lay - ty_bench, 0],   # 15: Topo interno direito
    [2.5*a - dx, cy, 0],             # 16: Ponta interna direita
    [2*a - dx, ty_bench, 0],         # 17: Base interna direita
    [a + dx, ty_bench, 0],           # 18: Base interna esquerda

    # --- CONTORNO FECHAMENTO RETANGULO EXTERNO (4 pontos) ---
    [0, 0, 0],                       # 19: Quina inferior esquerda
    [3*a, 0, 0],                     # 20: Quina inferior direita
    [3*a, Lay, 0],                   # 21: Quina superior direita
    [0, Lay, 0],                     # 22: Quina superior esquerda
    
]

lines = [
    # Externo
    [1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,11], [11,12], [12,1],
    # Interno
    [13,14], [14,15], [15,16], [16,17], [17,18], [18,13],
    # Retagulo Externo
    [12,19], [19,10], [9, 20], [20, 7], [6, 21], [21, 4], [3, 22], [22, 1],
]

plane = [
    # honeycomb
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
   
    # furo hexagonal
    [-13, -14, -15, -16, -17, -18],

    # preenchimento com material
    [13, 14, 15, 16, 17, 18],

    # fechamento retangulos
    [1, 26, 25, 2],
    [20, 19, 11, 10],
    [21, 8, 7, 22],
    [4, 24, 23, 5],
]

esize = 1

modeldata = {

   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'honeycomb',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'quad4',   #quad4 tria3
            'sizeelement': 2*esize,
            'meshmap': {'on': True,
                        'edge': 'all', #[[1,3], [2,4]], #'all'
                        # "numbernodes": [20, 10],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        # 'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat1, mat2, mat2, mat2, mat2, mat2],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo, geo, geo, geo, geo, geo],
    },
}
fea.Model(modeldata)

# inci = fea.getInci()
# tabmat = fea.getTabmat()
# print(inci[:,2])

# print('tabmat: ', tabmat)

bc_X0_XX = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'STEP': 1
    }

bc_X1_XX = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgex',
    'LOC': {'x': Lax, 'y': 999, 'z': 0},
    'STEP': 1
    }

bc_Y0_XX = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'STEP': 1
    }

bc_Y1_XX = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': Lay, 'z': 0},
    'STEP': 1
    }

bc_X0_YY = bc_X0_XX.copy()
bc_X0_YY['STEP'] = 2

bc_X1_YY = bc_X1_XX.copy()
bc_X1_YY['STEP'] = 2

bc_Y0_YY = bc_Y0_XX.copy()
bc_Y0_YY['STEP'] = 2

bc_Y1_YY = bc_Y1_XX.copy()
bc_Y1_YY['STEP'] = 2

bc_X0_XY = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'STEP': 3
    }

bc_X1_XY = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgex',
    'LOC': {'x': Lax, 'y': 999, 'z': 0},
    'STEP': 3
    }

bc_Y0_XY = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'STEP': 3
    }

bc_Y1_XY = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': Lay, 'z': 0},
    'STEP': 3
    }

strainzero = {
    'TYPE': 'strainzero',
    'VAL': [[1,0,0], [0,1,0], [0,0,1]],    # mm/s^2 
    'DOF': 'none',
    'DIR':'none',
    'MESHNODE': 0,
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [strainzero],
               "BOUNDCOND": [
                    bc_X0_XX, bc_X1_XX, bc_Y0_XX, bc_Y1_XX,
                    bc_X0_YY, bc_X1_YY, bc_Y0_YY, bc_Y1_YY,
                    bc_X0_XY, bc_X1_XY, bc_Y0_XY, bc_Y1_XY
                             ],
    },
}
fea.Physic(physicdata)

# loadaply = fea.getLoadApply()
# print('forca total [N]', np.sum(loadaply[:,2]))

# sys.exit()

# previewset = {'RENDER': {'filename': 'cell_preview', 'show': True, 'scale': 2, 'savepng': True, 'lines': False,
#                          'plottags': {'point': True}
#                          },
#             #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
#               }
# fea.PreviewAnalysis(previewset)

# sys.exit()
# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
             'RHOH': True,
             }
solverdata = fea.Solve(solverset)

#---------------------------- POST PROCESS  -----------------------------------#
postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'micro_honeycomb', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

# CH = solverdata["solution"]['CH']
# print('Homoge. Elastic Tensor\n', CH)

# RH = solverdata["solution"]['RHOH']
# print('Homoge. Density\n', RH)

# # Inversão para matriz de flexibilidade
# S = np.linalg.inv(CH)

# E_x = 1 / S[0,0]
# E_y = 1 / S[1,1]
# nu_xy = -S[0,1] / S[0,0]
# nu_yx = -S[0,1] / S[1,1]
# G_xy = 1 / S[2,2]

# print("E_x:", E_x)
# print("E_y:", E_y)
# print("v_xy:", nu_xy)
# print("v_yx:", nu_yx)
# print("G_xy:", G_xy)

# # Verificação de isotropia
# tol = 1e-2
# if abs(E_x - E_y) < tol and abs(nu_xy - nu_yx) < tol:
#     print("Material isotrópico")
# elif abs(nu_xy - nu_yx) < tol:
#     print("Material ortotrópico")
# else:
#     print("Material anisotrópico")





# # 2. Calcular a matriz de flexibilidade S (S = inv(C))
# S = np.linalg.inv(CH)

# # 3. Extrair as propriedades de engenharia
# # No estado plano de tensões:
# # S11 = 1/E1, S22 = 1/E2, S33 = 1/G12
# # S12 = -nu12/E1 => nu12 = -S12 * E1
# # S21 = -nu21/E2 => nu21 = -S21 * E2

# E1 = 1 / S[0, 0]
# E2 = 1 / S[1, 1]
# G12 = 1 / S[2, 2]

# nu12 = -S[0, 1] * E1
# # nu21 = -S[1, 0] * E2

# # 4. Exibir resultados
# print(f"--- Propriedades Calculadas ---")
# print(f"Módulo de Young E1:  {E1:10.2f}")
# # print(f"Módulo de Young E2:  {E2:10.2f}")
# print(f"Módulo de Cisalhamento G12: {G12:10.2f}")
# print(f"Poisson nu12:        {nu12:10.4f}")
# # print(f"Poisson nu21:        {nu21:10.4f}")

# # # Verificação de simetria (Consistência termodinâmica)
# # print(f"\nVerificação de simetria (nu12/E1 == nu21/E2):")
# # print(f"{nu12/E1:.2e} == {nu21/E2:.2e}")
