import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import PhononicCrystalPlaneBCPeriodic

import numpy as np
import matplotlib.pyplot as plt
from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(PhononicCrystalPlaneBCPeriodic)

mat_hard = {
    "NAME": "hard",
    "VXY": 0.34,
    "EXX": 20E3, #4.3E3    # MPa 
    "RHO": 2.0E-6, # 1.142E-6
    }

mat_soft = {
    "NAME": "soft",
    "VXY": 0.34,
    "EXX": 4E3, # 70E3    # MPa
    "RHO": 1.0E-6, # 2.7E-6
    }

geo = {"NAME": "geo1",
       "THICKN": 1}

# MODEL SET 
Lc = 6  # mm 


modeldata = {

    "MESH": {
        'TYPE': 'legacy',
        'LX': Lc,
        'LY': Lc,
        'NX': 10,
        'NY': 10,
        },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        # 'INTGAUSS': 8,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat_hard],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

# inci = fea.getInci()
# coord = fea.getCoord()
# tabmat = fea.getTabmat()
# tabgeo = fea.getTabgeo()
# intgauss = fea.getIntGauss()

# element_number = 0
# ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# print(np.array2string(ke, separator=', '))
# print(np.allclose(ke, ke.T, rtol=1e-05, atol=1e-08))

# element_number = 0
# me = fea.getElemMassConsistentMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# print(np.array2string(me, separator=', '))
# print(np.allclose(me, me.T, rtol=1e-05, atol=1e-08))

# kg = fea.getGlobalMatrix(inci, coord, tabmat, tabgeo, intgauss, SYMM=True)
# Kg_fem = kg['stiffness'].todense()
# Mg_fem = kg['mass'].todense()

# print(np.allclose(Kg_fem, Kg_fem.T, rtol=1e-05, atol=1e-08))

# plt.figure(2)
# plt.spy(kg, markersize=4)
# plt.show()

# import scipy as sp
# omega, _ = sp.linalg.eig(Kg_fem, Mg_fem)
# omega = np.sort(np.real(omega))  # Ordena os valores reais
# freqs = np.sqrt(omega[:10]) / (2 * np.pi)  # Convertendo para Hz

# print(freqs)

# sys.exit()

# pc_left = {
#     'TYPE': 'periplane', 
#     'DOF': 'left',
#     'DIR': 'line',
#     'TAG': 4,
#     }

# pc_right = {
#     'TYPE': 'periplane',  
#     'DOF': 'right',
#     'DIR': 'line',
#     'TAG': 2,
#     }

# pc_bottom = {
#     'TYPE': 'periplane',  
#     'DOF': 'bottom',
#     'DIR': 'line',
#     'TAG': 1,
#     }

# pc_top = {
#     'TYPE': 'periplane', 
#     'DOF': 'top',
#     'DIR': 'line',
#     'TAG': 3,
#     }

# pc_bottom_left = {
#     'TYPE': 'periplane',  
#     'DOF': 'bottom-left',
#     'DIR': 'point',
#     'TAG': 1,
#     }

# pc_bottom_right = {
#     'TYPE': 'periplane',  
#     'DOF': 'bottom-right',
#     'DIR': 'point',
#     'TAG': 2,
#     }

# pc_top_left = {
#     'TYPE': 'periplane',  
#     'DOF': 'top-left',
#     'DIR': 'point',
#     'TAG': 4,
#     }

# pc_top_right = {
#     'TYPE': 'periplane',  
#     'DOF': 'top-right',
#     'DIR': 'point',
#     'TAG': 3,
#     }


pc_left = {
    'TYPE': 'bloch', 
    'DOF': 'left',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }

pc_right = {
    'TYPE': 'bloch',  
    'DOF': 'right',
    'DIR': 'edgex',
    'LOC': {'x': Lc, 'y': 999, 'z': 0},
    }

pc_bottom = {
    'TYPE': 'bloch',  
    'DOF': 'bottom',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    }

pc_top = {
    'TYPE': 'bloch', 
    'DOF': 'top',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': Lc, 'z': 0},
    }

pc_bottom_left = {
    'TYPE': 'bloch',  
    'DOF': 'bottom-left',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 0, 'z': 0},
    }

pc_bottom_right = {
    'TYPE': 'bloch',  
    'DOF': 'bottom-right',
    'DIR': 'node',
    'LOC': {'x': Lc, 'y': 0, 'z': 0},
    }

pc_top_left = {
    'TYPE': 'bloch',  
    'DOF': 'top-left',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': Lc, 'z': 0},
    }

pc_top_right = {
    'TYPE': 'bloch',  
    'DOF': 'top-right',
    'DIR': 'node',
    'LOC': {'x': Lc, 'y': Lc, 'z': 0},
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [],
               "BOUNDCOND": [pc_left, pc_right, pc_bottom, pc_top, pc_bottom_left, pc_bottom_right, pc_top_left, pc_top_right],
    },
}
fea.Physic(physicdata)

# bcaply = fea.getBCApply()
# print(bcaply)

# freedof, fixedof, constdof = fea.getConstrains(bcaply)
# print('free',freedof)
# print('fixed',fixedof)
# print('const',constdof)
# print('forca total [N]', np.sum(loadaply[:,2]))


previewset = {'RENDER': {'filename': 'phc_preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)


# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'modes',
                        'start': 0,
                        'end': 8,
                        'step': 1},
             'SYMM':True,
             'IBZ': np.array([[0.0, 0.0],
                              [np.pi, 0.0],
                              [np.pi, np.pi],
                              [0.0, 0.0]]),
             }
solverdata = fea.Solve(solverset)

# PLOT DISPERCAO DE ONDAS
freqs_evo = solverdata['solution']['FREQ']
tot_steps = solverdata['solution']['IBZRANGE']

# Suponha que kx_values, Lc, pi, freqs e neig já tenham sido definidos
plt.figure()
plt.grid(True)
for idx in range(freqs_evo.shape[0]):
    plt.plot(np.arange(0, tot_steps[-1] + 1), freqs_evo[idx, :]/1e3, '.', mfc='none')
plt.xlabel(r"$\Re(kx)$", fontsize=9)
plt.ylabel('Frequência (kHz)', fontsize=9)
plt.xticks(tot_steps, ['O','A','B','O'], fontsize=9)
plt.yticks(fontsize=9)
plt.show()