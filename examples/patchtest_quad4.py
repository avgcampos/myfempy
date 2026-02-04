from myfempy import newAnalysis
from myfempy import SteadyStateLinear

fea = newAnalysis(SteadyStateLinear)

mat = {
    'NAME': 'material',
    'EXX': 1000,
    'VXX': 0.3
}

geo = {
    'NAME':'geo',
    'THICKN': 0.1,
}

nodes = [
    [1,     0.00,   0.00,   0.00],
    [2,     2.00,   0.00,   0.00],
    [3,     2.00,   3.00,   0.00],
    [4,     0.00,   2.00,   0.00],
    [5,     0.40,   0.40,   0.00],
    [6,     1.40,   0.60,   0.00],
    [7,     1.50,   2.00,   0.00],
    [8,     0.30,   1.60,   0.00],    
        ]

conec = [
    [1, 1, 1, 1, 2, 6, 5],
    [2, 1, 1, 2, 3, 7, 6],
    [3, 1, 1, 3, 4, 8, 7],
    [4, 1, 1, 1, 5, 8, 4],
    [5, 1, 1, 5, 6, 7, 8],
        ]

modeldata = {
    'MESH':{
        'TYPE':'manual',
        'COORD':nodes,
        'INCI':conec
    },

    'ELEMENT':{
        'TYPE':'structplane',
        'SHAPE':'quad4',
        'INTGAUSS':2,
    },

    'MATERIAL':{
        'MAT':'planestress',
        'TYPE':'isotropic',
        'PROPMAT':[mat]
    },

    'GEOMETRY':{
        'GEO':'thickness',
        'PROPGEO':[geo]
    }
}

fea.Model(modeldata)


previewset = {
    'RENDER':{
        'filename':'patchtest',
        'show':True,
        'scale':5,
        'savepng':True,
        'lines':True
    }
}

bcfix = {
    'TYPE':'fixed',
    'DOF':'full',
    'DIR':'node',
    'LOC':{'x':0, 'y':0, 'z':0},
}

def displ_patch_test(x,y):
    ux = 0.002*x
    uy = -0.0006*y
    return ux, uy

bcnh_node1_ux = {
    'TYPE':'displ',
    'DOF':'ux',
    'DIR':'node',
    'LOC':{'x':0.0, 'y':0, 'z':0},
    'VAL':[displ_patch_test(0.0,0.00)[0]]
} 

bcnh_node1_uy = {
    'TYPE':'displ',
    'DOF':'uy',
    'DIR':'node',
    'LOC':{'x':0.0, 'y':0, 'z':0},
    'VAL':[displ_patch_test(0.0,0.00)[1]]
} 

bcnh_node2_ux = {
    'TYPE':'displ',
    'DOF':'ux',
    'DIR':'node',
    'LOC':{'x':2.0, 'y':0, 'z':0},
    'VAL':[displ_patch_test(2.0,0.00)[0]]
} 

bcnh_node2_uy = {
    'TYPE':'displ',
    'DOF':'uy',
    'DIR':'node',
    'LOC':{'x':2.0, 'y':0, 'z':0},
    'VAL':[displ_patch_test(2.0,0.00)[1]]
} 

bcnh_node3_ux = {
    'TYPE':'displ',
    'DOF':'ux',
    'DIR':'node',
    'LOC':{'x':2.0, 'y':3.0,'z':0},
    'VAL':[displ_patch_test(2.0,3.0)[0]]
} 

bcnh_node3_uy = {
    'TYPE':'displ',
    'DOF':'uy',
    'DIR':'node',
    'LOC':{'x':2.0, 'y':3.0, 'z':0},
    'VAL':[displ_patch_test(2.0,3.0)[1]]
} 


bcnh_node4_ux = {
    'TYPE':'displ',
    'DOF':'ux',
    'DIR':'node',
    'LOC':{'x':0, 'y':2.0,'z':0},
    'VAL':[displ_patch_test(0,2.0)[0]]
} 

bcnh_node4_uy = {
    'TYPE':'displ',
    'DOF':'uy',
    'DIR':'node',
    'LOC':{'x':0, 'y':2.0, 'z':0},
    'VAL':[displ_patch_test(0,2.0)[1]]
} 

physicdata = {
    'PHYSIC':{
        'DOMAIN':'structural',
        'LOAD':[],
        'BOUNDCOND':[
                    # bcfix,
                    bcnh_node1_ux, bcnh_node1_uy,
                    bcnh_node2_ux, bcnh_node2_uy,
                    bcnh_node3_ux, bcnh_node3_uy,
                    bcnh_node4_ux, bcnh_node4_uy,]
    }
}

fea.Physic(physicdata)


fea.PreviewAnalysis(previewset)

solverset = {
    'STEPSET':{
        'type':'table',
        'start':0,
        'end':1,
        'step':1
    },
    'SYMM':False,
    'MP':False
}

solverdata = fea.Solve(solverset)

print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'PatchTest', 'savepng': True},
                "OUTPUT": {'log': True, 'get':{
                        'nelem': True,
                        'nnode': True,
                        'inci': True,
                        'coord':True,
                        'tabmat':True,
                        'tabgeo':True,
                        'boundcond_list':True,
                        'forces_list':True,
                    }
            }}
postprocdata = fea.PostProcess(postprocset)

print('STRESS_XX: ',postprocdata['STRESS_XX'])