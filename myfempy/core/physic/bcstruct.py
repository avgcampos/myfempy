from __future__ import annotations

import numpy as np

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import get_nodes_from_list

class BoundCondStruct(Structural):
    '''Structural Load Class <ConcreteClassService>'''

    def getBCApply(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 4))
        for bc_index in range(len(bclist)):
            if bclist[bc_index][0] == "fixed":
                bcl = bclist[bc_index]
                bcapp = BoundCondStruct.__BCFixed(modelinfo, bcl)
                boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
            elif bclist[bc_index][0] == "displ":
                bcl = bclist[bc_index]
                bcapp = BoundCondStruct.__BCDispl(modelinfo, bcl)
                boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
            else:
                pass
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply
               
    def __BCFixed(modelinfo, bclist):
        
        boncdnodeaply = np.zeros((1, 4))
       
        nodelist = bclist[2:]
        node_list_bc, dir_fc = get_nodes_from_list(nodelist, modelinfo['coord'], modelinfo['regions'])
        
        if bclist[1] == 'full':
            bcdof = 0
        else:
            bcdof = modelinfo['dofs']['d'][bclist[1]]
        
        for j in range(len(node_list_bc)):
            bcapp = np.array([[int(node_list_bc[j]), bcdof, 0.0, int(bclist[8])]])
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
                
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply
    
    def __BCDispl(modelinfo, bclist):
        
        boncdnodeaply = np.zeros((1, 4))
       
        nodelist = bclist[2:]
        node_list_bc, dir_fc = get_nodes_from_list(nodelist, modelinfo['coord'], modelinfo['regions'])
        
        bcdof = modelinfo['dofs']['d'][bclist[1]]
        
        for j in range(len(node_list_bc)):
            bcapp = np.array([[int(node_list_bc[j]), bcdof, float(bclist[7]), int(bclist[8])]])
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
                
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply