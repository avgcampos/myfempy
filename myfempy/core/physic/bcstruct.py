from __future__ import annotations

import numpy as np

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import get_nodes_from_list

class BoundCondStruct(Structural):
    '''Structural Load Class <ConcreteClassService>'''

    def getBCApply(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 2))
        for bc_index in range(len(bclist)):
            if bclist[bc_index][0] == "fixed":
                bcl = bclist[bc_index]
                bcapp = BoundCondStruct.getBCFixed(modelinfo, bcl)
                boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
            else:
                pass

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply
        
    def setNodes(bclist, coord, regions):
        return get_nodes_from_list(bclist, coord, regions)
    
    def setBCDof(bclist, node_list_bc):
        return Structural.setBCDof(bclist, node_list_bc)
               
    def getBCFixed(modelinfo, bclist):
        
        boncdnodeaply = np.zeros((1, 2))
       
        nodelist = bclist[2:]
        node_list_bc, dir_fc = BoundCondStruct.setNodes(nodelist, modelinfo['coord'], modelinfo['regions'])
        
        if bclist[1] == 'all':
            bcdof = 0
        else:
            bcdof = modelinfo['dofs']['d'][bclist[1]]
        for j in range(len(node_list_bc)):
            # bcdof = BoundCondStruct.setBCDof(bclist, node_list_bc[j])
            bcapp = np.array([[bcdof, int(node_list_bc[j])]])
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
                
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply