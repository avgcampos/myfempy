from __future__ import annotations

import numpy as np

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import get_nodes_from_list, search_edgex, search_edgey, search_edgez, search_surfxy, search_surfyz, search_surfzx, search_nodexyz

class BoundCondStruct(Structural):
    '''Structural Load Class <ConcreteClassService>'''

    def getBCApply(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 2))
        for bc_index in range(len(bclist)):
            if bclist[bc_index][0] == "fixed":
                bcl = bclist[bc_index]
                bcapp = BoundCondStruct.__BCFixed(modelinfo, bcl)
                boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
            else:
                pass

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply
        
    def setNodes(nodelist, coord, regions):
      return get_nodes_from_list(nodelist, coord, regions)
    
    def setBCDof(bclist, node_list_bc):
        if bclist[1] == "ux":
            bcdof = np.array([[1, int(node_list_bc)]])
        elif bclist[1] == "uy":
            bcdof = np.array([[2, int(node_list_bc)]])
        elif bclist[1] == "uz":
            bcdof = np.array([[3, int(node_list_bc)]])
        elif bclist[1] == "rx":
            bcdof = np.array([[4, int(node_list_bc)]])
        elif bclist[1] == "ry":
            bcdof = np.array([[5, int(node_list_bc)]])
        elif bclist[1] == "rz":
            bcdof = np.array([[6, int(node_list_bc)]])
        elif bclist[1] == "all":
            bcdof = np.array([[0, int(node_list_bc)]])
        else:
            print("input erro: bc_opt_dir don't defined")  
        return bcdof
               
    def __BCFixed(modelinfo, bclist):
        
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