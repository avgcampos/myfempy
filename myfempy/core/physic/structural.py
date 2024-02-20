from abc import ABC, abstractmethod

import numpy as np

from myfempy.core.utilities import search_edgex, search_edgey, search_edgez, search_surfxy, search_surfyz, search_surfzx, search_nodexyz


class Structural(ABC):
    '''Load API Class <ClassService>'''
    
    @abstractmethod 
    def getForceApply():
        pass

    @abstractmethod
    def getForceNodeLoadApply():
        pass
    
    @abstractmethod
    def getForceEdgeLoadApply():
        pass
    
    @abstractmethod
    def getForceSurfLoadApply():
        pass

    @abstractmethod
    def getBCApply():
        pass    
 
    @abstractmethod
    def getBCFixed():
        pass

    def setLoadDof(forcedof, node_list_fc):
        if forcedof == "fx":
            fcdof = 1 * np.ones_like(node_list_fc)
        elif forcedof == "fy":
            fcdof = 2 * np.ones_like(node_list_fc)
        elif forcedof == "fz":
            fcdof = 3 * np.ones_like(node_list_fc)
        elif forcedof == "tx":
            fcdof = 4 * np.ones_like(node_list_fc)
        elif forcedof == "ty":
            fcdof = 5 * np.ones_like(node_list_fc)
        elif forcedof == "tz":
            fcdof = 6 * np.ones_like(node_list_fc)
        elif forcedof == "masspoint":
            fcdof = 15 * np.ones_like(node_list_fc)
        elif forcedof == "spring2ground":
            fcdof = 16 * np.ones_like(node_list_fc)
        elif forcedof == "damper2ground":
            fcdof = 17 * np.ones_like(node_list_fc)
        else:
            print("input erro: force_opt_dir don't defined")

        return fcdof

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
        
    def setNodes(nodelist, coord, regions):
        tol = 1e-10
        # ----- SEEKERS WITH LOC -----
        if nodelist[0] == "lengthx":
            coord_0 = float(nodelist[2])
            coord_1 = float(nodelist[3])
            nodes = coord[np.where((coord[:, 1] >= coord_0)&(coord[:, 1] <= coord_1)), 0,][0]
            # nodesapply.append(nodes)
            dir_fc = "x"
        
        elif nodelist[0] == "lengthy":
            coord_0 = float(nodelist[2])
            coord_1 = float(nodelist[3])
            nodes = coord[np.where((coord[:, 2] >= coord_0)&(coord[:, 2] <= coord_1)), 0,][0]
            # nodesapply.append(nodes)
            dir_fc = "y"
        
        elif nodelist[0] == "lengthz":
            coord_0 = float(nodelist[2])
            coord_1 = float(nodelist[3])
            nodes = coord[np.where((coord[:, 3] >= coord_0)&(coord[:, 3] <= coord_1)), 0,][0]
            # nodesapply.append(nodes)
            dir_fc = "z"
        
        elif nodelist[0] == "edgex":
            edge_coordX = nodelist[1].astype(float) #float(nodelist[1])
            if int(nodelist[2]) == 999:
                dir_fc = "x_y"
                coord_fc = (coord[np.where(coord[:, 3] == float(nodelist[3])),:,])[0]
            elif int(nodelist[3]) == 999:
                dir_fc = "x_z"
                coord_fc = (coord[np.where(coord[:, 2] == float(nodelist[2])),:,])[0]
            nodes = search_edgex(edge_coordX, coord_fc, tol)
            # nodesapply.append(nodes)
        
        elif nodelist[0] == "edgey":
            edge_coordY = float(nodelist[2])
            if float(nodelist[1]) == 999:
                dir_fc = "y_x"
                coord_fc = (coord[np.where(coord[:, 3] == float(nodelist[3])),:,])[0]
            elif float(nodelist[3]) == 999:
                dir_fc = "y_z"
                coord_fc = (coord[np.where(coord[:, 1] == float(nodelist[1])),:,])[0]
            nodes = search_edgey(edge_coordY, coord_fc, tol)
            # nodesapply.append(nodes)
        
        elif nodelist[0] == "edgez":
            edge_coordZ = float(nodelist[3])
            if float(nodelist[1]) == 999:
                dir_fc = "z_x"
                coord_fc = (coord[np.where(coord[:, 2] == float(nodelist[2])),:,])[0]
            elif float(nodelist[2]) == 999:
                dir_fc = "z_x"
                coord_fc = (coord[np.where(coord[:, 1] == float(nodelist[1])), :,])[0]
            nodes = search_edgez(edge_coordZ, coord_fc, tol)
            # nodesapply.append(nodes)
        
        elif nodelist[0] == "surfxy":
            orthg_coordZ = float(nodelist[3])
            nodes = search_surfxy(orthg_coordZ, coord, tol)
            # nodesapply.append(nodes)
            dir_fc = "z"
        
        elif nodelist[0] == "surfyz":
            orthg_coordX = float(nodelist[1])
            nodes = search_surfyz(orthg_coordX, coord, tol)
            # nodesapply.append(nodes)
            dir_fc = "x"
        
        elif nodelist[0] == "surfzx":
            orthg_coordY = float(nodelist[2])
            nodes = search_surfzx(orthg_coordY, coord, tol)
            # nodesapply.append(nodes)
            dir_fc = "y"
        
        elif nodelist[0] == "node":
            node_coordX = float(nodelist[1])
            node_coordY = float(nodelist[2])
            node_coordZ = float(nodelist[3])
            nodes = search_nodexyz(node_coordX, node_coordY, node_coordZ, coord, tol)
            # nodesapply.append(nodes)
            dir_fc = "x"
            
        # ----- SEEKERS WITH TAG -----
        elif nodelist[0] == "point":
            nodes = regions[0][1][int(nodelist[-1]) - 1][1][:]
            dir_fc = "x"
            # nodesapply.append(nodes)
        
        elif nodelist[0] == "line":
            nodes = regions[1][1][int(nodelist[-1]) - 1][1][:]
            dir_fc = "x"
            # nodesapply.append(nodes)
        
        elif nodelist[0] == "plane":
            nodes = regions[2][1][int(nodelist[-1]) - 1][1][:]
            dir_fc = "x"
            # nodesapply.append(nodes)
        
        else:
            print("input erro: force_opt don't defined")
            
        return nodes, dir_fc
    

