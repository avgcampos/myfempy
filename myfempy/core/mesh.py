import os
from abc import ABC, abstractmethod

import numpy as np

from myfempy.core.utils import nodes_from_regions
from myfempy.io.iogmsh import get_gmsh_geo, get_gmsh_msh, meshid2gmshid


def getMesh(set_mesh):
        if set_mesh['TYPE'] == 'add':
            return MeshADD
    
        elif set_mesh['TYPE'] == 'legacy':
            if set_mesh['shape'] == 'quad4':
                return LegacyQuad4
            elif set_mesh['shape'] == 'tria3':
                return LegacyTria3
            else:
                pass
        
        elif set_mesh['TYPE'] == 'gmsh':
            if "meshimport" in set_mesh.keys():
                pass
            else:
                filename = set_mesh["filename"]
                get_gmsh_geo(filename, set_mesh)
                get_gmsh_msh(filename, set_mesh)
            return MeshGmsh
        else:
            pass

class Mesh(ABC):
    '''Element API Class <ClassService>'''
    
    @abstractmethod
    def getElementConection():
        pass

    @abstractmethod
    def getNodesCoord():
        pass

    @abstractmethod
    def getElementList():
        pass
        
    def getElementListMeshADD(conec, meshset, modeldata):
        elemlist = [[None] * 3]
        for ee in range(len(conec)):
            elemlist.append(
                [
                    int(conec[ee, 0]),
                    meshset,
                    modeldata["MATERIAL"]["PROPMAT"][conec[ee, 1]-1]["NAME"],
                    modeldata["GEOMETRY"]["PROPGEO"][conec[ee, 2]-1]["NAME"],
                    conec[ee, 3:].astype(int).tolist(),
                ]
            )
        elemlist = elemlist[1::][::]
        return elemlist
        
    def getElementListMeshLEGACY(conec, meshset, modeldata):
        elemlist = [[None] * 3]
        for ee in range(len(conec)):
            elemlist.append(
                [
                    int(conec[ee, 0]),
                    meshset,
                    modeldata["MATERIAL"]["PROPMAT"][0]["NAME"],
                    modeldata["GEOMETRY"]["PROPGEO"][0]["NAME"],
                    conec[ee, 1:].astype(int).tolist(),
                ]
            )
        elemlist = elemlist[1::][::]
        return elemlist
    
    def getElementListMeshGMSH(conec, meshset, modeldata):

        meshid = meshid2gmshid(str(meshset))
        elemlist = [[None] * 5]
        contelm = 0
        for ee in range(len(conec)):
            if int(conec[ee][1]) == meshid:
                contelm += 1
                elemlist.append(
                    [
                        contelm,
                        meshset,
                        modeldata["MATERIAL"]["PROPMAT"][int(conec[ee][3])-1]["NAME"],
                        modeldata["GEOMETRY"]["PROPGEO"][int(conec[ee][3])-1]["NAME"],
                        (np.array(conec[ee][5:])).astype(int).tolist(),
                    ]
                )
        
        elemlist = elemlist[1::][::]
        
        return elemlist
         


class MeshADD(Mesh):
    '''Mesh Add Class <ConcreteClassService>'''
    
    def getNodesCoord(set_mesh):
        return np.asarray(set_mesh['COORD'])
    
    def getElementConection(set_mesh):
        return np.asarray(set_mesh['INCI'])
    
    def getElementList(conec, elemtype, modeldata):
        return Mesh.getElementListMeshADD(conec, elemtype, modeldata)

class LegacyQuad4(Mesh):
    '''Mesh Quad Class <ConcreteClassService>'''
        
    def getElementConection(set_mesh):
        """get a quadrangular 4 nodes mesh
        
        (l)---------------(k)
         |                 |
         |                 |
         |      {1}        |
         |                 |
         |                 |
        (i)---------------(j)
            
        """
        
        nelx = set_mesh["NX"]
        nely = set_mesh["NY"]
        
        nnx = nelx + 1        
        nel = nelx * nely
                        
        conec = np.zeros((nel, 5))
        for i in range(1, nely + 1, 1):
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 0] = np.arange((i - 1) * nelx + 1, i * nelx + 1, 1).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 1] = np.arange((i - 1) * nnx + 1, i * nnx, 1).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 2] = np.arange((i - 1) * nnx + 2, i * nnx + 1, 1).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 3] = np.arange(i * nnx + 2, (i + 1) * nnx + 1, 1).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 4] = np.arange(i * nnx + 1, (i + 1) * nnx, 1).tolist()
        
        return conec

    def getNodesCoord(set_mesh):
        nelx = set_mesh["NX"]
        nely = set_mesh["NY"]
        lx = set_mesh["LX"]
        ly = set_mesh["LY"]
        
        nnx = nelx + 1
        nny = nely + 1
        nos = nnx * nny
                
        dx = np.arange(0, lx + lx / nelx, lx / nelx)
        dy = np.arange(0, ly + ly / nely, ly / nely)                
        xv, yv = np.meshgrid(dx, dy)
                
        coord = np.zeros((nos, 4))
        for i in range(1, nny + 1, 1):
            coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 0] = np.arange((i - 1) * nnx + 1, i * nnx + 1, 1).tolist()
            coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 1] = xv[i - 1, :]
            coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 2] = yv[i - 1, :]
        
        return coord

    def getElementList(conec, elemtype, modeldata):
        return Mesh.getElementListMeshLEGACY(conec, elemtype, modeldata)
    

class LegacyTria3(Mesh):
    '''Mesh Quad Class <ConcreteClassService>'''
        
    def getElementConection(set_mesh):
        """get a triagular 3 nodes mesh
  
        (k) 
        |\ 
        | \ 
        |  \
        |   \
        |    \
        |     \
        |  {1} \ 
        |       \
        |        \
        |         \
        |          \ 
        (i)--------(j)
        
        """
        
        nelx = set_mesh["NX"]
        nely = set_mesh["NY"]
        nel = nelx * nely * 2
                        
        conec = np.zeros((nel, 4))
        for i in range(1, nel, 2):
            linha = int(np.ceil(i / (2 * nelx)))
            y = 2 * linha - 1
            n1 = (i + y) / 2
            n2 = n1 + 1
            n3 = n2 + nelx
            conec[i - 1, 0] = i
            conec[i - 1, 1] = n1
            conec[i - 1, 2] = n2
            conec[i - 1, 3] = n3
        
        for i in range(2, nel + 1, 2):
            linha = int(np.ceil(i / (2 * nelx)))
            y = 2 * linha
            n1 = (i + y) / 2
            n2 = n1 + nelx
            n3 = n2 + 1
            conec[i - 1, 0] = i
            conec[i - 1, 1] = n1
            conec[i - 1, 2] = n2
            conec[i - 1, 3] = n3

        return conec

    def getNodesCoord(set_mesh):
        nelx = set_mesh["NX"]
        nely = set_mesh["NY"]
        lx = set_mesh["LX"]
        ly = set_mesh["LY"]
        
        nnx = nelx + 1
        nny = nely + 1
        nos = nnx * nny
                             
        coord = np.zeros((nos, 4))
        coord[0, 0] = 1
        for i in range(2, nos + 1):
            linha = int(np.ceil(i / (nelx + 1))) - 1
            coord[i - 1, 0] = i
            coord[i - 1, 1] = ((i - 1) - linha * (nelx + 1)) * (lx / nelx)
            coord[i - 1, 2] = linha * (ly / nely)
        return coord

    def getElementList(conec, elemtype, modeldata):
        return Mesh.getElementListMeshLEGACY(conec, elemtype, modeldata)



class MeshGmsh(Mesh):
    '''Mesh GMSH Class <ConcreteClassService>'''

    def getElementConection(set_mesh):
        conec, nodes = MeshGmsh.__setmesh_from_gmsh(set_mesh)
        return conec
    
    def getNodesCoord(set_mesh):
        conec, nodes = MeshGmsh.__setmesh_from_gmsh(set_mesh)
        return nodes
    
    def getElementList(conec, meshset, modeldata):
        return Mesh.getElementListMeshGMSH(conec, meshset, modeldata)
    
    def getRegionsList(conec):
        return MeshGmsh.__setregions(conec)
    
    def __setmesh_from_gmsh(set_mesh):
        
        if "meshimport" in set_mesh.keys():
            conec, nodes = MeshGmsh.__convert_from_msh1(
                os.getcwd() + "/" + set_mesh["meshimport"]
            )
        else:         
            conec, nodes = MeshGmsh.__convert_from_msh1(
                os.getcwd() + "/" + set_mesh["filename"]
            )
        return conec, nodes            

    def __convert_from_msh1(filename):
        file_imp = filename + ".msh1"
        with open(file_imp, "r") as file_object:
            file_object.readline()
            NNOD = int(file_object.readline())
            nodelist = [[None] * 4]
            for ii in range(0, NNOD):
                line = file_object.readline()
                lineaux = line.split()
                contstr = lineaux[0:4]
                nodelist.append(
                    [
                        float(contstr[0]),
                        float(contstr[1]),
                        float(contstr[2]),
                        float(contstr[3]),
                    ]
                )
            nodelist = nodelist[1::][::]
            file_object.readline()
            file_object.readline()
            NELM = int(file_object.readline())
            conec_elm = []
            for kk in range(0, NELM):
                line = file_object.readline()
                lineaux = line.split()
                conec_elm.append(list(map(float, lineaux[:])))
        
        return conec_elm, nodelist


    def __setregions(conec):

        reglist = dict()
        reglist["reglist"] = []
        contptn = 0
        contedg = 0
        contsfr = 0
        contvlr = 0
        vlr0 = 0
        sfr0 = 0
        edg0 = 0
        ptn0 = 0
        for ee in range(len(conec)):
            if int(conec[ee][1]) == 15:
                ptn1 = int(conec[ee][3])
                if ptn1 != ptn0:
                    contptn += 1
                reglist["reglist"].append(
                    {
                        "type": "point",
                        "list": str(contptn),
                        "nodes": (np.array(conec[ee][5:])).astype(int).tolist(),
                    }
                )
                ptn0 = ptn1
            elif int(conec[ee][1]) == 1:
                edg1 = int(conec[ee][3])
                if edg1 != edg0:
                    contedg += 1
                reglist["reglist"].append(
                    {
                        "type": "line",
                        "list": str(contedg),
                        "nodes": (np.array(conec[ee][5:])).astype(int).tolist(),
                    }
                )
                edg0 = edg1
            elif (int(conec[ee][1]) == 2) or (int(conec[ee][1]) == 3):
                sfr1 = int(conec[ee][3])
                if sfr1 != sfr0:
                    contsfr += 1
                reglist["reglist"].append(
                    {
                        "type": "plane",
                        "list": str(contsfr),
                        "nodes": (np.array(conec[ee][5:])).astype(int).tolist(),
                    }
                )
                sfr0 = sfr1
        keys = ["nodes"]
        list_point = (np.arange(1, contptn + 1).astype(str)).tolist()
        pointlist = {
            name.capitalize(): {key: [] for key in keys} for name in list_point
        }
        list_edge = (np.arange(1, contedg + 1).astype(str)).tolist()
        edgelist = {
            name.capitalize(): {key: [] for key in keys} for name in list_edge
        }
        list_surf = (np.arange(1, contsfr + 1).astype(str)).tolist()
        surflist = {
            name.capitalize(): {key: [] for key in keys} for name in list_surf
        }
        
        for rr in range(len(reglist["reglist"])):
            reg = reglist["reglist"][rr]
            if reg["type"] == "point":
                pointlist[reg["list"]]["nodes"].extend(reg["nodes"])
            elif reg["type"] == "line":
                edgelist[reg["list"]]["nodes"].extend(reg["nodes"])
            elif reg["type"] == "plane":
                surflist[reg["list"]]["nodes"].extend(reg["nodes"])
            else:
                pass
        
        regionlist = dict()
        regionlist["point"] = pointlist
        regionlist["line"] = edgelist
        regionlist["plane"] = surflist
        
        regions = nodes_from_regions(regionlist)
        
        return regions