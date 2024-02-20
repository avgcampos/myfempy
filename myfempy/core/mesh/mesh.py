import os
from abc import ABC, abstractmethod

import numpy as np
from myfempy.io.iogmsh import get_gmsh_geo, get_gmsh_msh, meshid2gmshid


def getMesh(set_mesh):
        if set_mesh['TYPE'] == 'add':
            return MeshADD
    
        elif set_mesh['TYPE'] == 'legacy':
            if set_mesh['shape'] == 'quad4':
                from myfempy.core.mesh.legacyquad4 import LegacyQuad4
                return LegacyQuad4
            elif set_mesh['shape'] == 'tria3':
                from myfempy.core.mesh.legacytria3 import LegacyTria3
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
            from myfempy.core.mesh.gmsh import MeshGmsh
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
    
