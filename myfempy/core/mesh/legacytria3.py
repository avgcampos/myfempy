
import numpy as np

from myfempy.core.mesh.mesh import Mesh


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
                        
        conec = np.zeros((nel, 4), dtype=np.int64)
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
                             
        coord = np.zeros((nos, 4), dtype=np.float64)
        coord[0, 0] = 1
        for i in range(2, nos + 1):
            linha = int(np.ceil(i / (nelx + 1))) - 1
            coord[i - 1, 0] = i
            coord[i - 1, 1] = ((i - 1) - linha * (nelx + 1)) * (lx / nelx)
            coord[i - 1, 2] = linha * (ly / nely)
        return coord

    def getElementList(conec, meshset, modeldata):
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