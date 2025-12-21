import os

import numpy as np

from myfempy.core.mesh.mesh import Mesh
from myfempy.core.utilities import nodes_from_regions
from myfempy.io.iogmsh import meshid2gmshid


class MeshGmsh(Mesh):
    """Mesh GMSH Class <ConcreteClassService>"""

    def getElementConection(set_mesh):
        conec, __ = MeshGmsh.__setmesh_from_gmsh(set_mesh)
        return conec

    def getNodesCoord(set_mesh):
        __, nodes = MeshGmsh.__setmesh_from_gmsh(set_mesh)
        return nodes

    def getElementList(conec, meshset, modeldata):
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
                        modeldata["MATERIAL"]["PROPMAT"][int(conec[ee][3]) - 1]["NAME"],
                        modeldata["GEOMETRY"]["PROPGEO"][int(conec[ee][3]) - 1]["NAME"],
                        (np.array(conec[ee][5:])).astype(int).tolist(),
                    ]
                )

        elemlist = elemlist[1::][::]

        return elemlist

    def getRegionsList(conec):
        return MeshGmsh.__setregions(conec)

    def __setmesh_from_gmsh(set_mesh):
        if "meshimport" in set_mesh.keys():
            conec, nodes = MeshGmsh.__convert_from_msh2(set_mesh["user_path"] + "/" + set_mesh["meshimport"]['object'])
        else:
            conec, nodes = MeshGmsh.__convert_from_msh2(
                set_mesh["user_path"] + "/" + set_mesh["filename"]
            )
        return conec, nodes

    # def __convert_from_msh1(filename):
    #     file_imp = filename + ".msh1"
    #     with open(file_imp, "r") as file_object:
    #         file_object.readline()
    #         NNOD = int(file_object.readline())
    #         nodelist = [[None] * 4]
    #         for ii in range(0, NNOD):
    #             line = file_object.readline()
    #             lineaux = line.split()
    #             contstr = lineaux[0:4]
    #             nodelist.append(
    #                 [
    #                     float(contstr[0]),
    #                     float(contstr[1]),
    #                     float(contstr[2]),
    #                     float(contstr[3]),
    #                 ]
    #             )
    #         nodelist = nodelist[1::][::]
    #         file_object.readline()
    #         file_object.readline()
    #         NELM = int(file_object.readline())
    #         conec_elm = []
    #         for kk in range(0, NELM):
    #             line = file_object.readline()
    #             lineaux = line.split()
    #             conec_elm.append(list(map(float, lineaux[:])))

    #     return conec_elm, nodelist    

    def __convert_from_msh2(filename):
        file_imp = filename + ".msh2"
        with open(file_imp, "r") as file_object:
            file_object.readline()
            file_object.readline()
            file_object.readline()
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
            elif (int(conec[ee][1]) == 1) or (int(conec[ee][1]) == 8):
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
        edgelist = {name.capitalize(): {key: [] for key in keys} for name in list_edge}
        list_surf = (np.arange(1, contsfr + 1).astype(str)).tolist()
        surflist = {name.capitalize(): {key: [] for key in keys} for name in list_surf}

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
