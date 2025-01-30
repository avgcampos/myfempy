from __future__ import annotations

import numpy as np

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import get_nodes_from_list


class BoundCondStruct(Structural):
    """Structural Load Class <ConcreteClassService>"""

    def getBCApply(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 4))
        if bclist["TYPE"] == "fixed":
            bcapp = BoundCondStruct.__BCFixed(modelinfo, bclist)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
        elif bclist["TYPE"] == "displ":
            bcapp = BoundCondStruct.__BCDispl(modelinfo, bclist)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
        elif bclist["TYPE"] == "csymm":
            bcapp = BoundCondStruct.__BCCS(modelinfo, bclist)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
        else:
            pass
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    def __BCFixed(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 4))

        nodelist = [
            bclist["DIR"],
            bclist["LOCX"],
            bclist["LOCY"],
            bclist["LOCZ"],
            bclist["TAG"],
        ]
        node_list_bc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )

        if bclist["DOF"] == "full":
            bcdof = 0
        else:
            bcdof = modelinfo["dofs"]["d"][bclist["DOF"]]

        for j in range(len(node_list_bc)):
            bcapp = np.array([[int(node_list_bc[j]), bcdof, 0.0, int(bclist["STEP"])]])
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    def __BCDispl(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 4))

        nodelist = [
            bclist["DIR"],
            bclist["LOCX"],
            bclist["LOCY"],
            bclist["LOCZ"],
            bclist["TAG"],
        ]
        node_list_bc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )

        bcdof = modelinfo["dofs"]["d"][bclist["DOF"]]

        for j in range(len(node_list_bc)):
            bcapp = np.array(
                [
                    [
                        int(node_list_bc[j]),
                        bcdof,
                        float(bclist["VAL"]),
                        int(bclist["STEP"]),
                    ]
                ]
            )
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    def __BCCS(modelinfo, bclist):
        boncdnodeaply = np.zeros((1, 4))

        nodelist = [
            bclist["DIR"],
            bclist["LOCX"],
            bclist["LOCY"],
            bclist["LOCZ"],
            bclist["TAG"],
        ]
        node_list_bc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )

        if bclist["DOF"] == "left":
            bcdof = 11
        elif bclist["DOF"] == "right":
            bcdof = 12
        else:
            bcdof = 0

        for j in range(len(node_list_bc)):
            bcapp = np.array([[int(node_list_bc[j]), bcdof, 0.0, int(bclist["STEP"])]])
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)

        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply
