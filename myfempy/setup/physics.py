from __future__ import annotations

import numpy as np

from myfempy.core.utilities import (get_elemen_from_nodelist,
                                    get_nodes_from_list)


class SetPhysics:
    """Set Physics Class <ClassOrder>"""

    def __init__(self, Model, Loads, BoundCond):
        self.loads = Loads
        self.boundcond = BoundCond
        self.model = Model
        # self.domain = Domain

    def getForceList(self, physicdata):
        return SetPhysics.setForceList(self, physicdata)

    def setForceList(self, physicdata):
        forcelist = SetPhysics.__forcelist(physicdata["LOAD"])
        self.forcelist = forcelist
        return forcelist

    def getBoundCondList(self, physicdata):
        return SetPhysics.setBoundCondList(self, physicdata)

    def setBoundCondList(self, physicdata):
        boundlist = SetPhysics.__boundcondlist(physicdata["BOUNDCOND"])
        self.boundcondlist = boundlist
        return boundlist

    def getLoadApply(self, physicdata, modelinfo):
        flist = SetPhysics.setForceList(self, physicdata)
        forcenodeaply = np.zeros((1, 4))
        for step in range(int(len(np.unique([flist[:, 8]])))):
            forcelist = flist[np.where(flist[:, 8] == str(step + 1))[0], :]
            fapp = self.loads.getLoadApply(self.model, modelinfo, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def getBoundCondApply(self, physicdata, modelinfo):
        bclist = SetPhysics.setBoundCondList(self, physicdata)
        boncdnodeaply = np.zeros((1, 4))
        for step in range(int(len(np.unique([bclist[:, 8]])))):
            listapply = bclist[np.where(bclist[:, 8] == str(step))[0], :]
            bcapp = self.boundcond.getBCApply(modelinfo, listapply)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
        boncdnodeaply = boncdnodeaply[1::][::]
        return boncdnodeaply

    def getNodeList(self, domain_nodelist, coord, regions):
        return get_nodes_from_list(domain_nodelist, coord, regions)

    def getElementList(self, inci, nodelist):
        return get_elemen_from_nodelist(inci, nodelist)

    # -----------------------------------------------
    # privates methods
    def __forcelist(forcelist):
        """force set"""
        nforc = len(forcelist)
        flist = np.zeros((1, 9))
        for fl in range(nforc):
            fap = forcelist[fl]
            for fs in range(len(fap["VAL"])):
                if "LOC" in fap.keys():
                    linearray = np.array(
                        [
                            fap["TYPE"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            fap["LOC"]["x"],
                            fap["LOC"]["y"],
                            fap["LOC"]["z"],
                            0,
                            int(fs + 1),
                        ]
                    )
                # else:
                elif "TAG" in fap.keys():
                    linearray = np.array(
                        [
                            fap["TYPE"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            0.0,
                            0.0,
                            0.0,
                            fap["TAG"],
                            int(fs + 1),
                        ]
                    )
                else:
                    linearray = np.array(
                        [
                            fap["TYPE"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            0.0,
                            0.0,
                            0.0,
                            0,
                            int(fs + 1),
                        ]
                    )
                flist = np.append(flist, [linearray], axis=0)
        flist = flist[1::][::]
        return flist

    def __boundcondlist(boundcondlist):
        """boundary conditions set"""
        nbound = len(boundcondlist)
        blist = np.zeros((1, 9))
        for bl in range(nbound):
            bap = boundcondlist[bl]
            if "VAL" in bap.keys():
                for bs in range(len(bap["VAL"])):
                    if "LOC" in bap.keys():
                        linearray = np.array(
                            [
                                bap["TYPE"],
                                bap["DOF"],
                                bap["DIR"],
                                bap["LOC"]["x"],
                                bap["LOC"]["y"],
                                bap["LOC"]["z"],
                                0,
                                bap["VAL"][bs],
                                int(bs + 1),
                            ]
                        )
                    elif "TAG" in bap.keys():
                        linearray = np.array(
                            [
                                bap["TYPE"],
                                bap["DOF"],
                                bap["DIR"],
                                0.0,
                                0.0,
                                0.0,
                                bap["TAG"],
                                bap["VAL"][bs],
                                int(bs + 1),
                            ]
                        )
                    else:
                        linearray = np.array(
                            [
                                bap["TYPE"],
                                bap["DOF"],
                                bap["DIR"],
                                0.0,
                                0.0,
                                0.0,
                                0,
                                bap["VAL"][bs],
                                int(bs + 1),
                            ]
                        )
                    blist = np.append(blist, [linearray], axis=0)
            else:
                if "LOC" in bap.keys():
                    linearray = np.array(
                        [
                            bap["TYPE"],
                            bap["DOF"],
                            bap["DIR"],
                            bap["LOC"]["x"],
                            bap["LOC"]["y"],
                            bap["LOC"]["z"],
                            0,
                            0.0,
                            0,
                        ]
                    )
                elif "TAG" in bap.keys():
                    linearray = np.array(
                        [
                            bap["TYPE"],
                            bap["DOF"],
                            bap["DIR"],
                            0.0,
                            0.0,
                            0.0,
                            bap["TAG"],
                            0.0,
                            0,
                        ]
                    )
                else:
                    linearray = np.array(
                        [
                            bap["TYPE"],
                            bap["DOF"],
                            bap["DIR"],
                            0.0,
                            0.0,
                            0.0,
                            0,
                            0.0,
                            0,
                        ]
                    )
                blist = np.append(blist, [linearray], axis=0)
        blist = blist[1::][::]
        return blist
